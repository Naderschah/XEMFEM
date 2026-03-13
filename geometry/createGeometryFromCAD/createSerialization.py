#!/usr/bin/env python3
import argparse
import os
import json
import math
import glob
import re
import traceback
from datetime import datetime
from collections import Counter

import ezdxf
from ezdxf import edgesmith, edgeminer
from ezdxf.math import Vec3, OCS

# -------------------------
# Config
# -------------------------
ROOT_DIR = "/work/geometry/createGeometryFromCAD/DXF_slices_parts"
DXF_GLOB = "**/*.dxf"
GAP_TOL  = 1e-1
MAX_FORCE_MERGE_DIST = 1e-1
ENABLE_MANUAL_LINE_CLOSURE = False

# SHAPER Sketch.addArc(..., direction) convention:
# From your observation, the "inverted CCW" worked for most => True probably means "CW" (i.e. not-CCW).
TRUE_MEANS_CCW = False

FAILED_FILES_LOG   = os.path.join(os.path.abspath(ROOT_DIR), "failed_files.log")
FAILED_OBJECTS_LOG = os.path.join(os.path.abspath(ROOT_DIR), "failed_objects.log")


def _parse_args():
    p = argparse.ArgumentParser(
        description="Serialize DXF components into JSON sketch instruction files."
    )
    p.add_argument(
        "--root-dir",
        default=ROOT_DIR,
        help="Root directory searched by --dxf-glob (default: %(default)s).",
    )
    p.add_argument(
        "--dxf-glob",
        default=DXF_GLOB,
        help="Glob relative to --root-dir used to find DXF files (default: %(default)s).",
    )
    p.add_argument(
        "--dxf-file",
        action="append",
        default=[],
        help="Specific DXF file to process; can be passed multiple times. If set, glob search is skipped.",
    )
    p.add_argument(
        "--path-regex",
        default=None,
        help="Regex filter applied to absolute DXF path before processing.",
    )
    p.add_argument(
        "--gap-tol",
        type=float,
        default=GAP_TOL,
        help="Endpoint gap tolerance used in contour extraction (default: %(default)s).",
    )
    return p.parse_args()


# -------------------------
# Logging
# -------------------------
def _reset_log(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        ts = datetime.now().isoformat(timespec="seconds")
        f.write(f"{ts} [RUN_START]\n")

def _log_line(path, msg):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ts = datetime.now().isoformat(timespec="seconds")
    with open(path, "a", encoding="utf-8") as f:
        f.write(f"{ts} {msg}\n")

def _log_exception(path, header, exc):
    _log_line(
        path,
        f"{header} exc_type={type(exc).__name__} exc_msg={str(exc)}"
    )
    tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)).rstrip()
    _log_line(path, f"{header}\n{tb}\n---")

def _infer_no_components_reason(entity_types):
    types = set(entity_types.keys())
    if not types:
        return "modelspace had no supported entities"
    if types.issubset({"LINE"}):
        return "only LINE entities; no closed contour found within GAP_TOL"
    if types.issubset({"ARC", "LINE"}):
        return "ARC/LINE entities did not form a closed contour"
    if types.issubset({"SPLINE"}):
        return "SPLINE entities were open/unlinked or unsupported for closed extraction"
    if "SPLINE" in types:
        return "mixed entities including SPLINE did not resolve into closed/chainable contour"
    return "entities did not resolve into closed/chainable contour"


# -------------------------
# Small math helpers
# -------------------------
def _vxy(p, nd=12):
    p = Vec3(p)
    return (round(p.x, nd), round(p.y, nd))

def _cross2(ax, ay, bx, by):
    return ax * by - ay * bx

def _traversal_is_ccw_xy(x1, y1, x2, y2, cx, cy):
    return _cross2(x1 - cx, y1 - cy, x2 - cx, y2 - cy) > 0.0

def _shaper_direction_from_traversal_ccw(traversal_is_ccw: bool) -> bool:
    return traversal_is_ccw if TRUE_MEANS_CCW else (not traversal_is_ccw)

def _norm_deg(a):
    a = float(a) % 360.0
    return a + 360.0 if a < 0 else a

def _ang_diff(a, b):
    d = (a - b) % 360.0
    return min(d, 360.0 - d)


# -------------------------
# Robust arc direction using ARC start/end angles (handles major arcs)
# -------------------------
def arc_direction_flag_for_shaper(arc_entity, x1, y1, x2, y2, tol_deg=1e-3):
    ext = Vec3(getattr(arc_entity.dxf, "extrusion", (0, 0, 1)))
    ocs = OCS(ext)

    c_wcs = Vec3(arc_entity.dxf.center)
    c = ocs.from_wcs(c_wcs)

    p1 = ocs.from_wcs(Vec3(x1, y1, 0.0))
    p2 = ocs.from_wcs(Vec3(x2, y2, 0.0))

    a1 = _norm_deg(math.degrees(math.atan2(p1.y - c.y, p1.x - c.x)))
    a2 = _norm_deg(math.degrees(math.atan2(p2.y - c.y, p2.x - c.x)))

    sa = _norm_deg(arc_entity.dxf.start_angle)
    ea = _norm_deg(arc_entity.dxf.end_angle)

    m_fwd = _ang_diff(a1, sa) + _ang_diff(a2, ea)
    m_rev = _ang_diff(a1, ea) + _ang_diff(a2, sa)

    if abs(m_fwd - m_rev) < tol_deg:
        traversal_is_ccw = _traversal_is_ccw_xy(x1, y1, x2, y2, float(c_wcs.x), float(c_wcs.y))
    else:
        traversal_is_ccw = (m_fwd < m_rev)  # DXF: CCW start->end

    return _shaper_direction_from_traversal_ccw(traversal_is_ccw)


# -------------------------
# DXF spline utilities
# -------------------------
def _knot_multiplicities(full_knots, tol=1e-12):
    """
    DXF stores the full knot vector; SHAPER wants unique knots + multiplicities.
    """
    if not full_knots:
        return [], []
    uniq = [float(full_knots[0])]
    mults = [1]
    for k in full_knots[1:]:
        k = float(k)
        if abs(k - uniq[-1]) <= tol:
            mults[-1] += 1
        else:
            uniq.append(k)
            mults.append(1)
    return uniq, mults

def _dist2_xy(a, b):
    return (a.x - b.x) ** 2 + (a.y - b.y) ** 2

def _reverse_bspline_payload_inplace(payload_dict):
    """
    Reverse a non-periodic (clamped/open) spline definition.
    Do NOT use for periodic splines.
    """
    poles = payload_dict.get("poles", [])
    weights = payload_dict.get("weights", [])
    knots = payload_dict.get("knots", [])
    mults = payload_dict.get("mults", payload_dict.get("multiplicities", []))

    payload_dict["poles"] = list(reversed(poles))
    if weights:
        payload_dict["weights"] = list(reversed(weights))

    if knots:
        kmin = float(knots[0])
        kmax = float(knots[-1])
        payload_dict["knots"] = [float((kmin + kmax) - k) for k in reversed(knots)]
    if mults:
        payload_dict["mults"] = list(reversed(mults))


def _spline_flags(e) -> int:
    return int(getattr(e.dxf, "flags", 0))

def _spline_is_closed(e) -> bool:
    # DXF: bit 1 = closed
    return bool(_spline_flags(e) & 1)

def _spline_is_periodic(e) -> bool:
    # DXF: bit 2 = periodic
    return bool(_spline_flags(e) & 2)


# -------------------------
# Ellipse helpers (for ELLIPSE entity and full ellipse export)
# -------------------------
def _ellipse_params_from_dxf_ellipse(e):
    c = Vec3(e.dxf.center)
    maj = Vec3(e.dxf.major_axis)
    a = maj.magnitude
    b = a * float(e.dxf.ratio)
    phi_deg = math.degrees(math.atan2(maj.y, maj.x))
    return (float(c.x), float(c.y), float(a), float(b), float(phi_deg))

def _is_full_ellipse(e, tol=1e-9):
    sp = float(getattr(e.dxf, "start_param", 0.0))
    ep = float(getattr(e.dxf, "end_param", 0.0))
    two_pi = 2.0 * math.pi
    span = (ep - sp) % two_pi
    return abs(span - two_pi) < tol or abs(span) < tol


# -------------------------
# Closed primitive exporters (because EdgeSmith ignores closed shapes)
# -------------------------
def _circle_to_component(e):
    c = Vec3(e.dxf.center)
    r = float(e.dxf.radius)
    cx, cy = float(c.x), float(c.y)

    # Split into two arcs to satisfy vertex-loop format.
    x0, y0 = cx + r, cy
    x1, y1 = cx - r, cy

    direction = _shaper_direction_from_traversal_ccw(True)  # CCW traversal
    return {"pts": [
        ["arc", float(x0), float(y0), cx, cy, direction],
        ["arc", float(x1), float(y1), cx, cy, direction],
    ]}

def _bulge_to_instr(x1, y1, x2, y2, bulge):
    if abs(bulge) < 1e-15:
        return ["line", x1, y1]

    theta = 4.0 * math.atan(bulge)
    chord = math.hypot(x2 - x1, y2 - y1)
    if chord < 1e-15:
        return ["line", x1, y1]

    r = chord / (2.0 * math.sin(abs(theta) / 2.0))
    mx, my = (x1 + x2) * 0.5, (y1 + y2) * 0.5
    h = math.sqrt(max(r * r - (chord * 0.5) ** 2, 0.0))

    dx, dy = (x2 - x1) / chord, (y2 - y1) / chord
    nx, ny = -dy, dx  # left normal
    if bulge < 0.0:
        nx, ny = -nx, -ny

    cx, cy = mx + nx * h, my + ny * h

    traversal_is_ccw = bulge > 0.0
    direction = _shaper_direction_from_traversal_ccw(traversal_is_ccw)
    return ["arc", x1, y1, float(cx), float(cy), direction]

def _lwpolyline_closed_to_component(e):
    if not bool(getattr(e, "closed", False)):
        return None
    pts = list(e.get_points("xyb"))  # (x,y,bulge)
    if not pts:
        return None

    out = []
    n = len(pts)
    for i in range(n):
        x1, y1, b = pts[i]
        x2, y2, _ = pts[(i + 1) % n]
        out.append(_bulge_to_instr(float(x1), float(y1), float(x2), float(y2), float(b)))
    return {"pts": out}

def _polyline_closed_to_component(e):
    if not bool(getattr(e, "is_closed", False)):
        return None
    verts = list(getattr(e, "vertices", []))
    if not verts:
        return None

    out = []
    pts = []
    for v in verts:
        p = Vec3(v.dxf.location)
        b = float(getattr(v.dxf, "bulge", 0.0))
        pts.append((float(p.x), float(p.y), b))

    n = len(pts)
    for i in range(n):
        x1, y1, b = pts[i]
        x2, y2, _ = pts[(i + 1) % n]
        out.append(_bulge_to_instr(x1, y1, x2, y2, b))
    return {"pts": out}

def _full_ellipse_to_component(e):
    if not _is_full_ellipse(e):
        return None
    cx, cy, a, b, phi_deg = _ellipse_params_from_dxf_ellipse(e)
    phi = math.radians(phi_deg)
    x0, y0 = cx + a * math.cos(phi), cy + a * math.sin(phi)
    x1, y1 = cx - a * math.cos(phi), cy - a * math.sin(phi)

    inversed = False
    return {"pts": [
        ["ellipse", float(x0), float(y0), cx, cy, a, b, phi_deg, inversed],
        ["ellipse", float(x1), float(y1), cx, cy, a, b, phi_deg, inversed],
    ]}

def _closed_or_periodic_spline_to_component(e):
    """
    Export CLOSED/PERIODIC SPLINE as a SINGLE instruction component.

    This is required because EdgeSmith ignores closed/periodic splines.
    Your format supports this because n=1 works for 'spline' (it doesn't use x2,y2),
    and your make_shape() will add exactly one spline feature.
    """
    if e.dxftype() != "SPLINE":
        return None

    flags = _spline_flags(e)
    periodic = _spline_is_periodic(e)
    closed = _spline_is_closed(e) or periodic

    if not closed:
        return None

    degree = int(e.dxf.degree)
    uniq_knots, mults = _knot_multiplicities(list(e.knots))
    poles = [(float(Vec3(p).x), float(Vec3(p).y)) for p in e.control_points]
    weights = [float(w) for w in getattr(e, "weights", [])]

    if not poles:
        return None

    payload_dict = {
        "degree": degree,
        "poles": poles,
        "weights": weights,
        "knots": uniq_knots,
        "mults": mults,
        "periodic": bool(periodic),
        "flags": int(flags),
    }

    # single-vertex contour: use first pole as the "vertex" coordinate
    x0, y0 = float(poles[0][0]), float(poles[0][1])
    return {"pts": [["spline", x0, y0, payload_dict]]}


# -------------------------
# Loop ordering/orientation for open-edge loops
# -------------------------
def order_and_orient_closed_loop(loop_edges, gap_tol=1e-6):
    if len(loop_edges) < 2:
        return list(loop_edges)

    remaining = list(loop_edges)
    ordered = [remaining.pop(0)]
    cur_end = ordered[0].end
    start_start = ordered[0].start

    while remaining:
        nxt_idx = None
        nxt_edge = None
        for j, e in enumerate(remaining):
            if edgeminer.isclose(cur_end, e.start, gap_tol=gap_tol):
                nxt_idx, nxt_edge = j, e
                break
            if edgeminer.isclose(cur_end, e.end, gap_tol=gap_tol):
                nxt_idx, nxt_edge = j, e.reversed()
                break
        if nxt_edge is None:
            raise RuntimeError("Could not order loop: gap/junction or wrong loop selection.")
        remaining.pop(nxt_idx)
        ordered.append(nxt_edge)
        cur_end = nxt_edge.end

    if not edgeminer.isclose(cur_end, start_start, gap_tol=gap_tol):
        raise RuntimeError("Ordered chain did not close.")
    return ordered

# Fallback some closed splines are not recognized as closed
def _snap_key(xy, tol):
    # quantize to tolerance grid for stable node ids
    return (int(round(xy[0] / tol)), int(round(xy[1] / tol)))

def _build_endpoint_graph(entities, eps, tol):
    """
    Returns:
      nodes: dict node_id -> representative (x,y)
      adj:   dict node_id -> list of (neighbor_node_id, entity, a_xy, b_xy, flipped)
    """
    nodes = {}          # node_id -> (x,y) representative
    adj = {}            # node_id -> list of edges

    def add_node(xy):
        nid = _snap_key(xy, tol)
        if nid not in nodes:
            nodes[nid] = xy
            adj[nid] = []
        return nid

    for e in entities:
        a, b = eps[e]
        na = add_node(a)
        nb = add_node(b)
        # store both directions; flipped indicates whether traversal matches (a->b) or (b->a)
        adj[na].append((nb, e, a, b, False))
        adj[nb].append((na, e, b, a, True))
    return nodes, adj

def _extract_paths_from_graph(adj):
    """
    Extract one or more paths/cycles from an undirected multigraph.
    Output: list of paths, where each path is a list of directed edge records:
      (entity, start_xy, end_xy)
    """
    # track used edges by (entity_id, from_node, to_node) directed
    used = set()

    def edge_id(e, u, v):
        return (id(e), u, v)

    def degree(n):
        return len(adj.get(n, []))

    nodes = list(adj.keys())

    # Prefer starting at endpoints (degree 1) to get open paths first; cycles have all degree 2.
    start_nodes = [n for n in nodes if degree(n) == 1] + [n for n in nodes if degree(n) != 0 and degree(n) != 1]

    paths = []

    for start in start_nodes:
        for (nbr, e, a, b, flipped) in adj[start]:
            if edge_id(e, start, nbr) in used:
                continue

            path = []
            u = start
            v = nbr
            cur_entity = e
            cur_a = a
            cur_b = b

            # walk forward
            while True:
                used.add(edge_id(cur_entity, u, v))
                used.add(edge_id(cur_entity, v, u))
                path.append((cur_entity, cur_a, cur_b))

                # choose next edge from node v that is unused
                next_edge = None
                for (nbr2, e2, a2, b2, flipped2) in adj[v]:
                    if edge_id(e2, v, nbr2) in used:
                        continue
                    next_edge = (v, nbr2, e2, a2, b2)
                    break

                if next_edge is None:
                    break

                u, v, cur_entity, cur_a, cur_b = next_edge

                # stop if we closed a cycle back to start and there is no unused branch preference
                if v == start:
                    # check if any unused incident edges remain at start; if yes, keep going in another path later
                    break

            if path:
                paths.append(path)

    # There may be isolated cycles not reached above (all nodes degree 2)
    for start in nodes:
        # find any unused incident edge
        incident = None
        for (nbr, e, a, b, flipped) in adj[start]:
            if edge_id(e, start, nbr) not in used:
                incident = (start, nbr, e, a, b)
                break
        if incident is None:
            continue

        # walk cycle
        u, v, cur_entity, cur_a, cur_b = incident
        cycle = []
        while True:
            used.add(edge_id(cur_entity, u, v))
            used.add(edge_id(cur_entity, v, u))
            cycle.append((cur_entity, cur_a, cur_b))

            # advance
            next_edge = None
            for (nbr2, e2, a2, b2, flipped2) in adj[v]:
                if edge_id(e2, v, nbr2) in used:
                    continue
                next_edge = (v, nbr2, e2, a2, b2)
                break

            if next_edge is None:
                break
            u, v, cur_entity, cur_a, cur_b = next_edge
            if v == start:
                break

        if cycle:
            paths.append(cycle)

    return paths
def _any_spline_to_component(e):
    """
    Export ANY SPLINE (open/closed/periodic) as a SINGLE instruction component.
    Used as fallback when EdgeSmith produces no edges.
    """
    if e.dxftype() != "SPLINE":
        return None

    flags = _spline_flags(e)
    periodic = _spline_is_periodic(e)

    degree = int(e.dxf.degree)
    uniq_knots, mults = _knot_multiplicities(list(e.knots))
    poles = [(float(Vec3(p).x), float(Vec3(p).y)) for p in e.control_points]
    weights = [float(w) for w in getattr(e, "weights", [])]

    if not poles:
        return None

    payload_dict = {
        "degree": degree,
        "poles": poles,
        "weights": weights,
        "knots": uniq_knots,
        "mults": mults,
        "periodic": bool(periodic),
        "flags": int(flags),
    }

    # single-instruction contour: x,y can be any stable reference (use first pole)
    x0, y0 = float(poles[0][0]), float(poles[0][1])
    return {"pts": [["spline", x0, y0, payload_dict]]}


# -------------------------
# Convert an ordered edge chain into your {"pts":[...]} format
# -------------------------
def _chain_to_component(chain):
    pts = []
    for ed in chain:
        payload = getattr(ed, "payload", None)
        if payload is None:
            raise RuntimeError("Edge has no payload.")

        x1, y1 = _vxy(ed.start)
        x2, y2 = _vxy(ed.end)
        t = payload.dxftype()

        if t == "LINE":
            pts.append(["line", x1, y1])

        elif t == "ARC":
            c = Vec3(payload.dxf.center)
            cx, cy = float(c.x), float(c.y)
            direction = arc_direction_flag_for_shaper(payload, x1, y1, x2, y2)
            pts.append(["arc", x1, y1, cx, cy, direction])

        elif t == "ELLIPSE":
            cx, cy, a, b, phi_deg = _ellipse_params_from_dxf_ellipse(payload)
            inversed = not _traversal_is_ccw_xy(x1, y1, x2, y2, cx, cy)
            pts.append(["ellipse", x1, y1, cx, cy, a, b, phi_deg, inversed])

        elif t == "SPLINE":
            degree = int(payload.dxf.degree)
            uniq_knots, mults = _knot_multiplicities(list(payload.knots))
            poles = [(float(Vec3(p).x), float(Vec3(p).y)) for p in payload.control_points]
            weights = [float(w) for w in getattr(payload, "weights", [])]
            flags = int(payload.dxf.flags)

            periodic = _spline_is_periodic(payload)
            # DO NOT reject periodic splines anymore
            payload_dict = {
                "degree": degree,
                "poles": poles,
                "weights": weights,
                "knots": uniq_knots,
                "mults": mults,
                "periodic": bool(periodic),
                "flags": int(flags),
            }

            # Only attempt orientation fix for non-periodic (clamped/open) splines.
            if (not periodic) and poles:
                sp0 = Vec3(poles[0][0], poles[0][1], 0.0)
                sp1 = Vec3(poles[-1][0], poles[-1][1], 0.0)
                e0 = Vec3(ed.start.x, ed.start.y, 0.0)
                e1 = Vec3(ed.end.x, ed.end.y, 0.0)
                if (_dist2_xy(sp0, e1) + _dist2_xy(sp1, e0)) < (_dist2_xy(sp0, e0) + _dist2_xy(sp1, e1)):
                    _reverse_bspline_payload_inplace(payload_dict)

            pts.append(["spline", x1, y1, payload_dict])

        else:
            raise NotImplementedError(f"Unsupported payload type: {t}")

    return {"pts": pts}


def _component_start_xy(instr):
    if not isinstance(instr, (list, tuple)) or len(instr) < 3:
        return None
    try:
        return float(instr[1]), float(instr[2])
    except Exception:
        return None


def _normalize_component_pts(component, tol=1e-9):
    """
    Drop degenerate segments caused by consecutive duplicate start vertices.
    In this representation, entry i defines the segment from pts[i] to pts[i+1],
    so if those two starts coincide, the current entry is a zero-length segment.
    """
    pts = list(component.get("pts", []))
    if len(pts) < 2:
        return component, 0

    removed = 0
    while len(pts) >= 2:
        drop = set()
        n = len(pts)
        for i in range(n):
            a = _component_start_xy(pts[i])
            b = _component_start_xy(pts[(i + 1) % n])
            if a is None or b is None:
                continue
            if _sqdist_xy(a, b) <= tol * tol:
                drop.add(i)

        if not drop:
            break

        pts = [p for i, p in enumerate(pts) if i not in drop]
        removed += len(drop)
        if len(pts) < 2:
            break

    out = dict(component)
    out["pts"] = pts
    return out, removed


def _component_is_strictly_negative_x(component, tol=1e-9):
    """
    Conservative filter: return True only if the serialized contour is provably
    entirely on the x < 0 side using per-segment bounds.
    """
    pts = list(component.get("pts", []))
    if not pts:
        return False

    n = len(pts)
    for i, instr in enumerate(pts):
        kind = str(instr[0]).lower() if instr else ""
        a = _component_start_xy(instr)
        b = _component_start_xy(pts[(i + 1) % n])
        if a is None or b is None:
            return False

        if max(a[0], b[0]) >= -tol:
            return False

        if kind == "line":
            continue

        if kind == "spline":
            if len(instr) < 4 or not isinstance(instr[3], dict):
                return False
            poles = instr[3].get("poles", [])
            if not poles:
                return False
            if max(float(p[0]) for p in poles) >= -tol:
                return False
            continue

        if kind == "arc":
            if len(instr) < 5:
                return False
            cx, cy = float(instr[3]), float(instr[4])
            r = math.hypot(float(a[0]) - cx, float(a[1]) - cy)
            if (cx + r) >= -tol:
                return False
            continue

        if kind == "ellipse":
            if len(instr) < 8:
                return False
            cx = float(instr[3])
            a_len = float(instr[5])
            b_len = float(instr[6])
            phi = math.radians(float(instr[7]))
            xmax_offset = math.sqrt((a_len * math.cos(phi)) ** 2 + (b_len * math.sin(phi)) ** 2)
            if (cx + xmax_offset) >= -tol:
                return False
            continue

        # Unknown curve type: keep it rather than risk dropping valid geometry.
        return False

    return True
def _entity_endpoints_xy(e):
    t = e.dxftype()

    if t == "LINE":
        return _vxy(e.dxf.start), _vxy(e.dxf.end)

    if t == "ARC":
        c = Vec3(e.dxf.center)
        r = float(e.dxf.radius)
        a1 = math.radians(float(e.dxf.start_angle))
        a2 = math.radians(float(e.dxf.end_angle))
        p1 = (c.x + r * math.cos(a1), c.y + r * math.sin(a1), 0.0)
        p2 = (c.x + r * math.cos(a2), c.y + r * math.sin(a2), 0.0)
        return _vxy(p1), _vxy(p2)

    if t == "ELLIPSE":
        c = Vec3(e.dxf.center)
        maj = Vec3(e.dxf.major_axis)
        a = maj.magnitude
        b = a * float(e.dxf.ratio)
        phi = math.atan2(maj.y, maj.x)
        sp = float(e.dxf.start_param)
        ep = float(e.dxf.end_param)

        def pt(param):
            x = a * math.cos(param)
            y = b * math.sin(param)
            xr = math.cos(phi) * x - math.sin(phi) * y + c.x
            yr = math.sin(phi) * x + math.cos(phi) * y + c.y
            return (xr, yr, 0.0)

        return _vxy(pt(sp)), _vxy(pt(ep))

    if t == "SPLINE":
        # Prefer fit points if they exist (often cleaner in CAD exports)
        fit = list(getattr(e, "fit_points", []))
        if len(fit) >= 2:
            return _vxy(fit[0]), _vxy(fit[-1])

        cps = list(getattr(e, "control_points", []))
        if len(cps) >= 2:
            return _vxy(cps[0]), _vxy(cps[-1])

        return None

    return None

def _isclose_xy(a, b, tol):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 <= tol * tol


def _reverse_entity_instruction(entity, instr):
    """
    Return a reversed instruction for chaining.
    instr is what you push into pts: ["line", x,y], ["arc", x,y,cx,cy,dir], ["ellipse", ...], ["spline", x,y,payload]
    For reversal we flip direction flags where needed and reverse spline payload if non-periodic.
    """
    kind = instr[0]
    if kind == "line":
        return instr[:]  # nothing stored besides start vertex; reversal handled by chain order

    if kind == "arc":
        # reverse direction
        out = instr[:]
        out[5] = (not bool(out[5]))
        return out

    if kind == "ellipse":
        out = instr[:]
        out[8] = (not bool(out[8]))  # inversed toggled
        return out

    if kind == "spline":
        out = ["spline", instr[1], instr[2], dict(instr[3])]
        payload = out[3]
        if not bool(payload.get("periodic", False)):
            _reverse_bspline_payload_inplace(payload)
        return out

    return instr[:]


def _entity_to_vertex_instr(e, start_xy, end_xy):
    """
    Create your per-vertex instruction for the segment from start_xy to end_xy,
    including payload where needed. This mirrors your _chain_to_component logic.
    """
    x1, y1 = start_xy
    x2, y2 = end_xy
    t = e.dxftype()

    if t == "LINE":
        return ["line", x1, y1]

    if t == "ARC":
        c = Vec3(e.dxf.center)
        cx, cy = float(c.x), float(c.y)
        direction = arc_direction_flag_for_shaper(e, x1, y1, x2, y2)
        return ["arc", x1, y1, cx, cy, direction]

    if t == "ELLIPSE":
        cx, cy, a, b, phi_deg = _ellipse_params_from_dxf_ellipse(e)
        inversed = not _traversal_is_ccw_xy(x1, y1, x2, y2, cx, cy)
        return ["ellipse", x1, y1, cx, cy, a, b, phi_deg, inversed]

    if t == "SPLINE":
        degree = int(e.dxf.degree)
        uniq_knots, mults = _knot_multiplicities(list(e.knots))
        poles = [(float(Vec3(p).x), float(Vec3(p).y)) for p in e.control_points]
        weights = [float(w) for w in getattr(e, "weights", [])]
        flags = int(e.dxf.flags)
        periodic = _spline_is_periodic(e)

        payload_dict = {
            "degree": degree,
            "poles": poles,
            "weights": weights,
            "knots": uniq_knots,
            "mults": mults,
            "periodic": bool(periodic),
            "flags": int(flags),
        }

        # Orient spline payload to match evaluated endpoints
        if (not periodic) and poles:
            sp0 = Vec3(poles[0][0], poles[0][1], 0.0)
            sp1 = Vec3(poles[-1][0], poles[-1][1], 0.0)
            e0 = Vec3(x1, y1, 0.0)
            e1 = Vec3(x2, y2, 0.0)
            if (_dist2_xy(sp0, e1) + _dist2_xy(sp1, e0)) < (_dist2_xy(sp0, e0) + _dist2_xy(sp1, e1)):
                _reverse_bspline_payload_inplace(payload_dict)

        return ["spline", x1, y1, payload_dict]

    raise NotImplementedError(f"Unsupported entity in fallback chaining: {t}")


def _group_open_curves_by_endpoints(entities, tol):
    """
    Build connected components by endpoint snapping.
    """
    eps = {}
    for e in entities:
        ee = _entity_endpoints_xy(e)
        if ee is None:
            continue
        eps[e] = ee

    remaining = set(eps.keys())
    groups = []

    while remaining:
        seed = remaining.pop()
        group = [seed]
        frontier = [seed]

        while frontier:
            cur = frontier.pop()
            a0, a1 = eps[cur]
            for other in list(remaining):
                b0, b1 = eps[other]
                if (_isclose_xy(a0, b0, tol) or _isclose_xy(a0, b1, tol) or
                    _isclose_xy(a1, b0, tol) or _isclose_xy(a1, b1, tol)):
                    remaining.remove(other)
                    group.append(other)
                    frontier.append(other)

        groups.append(group)

    return groups, eps


def _order_and_orient_entities_into_chain(group, eps, tol):
    """
    Order entities head-to-tail. Returns list of (entity, start_xy, end_xy) in chain order.
    """
    if not group:
        return []

    # pick an arbitrary start
    chain = []
    used = set()

    cur = group[0]
    used.add(cur)
    s, t = eps[cur]
    chain.append((cur, s, t))
    cur_end = t

    while len(used) < len(group):
        found = False
        for e in group:
            if e in used:
                continue
            a, b = eps[e]
            if _isclose_xy(cur_end, a, tol):
                chain.append((e, a, b))
                cur_end = b
                used.add(e)
                found = True
                break
            if _isclose_xy(cur_end, b, tol):
                chain.append((e, b, a))  # reversed
                cur_end = a
                used.add(e)
                found = True
                break
        if not found:
            # Not a single chain (branching or gap). Return what we have; caller can decide to skip/log.
            break

    return chain
def _curve_endpoints(e):
    t = e.dxftype()
    if t == "SPLINE":
        cps = list(e.control_points)
        if len(cps) < 2:
            return None
        return Vec3(cps[0]), Vec3(cps[-1])

    if t == "LINE":
        return Vec3(e.dxf.start), Vec3(e.dxf.end)

    if t == "ARC":
        c = Vec3(e.dxf.center)
        r = float(e.dxf.radius)
        a1 = math.radians(e.dxf.start_angle)
        a2 = math.radians(e.dxf.end_angle)
        p1 = Vec3(c.x + r * math.cos(a1), c.y + r * math.sin(a1), 0)
        p2 = Vec3(c.x + r * math.cos(a2), c.y + r * math.sin(a2), 0)
        return p1, p2

    return None
def _group_entities_by_connectivity(entities, tol=1e-6):
    groups = []
    unused = set(entities)

    def close(a, b):
        return (a - b).magnitude <= tol

    while unused:
        seed = unused.pop()
        ep = _curve_endpoints(seed)
        if ep is None:
            groups.append([seed])
            continue

        group = [seed]
        changed = True

        while changed:
            changed = False
            for e in list(unused):
                eep = _curve_endpoints(e)
                if eep is None:
                    continue
                for g in group:
                    gep = _curve_endpoints(g)
                    if gep is None:
                        continue
                    if (close(eep[0], gep[0]) or close(eep[0], gep[1]) or
                        close(eep[1], gep[0]) or close(eep[1], gep[1])):
                        unused.remove(e)
                        group.append(e)
                        changed = True
                        break

        groups.append(group)

    return groups
def _sqdist_xy(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return dx * dx + dy * dy


def _endpoints_for_chainables(chainables):
    """
    Returns:
      eps: dict entity -> [start_xy, end_xy]  (mutable list to allow snapping updates)
      endpoints: list of (entity, end_idx) where end_idx in {0,1}
    """
    eps = {}
    endpoints = []
    for e in chainables:
        ee = _entity_endpoints_xy(e)
        if ee is None:
            continue
        eps[e] = [ee[0], ee[1]]
        endpoints.append((e, 0))
        endpoints.append((e, 1))
    return eps, endpoints


def _bruteforce_snap_endpoints(eps, endpoints, max_snap_dist):
    """
    Greedily snap the closest two endpoints until exhausted.
    Snaps by replacing both endpoints with their midpoint.

    max_snap_dist: hard cap; if the closest pair exceeds this, stop.
    """
    max_d2 = max_snap_dist * max_snap_dist

    # active endpoints (each endpoint can be updated multiple times)
    active = list(endpoints)

    while len(active) >= 2:
        best_i = best_j = None
        best_d2 = None

        # O(N^2) brute force
        for i in range(len(active) - 1):
            ei, ki = active[i]
            pi = eps[ei][ki]
            for j in range(i + 1, len(active)):
                ej, kj = active[j]
                pj = eps[ej][kj]
                d2 = _sqdist_xy(pi, pj)
                if best_d2 is None or d2 < best_d2:
                    best_d2 = d2
                    best_i, best_j = i, j

        if best_d2 is None or best_d2 > max_d2:
            break

        e1, k1 = active[best_i]
        e2, k2 = active[best_j]
        p1 = eps[e1][k1]
        p2 = eps[e2][k2]

        mid = ((p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5)

        eps[e1][k1] = mid
        eps[e2][k2] = mid

        # remove the snapped endpoints from further matching
        # remove higher index first
        for idx in sorted([best_i, best_j], reverse=True):
            active.pop(idx)

    return eps


def _order_closed_cycle_greedy(entities, eps, tol):
    """
    Order entities into one closed cycle by greedy endpoint matching.
    Assumes the contour is closed and non-branching.

    Returns list of (entity, start_xy, end_xy) in traversal order.
    """
    if not entities:
        return []

    rem = list(entities)
    first = rem.pop(0)
    s0, t0 = eps[first][0], eps[first][1]
    chain = [(first, s0, t0)]
    cur_end = t0

    while rem:
        found = False
        for i, e in enumerate(rem):
            a, b = eps[e][0], eps[e][1]
            if _isclose_xy(cur_end, a, tol):
                chain.append((e, a, b))
                cur_end = b
                rem.pop(i)
                found = True
                break
            if _isclose_xy(cur_end, b, tol):
                chain.append((e, b, a))  # reverse
                cur_end = a
                rem.pop(i)
                found = True
                break
        if not found:
            raise RuntimeError("Bruteforce snap produced endpoints, but ordering still failed.")

    # closure check (cycle)
    if not _isclose_xy(cur_end, chain[0][1], tol):
        raise RuntimeError("Ordered chain did not close after bruteforce snapping.")

    return chain

def _open_chain_with_manual_closure(entities, eps, tol):
    """
    Last-resort fallback:
    - order a single open chain
    - append one synthetic line from chain end back to chain start
    Returns {"pts":[...]} or None.
    """
    chain = _order_and_orient_entities_into_chain(entities, eps, tol)
    if not chain or len(chain) != len(entities):
        return None

    pts = []
    for ent, a, b in chain:
        try:
            pts.append(_entity_to_vertex_instr(ent, a, b))
        except Exception:
            return None

    start_xy = chain[0][1]
    end_xy = chain[-1][2]
    if not _isclose_xy(start_xy, end_xy, tol):
        # make_shape() interprets each item as a segment start to next vertex;
        # this extra vertex creates an explicit closing segment end->start.
        pts.append(["line", float(end_xy[0]), float(end_xy[1])])

    return {"pts": pts} if pts else None

def _fallback_components_from_chainables(
    chainables,
    gap_tol,
    allow_manual_line_closure=None,
    max_force_merge_dist_override=None,
):
    """
    Build components from chainable entities by endpoint snapping + closed-cycle ordering.
    Tries progressively larger snap tolerances to tolerate endpoint drift from CAD/DXF export.
    Returns:
      (components, fail_reason)
      - components: list (possibly empty)
      - fail_reason: None if components were produced, else a short reason string
    """
    # Cap endpoint snapping to avoid over-merging distant/open contours.
    merge_cap = (
        float(MAX_FORCE_MERGE_DIST)
        if max_force_merge_dist_override is None
        else float(max_force_merge_dist_override)
    )
    base_cap = min(float(gap_tol), merge_cap)
    snap_caps = []
    for c in (base_cap, merge_cap):
        c = float(max(base_cap, c))
        if c not in snap_caps:
            snap_caps.append(c)

    last_error = None
    do_manual_line_closure = (
        bool(ENABLE_MANUAL_LINE_CLOSURE)
        if allow_manual_line_closure is None
        else bool(allow_manual_line_closure)
    )
    for snap_cap in snap_caps:
        eps, endpoints = _endpoints_for_chainables(chainables)
        if not eps:
            return [], "no_endpoints"

        eps = _bruteforce_snap_endpoints(eps, endpoints, max_snap_dist=snap_cap)
        attempt_components = []
        ok = True

        # Peel connected sets and require each set to form a closed cycle.
        remaining = list(eps.keys())
        while remaining:
            seed = remaining[0]
            seed_set = [seed]
            remaining.pop(0)
            changed = True
            while changed:
                changed = False
                for e in list(remaining):
                    ea0, ea1 = eps[e][0], eps[e][1]
                    for g in seed_set:
                        ga0, ga1 = eps[g][0], eps[g][1]
                        if (ea0 == ga0 or ea0 == ga1 or ea1 == ga0 or ea1 == ga1):
                            seed_set.append(e)
                            remaining.remove(e)
                            changed = True
                            break

            try:
                chain = _order_closed_cycle_greedy(seed_set, eps, tol=snap_cap)
                pts = []
                for (ent, a, b) in chain:
                    try:
                        pts.append(_entity_to_vertex_instr(ent, a, b))
                    except Exception as e:
                        last_error = f"entity_to_vertex_failed:{type(e).__name__}"
                        ok = False
                        break
                if not ok:
                    break
                if pts:
                    attempt_components.append({"pts": pts})
            except Exception as e:
                last_error = f"order_or_close_failed:{type(e).__name__}"
                if do_manual_line_closure:
                    # Optional last fallback: force-close a single open chain with one line.
                    comp = _open_chain_with_manual_closure(seed_set, eps, tol=snap_cap)
                    if comp is None:
                        last_error = "manual_line_closure_failed"
                        ok = False
                        break
                    attempt_components.append(comp)
                else:
                    # Keep open contours open: do not inject synthetic closing lines.
                    last_error = "open_chain_no_manual_closure"
                    ok = False
                    break

        if ok and attempt_components:
            return attempt_components, None

    if last_error is None:
        last_error = "no_closed_cycles_after_snapping"
    return [], last_error

# -------------------------
# Main extraction: closed primitives + looped open-edge networks
# -------------------------
def _entities_to_components(entities, gap_tol=GAP_TOL):
    components = []
    open_entities = []
    diag = {
        "input_entities": len(entities),
        "input_types": dict(Counter([e.dxftype() for e in entities])),
        "closed_primitive_components": 0,
        "open_entities": 0,
        "edges_count": 0,
        "loop_components": 0,
        "fallback_attempted": False,
        "chainables_count": 0,
        "fallback_components": 0,
        "fallback_reason": None,
        "drop_reason": None,
    }

    for e in entities:
        t = e.dxftype()

        # Closed primitives EdgeSmith won't give edges for:
        if t == "CIRCLE":
            components.append(_circle_to_component(e))
            diag["closed_primitive_components"] += 1
            continue

        if t == "LWPOLYLINE" and bool(getattr(e, "closed", False)):
            comp = _lwpolyline_closed_to_component(e)
            if comp is not None and comp.get("pts"):
                components.append(comp)
                diag["closed_primitive_components"] += 1
                continue

        if t == "POLYLINE" and bool(getattr(e, "is_closed", False)):
            comp = _polyline_closed_to_component(e)
            if comp is not None and comp.get("pts"):
                components.append(comp)
                diag["closed_primitive_components"] += 1
                continue

        if t == "ELLIPSE" and _is_full_ellipse(e):
            comp = _full_ellipse_to_component(e)
            if comp is not None and comp.get("pts"):
                components.append(comp)
                diag["closed_primitive_components"] += 1
                continue

        # Closed/periodic SPLINEs are ignored by EdgeSmith -> export directly as single-object contour
        if t == "SPLINE" and (_spline_is_closed(e) or _spline_is_periodic(e)):
            comp = _closed_or_periodic_spline_to_component(e)
            if comp is not None and comp.get("pts"):
                components.append(comp)
                diag["closed_primitive_components"] += 1
                continue

        # Everything else goes through edgesmith/edgeminer
        open_entities.append(e)

    diag["open_entities"] = len(open_entities)
    # Build edges for loop detection
    edges = list(edgesmith.edges_from_entities_2d(open_entities))
    diag["edges_count"] = len(edges)

    components_from_loops = 0
    if edges:
        deposit = edgeminer.Deposit(edges, gap_tol=gap_tol)
        loops = list(edgeminer.find_all_loops(deposit))
        if not loops:
            loop = edgeminer.find_loop(deposit)
            if loop is not None:
                loops = [loop]

        for loop_edges in loops:
            chain = order_and_orient_closed_loop(loop_edges, gap_tol=gap_tol)
            comp = _chain_to_component(chain)
            if comp is not None and comp.get("pts"):
                components.append(comp)
                components_from_loops += 1
                diag["loop_components"] += 1

    # FALLBACK: if no loop-based components were produced, still export any SPLINE geometry
    # (common when edges exist but do not form closed loops, while splines are present)
    if components_from_loops == 0:
        chainables = [e for e in open_entities if e.dxftype() in {"SPLINE", "LINE", "ARC", "ELLIPSE"}]
        diag["fallback_attempted"] = bool(chainables)
        diag["chainables_count"] = len(chainables)
        if chainables:
            fallback_components, fallback_reason = _fallback_components_from_chainables(
                chainables,
                gap_tol=gap_tol,
                allow_manual_line_closure=False,
            )
            components.extend(fallback_components)
            diag["fallback_components"] = len(fallback_components)
            diag["fallback_reason"] = fallback_reason

    if not components:
        if diag["fallback_attempted"]:
            diag["drop_reason"] = diag["fallback_reason"] or "fallback_failed"
        elif diag["edges_count"] > 0:
            diag["drop_reason"] = "no_loops_found"
        elif diag["open_entities"] > 0:
            diag["drop_reason"] = "no_edges_from_open_entities"
        else:
            diag["drop_reason"] = "no_supported_entities"

    return components, diag


def dxf_to_components(msp, gap_tol=GAP_TOL):
    # Expand INSERTs to visible geometry
    entities = []
    for e in msp:
        if e.dxftype() == "INSERT":
            entities.extend(list(e.virtual_entities()))
        else:
            entities.append(e)

    # Process layer-wise so open-chain fallback can run for layers that would
    # otherwise be masked by successful loop extraction in other layers.
    by_layer = {}
    for e in entities:
        layer = str(getattr(getattr(e, "dxf", None), "layer", "0"))
        by_layer.setdefault(layer, []).append(e)

    component_records = []
    layer_reports = []
    for layer in sorted(by_layer.keys()):
        layer_entities = by_layer[layer]
        comps, diag = _entities_to_components(layer_entities, gap_tol=gap_tol)
        for comp in comps:
            component_records.append(
                {
                    "component": comp,
                    "layer": layer,
                }
            )
        layer_reports.append(
            {
                "layer": layer,
                "component_count": len(comps),
                "diag": diag,
            }
        )

    return component_records, entities, layer_reports, by_layer


# -------------------------
# File IO helpers
# -------------------------
def _dxf_to_base_json_path(dxf_path):
    dxf_path = os.path.abspath(dxf_path)
    d = os.path.dirname(dxf_path)
    stem = os.path.splitext(os.path.basename(dxf_path))[0]
    return os.path.join(d, stem.replace(".", "").replace("_", "").replace("Feature", "") + ".json")

def _indexed_dump_paths(base_json_path, count):
    base_json_path = os.path.abspath(base_json_path)
    root, ext = os.path.splitext(base_json_path)
    if ext.lower() != ".json":
        ext = ".json"
    if count <= 1:
        return [root + ext]
    return [f"{root}_{i}{ext}" for i in range(count)]


def _layer_group_index(layer_name):
    """
    Extract numeric sub-index from DXF layer names like:
      ...__G056_CLOSED, ...__G072_OPEN, ...__G001_WIRE_FAIL
    Returns int or None.
    """
    layer_name = str(layer_name or "")
    m = re.search(r"__G(\d+)_(?:OPEN|CLOSED|WIRE_FAIL)\b", layer_name)
    if not m:
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None


def _indexed_dump_paths_by_layer(base_json_path, component_records):
    """
    Name outputs by DXF group index when available:
      <base>_<idx>.json
    If multiple components share an idx, append sub-index:
      <base>_<idx>_1.json, <base>_<idx>_2.json, ...
    If idx is unavailable, fall back to a sequential index.
    """
    base_json_path = os.path.abspath(base_json_path)
    root, ext = os.path.splitext(base_json_path)
    if ext.lower() != ".json":
        ext = ".json"

    used = set()
    fallback_counter = 0
    out = []
    for rec in component_records:
        layer = rec.get("layer")
        idx = _layer_group_index(layer)
        if idx is None:
            stem = f"{root}_{fallback_counter}"
            fallback_counter += 1
        else:
            stem = f"{root}_{idx}"

        path = f"{stem}{ext}"
        if path in used:
            sub = 1
            while True:
                cand = f"{stem}_{sub}{ext}"
                if cand not in used:
                    path = cand
                    break
                sub += 1

        used.add(path)
        out.append(path)

    return out


def _layer_is_open_group(layer_name):
    return re.search(r"__G\d+_OPEN\b", str(layer_name or "")) is not None


def _choose_canonical_layer(source_layers):
    """
    Choose a stable layer label for rescued components merged from multiple OPEN groups.
    Prefer the smallest numeric G-index when available.
    """
    if not source_layers:
        return "0"

    keyed = []
    for layer in source_layers:
        idx = _layer_group_index(layer)
        keyed.append(((idx if idx is not None else 10**9), str(layer), layer))
    keyed.sort()
    return keyed[0][2]


def _rescue_open_layers_by_endpoints(layer_reports, by_layer, gap_tol):
    """
    Second-pass rescue for dropped OPEN groups:
    - pool chainable entities from dropped OPEN layers
    - group by endpoint connectivity
    - try to close each group using the existing fallback solver

    Returns:
      rescued_records: [{"component":..., "layer":...}, ...]
      rescue_reports: list[dict]
      rescued_source_layers: set[str]
    """
    candidate_layers = []
    for r in layer_reports:
        if int(r.get("component_count", 0)) != 0:
            continue
        layer = r.get("layer", "")
        diag = r.get("diag", {})
        reason = str(diag.get("drop_reason", ""))
        if _layer_is_open_group(layer) and reason.startswith("open_chain"):
            candidate_layers.append(layer)

    if not candidate_layers:
        return [], [], set()

    rescued_records = []
    rescue_reports = []
    rescued_source_layers = set()

    # Pass 1: rescue each dropped OPEN layer independently.
    layer_chainables_map = {}
    unresolved_layers = set()
    for layer in sorted(candidate_layers):
        layer_chainables = []
        for e in by_layer.get(layer, []):
            if e.dxftype() not in {"SPLINE", "LINE", "ARC", "ELLIPSE"}:
                continue
            if _entity_endpoints_xy(e) is None:
                continue
            layer_chainables.append(e)
        layer_chainables_map[layer] = layer_chainables

        if not layer_chainables:
            rescue_reports.append(
                {
                    "status": "failed",
                    "method": "per_layer",
                    "layer": layer,
                    "entity_count": 0,
                    "reason": "no_chainables",
                }
            )
            unresolved_layers.add(layer)
            continue

        comps, reason = _fallback_components_from_chainables(
            layer_chainables,
            gap_tol=gap_tol,
            allow_manual_line_closure=True,
        )
        if comps:
            # Reject degenerate single-line "rescues" here; defer to pooled rescue.
            weak = any(len(c.get("pts", [])) < 3 for c in comps)
            if weak:
                rescue_reports.append(
                    {
                        "status": "weak",
                        "method": "per_layer",
                        "layer": layer,
                        "entity_count": len(layer_chainables),
                        "reason": "weak_component_too_small",
                    }
                )
                unresolved_layers.add(layer)
                continue
            for comp in comps:
                rescued_records.append(
                    {
                        "component": comp,
                        "layer": layer,
                        "rescued": True,
                        "source_layers": [layer],
                        "method": "per_layer",
                    }
                )
            rescue_reports.append(
                {
                    "status": "rescued",
                    "method": "per_layer",
                    "layer": layer,
                    "entity_count": len(layer_chainables),
                    "component_count": len(comps),
                }
            )
            rescued_source_layers.add(layer)
        else:
            rescue_reports.append(
                {
                    "status": "failed",
                    "method": "per_layer",
                    "layer": layer,
                    "entity_count": len(layer_chainables),
                    "reason": reason or "rescue_failed",
                }
            )
            unresolved_layers.add(layer)

    # Pass 1b: targeted orphan pairing (e.g. single detached line layer with its nearest main layer).
    rescue_bridge_cap = float(max(3.5, 10.0 * gap_tol, float(MAX_FORCE_MERGE_DIST)))
    tiny_layers = [
        l for l in sorted(unresolved_layers)
        if 0 < len(layer_chainables_map.get(l, [])) <= 2
    ]
    for tiny in tiny_layers:
        if tiny not in unresolved_layers:
            continue
        tiny_entities = layer_chainables_map.get(tiny, [])
        tiny_pts = []
        for e in tiny_entities:
            ee = _entity_endpoints_xy(e)
            if ee is not None:
                tiny_pts.extend([ee[0], ee[1]])
        if not tiny_pts:
            continue

        best_layer = None
        best_d2 = None
        for other in sorted(unresolved_layers):
            if other == tiny:
                continue
            other_entities = layer_chainables_map.get(other, [])
            if len(other_entities) <= 2:
                continue
            other_pts = []
            for e in other_entities:
                ee = _entity_endpoints_xy(e)
                if ee is not None:
                    other_pts.extend([ee[0], ee[1]])
            if not other_pts:
                continue

            d2 = min(_sqdist_xy(a, b) for a in tiny_pts for b in other_pts)
            if best_d2 is None or d2 < best_d2:
                best_d2 = d2
                best_layer = other

        if best_layer is None:
            rescue_reports.append(
                {
                    "status": "failed",
                    "method": "pair",
                    "tiny_layer": tiny,
                    "reason": "no_pair_candidate",
                }
            )
            continue

        if best_d2 is None or best_d2 > (rescue_bridge_cap * rescue_bridge_cap):
            rescue_reports.append(
                {
                    "status": "failed",
                    "method": "pair",
                    "tiny_layer": tiny,
                    "target_layer": best_layer,
                    "distance": (math.sqrt(best_d2) if best_d2 is not None else None),
                    "reason": "pair_too_far",
                }
            )
            continue

        combo = list(tiny_entities) + list(layer_chainables_map.get(best_layer, []))
        comps, reason = _fallback_components_from_chainables(
            combo,
            gap_tol=gap_tol,
            allow_manual_line_closure=True,
            max_force_merge_dist_override=rescue_bridge_cap,
        )
        if comps:
            canonical_layer = best_layer
            for comp in comps:
                rescued_records.append(
                    {
                        "component": comp,
                        "layer": canonical_layer,
                        "rescued": True,
                        "source_layers": [best_layer, tiny],
                        "method": "pair",
                    }
                )
            rescue_reports.append(
                {
                    "status": "rescued",
                    "method": "pair",
                    "tiny_layer": tiny,
                    "target_layer": best_layer,
                    "distance": math.sqrt(best_d2),
                    "component_count": len(comps),
                }
            )
            rescued_source_layers.add(tiny)
            rescued_source_layers.add(best_layer)
            unresolved_layers.discard(tiny)
            unresolved_layers.discard(best_layer)
        else:
            rescue_reports.append(
                {
                    "status": "failed",
                    "method": "pair",
                    "tiny_layer": tiny,
                    "target_layer": best_layer,
                    "distance": math.sqrt(best_d2),
                    "reason": reason or "pair_rescue_failed",
                }
            )

    # Pass 2: pooled connectivity rescue for unresolved OPEN layers.
    # Build pooled groups only from unresolved layers to avoid over-merging.
    pooled_chainables = []
    entity_layer = {}
    for layer in sorted(unresolved_layers):
        for e in by_layer.get(layer, []):
            if e.dxftype() not in {"SPLINE", "LINE", "ARC", "ELLIPSE"}:
                continue
            if _entity_endpoints_xy(e) is None:
                continue
            pooled_chainables.append(e)
            entity_layer[e] = layer

    if pooled_chainables and unresolved_layers:
        connect_tol = float(max(gap_tol, 1e-6))
        conn_groups, _ = _group_open_curves_by_endpoints(pooled_chainables, tol=connect_tol)

        for gi, group in enumerate(conn_groups):
            source_layers = sorted({entity_layer[e] for e in group})
            comps, reason = _fallback_components_from_chainables(
                group,
                gap_tol=gap_tol,
                allow_manual_line_closure=True,
                max_force_merge_dist_override=float(MAX_FORCE_MERGE_DIST),
            )
            if comps:
                canonical_layer = _choose_canonical_layer(source_layers)
                # Replace any per-layer rescue records participating in this pooled merge.
                rescued_records = [
                    rec for rec in rescued_records
                    if not (rec.get("method") == "per_layer" and rec.get("layer") in set(source_layers))
                ]
                for comp in comps:
                    rescued_records.append(
                        {
                            "component": comp,
                            "layer": canonical_layer,
                            "rescued": True,
                            "source_layers": list(source_layers),
                            "method": "pooled",
                        }
                    )
                rescue_reports.append(
                    {
                        "status": "rescued",
                        "method": "pooled",
                        "group_index": gi,
                        "source_layers": source_layers,
                        "canonical_layer": canonical_layer,
                        "entity_count": len(group),
                        "component_count": len(comps),
                    }
                )
                for l in source_layers:
                    rescued_source_layers.add(l)
                unresolved_layers.difference_update(source_layers)
            else:
                rescue_reports.append(
                    {
                        "status": "failed",
                        "method": "pooled",
                        "group_index": gi,
                        "source_layers": source_layers,
                        "entity_count": len(group),
                        "reason": reason or "rescue_failed",
                    }
                )

    return rescued_records, rescue_reports, rescued_source_layers


# -------------------------
# Processing
# -------------------------
def process_one_file(dxf_path):
    try:
        doc = ezdxf.readfile(dxf_path)
        msp = doc.modelspace()
    except Exception as e:
        _log_exception(FAILED_FILES_LOG, header=f"[ERROR] readfile file={dxf_path}", exc=e)
        print(f"[ERROR] readfile: {dxf_path}: {e}")
        return

    try:
        component_records, all_entities, layer_reports, by_layer = dxf_to_components(
            msp, gap_tol=GAP_TOL
        )
    except Exception as e:
        _log_exception(FAILED_FILES_LOG, header=f"[ERROR] extract file={dxf_path}", exc=e)
        print(f"[ERROR] extract: {dxf_path}: {e}")
        return

    rescued_records, rescue_reports, rescued_source_layers = _rescue_open_layers_by_endpoints(
        layer_reports=layer_reports,
        by_layer=by_layer,
        gap_tol=GAP_TOL,
    )
    rescued_total = len(rescued_records)
    rescued_group_count = sum(1 for rr in rescue_reports if rr.get("status") == "rescued")
    print(
        f"[RESCUE_OBJECTS] file={dxf_path} rescued_objects={rescued_total} "
        f"rescued_groups={rescued_group_count} rescue_attempt_groups={len(rescue_reports)}"
    )
    _log_line(
        FAILED_OBJECTS_LOG,
        f"[RESCUE_OBJECTS] file={dxf_path} rescued_objects={rescued_total} "
        f"rescued_groups={rescued_group_count} rescue_attempt_groups={len(rescue_reports)}",
    )
    if rescue_reports:
        rescued_total_from_reports = sum(
            int(rr.get("component_count", 0))
            for rr in rescue_reports
            if rr.get("status") == "rescued"
        )
        print(
            f"[RESCUE_SUMMARY] file={dxf_path} groups={len(rescue_reports)} "
            f"rescued_components={rescued_total_from_reports}"
        )
        for rr in rescue_reports:
            method = rr.get("method", "pooled")
            if method == "per_layer":
                if rr.get("status") == "rescued":
                    msg = (
                        f"[RESCUE_LAYER] file={dxf_path} layer={rr['layer']} status=rescued "
                        f"entities={rr['entity_count']} components={rr['component_count']}"
                    )
                elif rr.get("status") == "weak":
                    msg = (
                        f"[RESCUE_LAYER] file={dxf_path} layer={rr['layer']} status=weak "
                        f"entities={rr['entity_count']} reason={rr.get('reason','unknown')}"
                    )
                else:
                    msg = (
                        f"[RESCUE_LAYER] file={dxf_path} layer={rr['layer']} status=failed "
                        f"entities={rr['entity_count']} reason={rr.get('reason','unknown')}"
                    )
            elif method == "pair":
                if rr.get("status") == "rescued":
                    msg = (
                        f"[RESCUE_PAIR] file={dxf_path} tiny={rr['tiny_layer']} "
                        f"target={rr['target_layer']} status=rescued "
                        f"distance={rr.get('distance')} components={rr['component_count']}"
                    )
                else:
                    msg = (
                        f"[RESCUE_PAIR] file={dxf_path} tiny={rr.get('tiny_layer')} "
                        f"target={rr.get('target_layer')} status=failed "
                        f"distance={rr.get('distance')} reason={rr.get('reason','unknown')}"
                    )
            elif rr.get("status") == "rescued":
                msg = (
                    f"[RESCUE_GROUP] file={dxf_path} group={rr['group_index']} status=rescued "
                    f"entities={rr['entity_count']} components={rr['component_count']} "
                    f"canonical_layer={rr['canonical_layer']} source_layers={rr['source_layers']}"
                )
            else:
                msg = (
                    f"[RESCUE_GROUP] file={dxf_path} group={rr['group_index']} status=failed "
                    f"entities={rr['entity_count']} reason={rr.get('reason','unknown')} "
                    f"source_layers={rr['source_layers']}"
                )
            print(msg)
            _log_line(FAILED_OBJECTS_LOG, msg)

    if rescued_records:
        component_records.extend(rescued_records)

    dropped_layers = [
        r for r in layer_reports
        if int(r.get("component_count", 0)) == 0 and r.get("layer") not in rescued_source_layers
    ]
    if dropped_layers:
        print(
            f"[DROP_SUMMARY] file={dxf_path} dropped_layers={len(dropped_layers)} "
            f"total_layers={len(layer_reports)}"
        )
        for r in dropped_layers:
            layer = r.get("layer", "<unknown>")
            diag = r.get("diag", {})
            reason = diag.get("drop_reason", "unknown")
            chainables = int(diag.get("chainables_count", 0))
            fallback_reason = diag.get("fallback_reason")
            types = diag.get("input_types", {})
            msg = (
                f"[DROP_LAYER] file={dxf_path} layer={layer} reason={reason} "
                f"entities={diag.get('input_entities', 0)} chainables={chainables}"
            )
            if fallback_reason:
                msg += f" fallback_reason={fallback_reason}"
            msg += f" types={types}"
            print(msg)
            _log_line(FAILED_OBJECTS_LOG, msg)

    if not component_records:
        types = Counter([e.dxftype() for e in all_entities])
        reason = _infer_no_components_reason(types)
        _log_line(
            FAILED_FILES_LOG,
            "[FAILED_FILE] no components produced "
            f"file={dxf_path} "
            f"reason={reason} "
            f"gap_tol={GAP_TOL} "
            f"entity_types={dict(types)}"
        )
        print(f"[SKIP] no components: {dxf_path}")
        return

    base_json = _dxf_to_base_json_path(dxf_path)
    os.makedirs(os.path.dirname(base_json), exist_ok=True)
    paths = _indexed_dump_paths_by_layer(base_json, component_records)

    # Remove stale outputs from previous naming/indexing schemes for this DXF base.
    root, ext = os.path.splitext(base_json)
    stale_patterns = [f"{root}{ext}", f"{root}_*{ext}"]
    stale_removed = 0
    for pat in stale_patterns:
        for oldp in glob.glob(pat):
            try:
                os.remove(oldp)
                stale_removed += 1
            except Exception as e:
                _log_exception(
                    FAILED_OBJECTS_LOG,
                    header=f"[WARN] failed_to_remove_stale file={dxf_path} stale={oldp}",
                    exc=e,
                )
                print(f"[WARN] failed to remove stale output {oldp}: {e}")

    print(f"[CLEAN_OUTPUTS] file={dxf_path} removed_stale={stale_removed}")
    _log_line(FAILED_OBJECTS_LOG, f"[CLEAN_OUTPUTS] file={dxf_path} removed_stale={stale_removed}")

    written = 0
    rescued_written = 0
    skipped_negative_x = 0
    for i, (rec, path) in enumerate(zip(component_records, paths)):
        comp = rec.get("component", {})
        layer = rec.get("layer", "<unknown>")
        try:
            comp, removed_zero_segments = _normalize_component_pts(comp)
            if removed_zero_segments:
                msg = (
                    f"[NORMALIZE_COMPONENT] file={dxf_path} layer={layer} out={path} "
                    f"removed_zero_segments={removed_zero_segments}"
                )
                print(msg)
                _log_line(FAILED_OBJECTS_LOG, msg)
            if not comp.get("pts"):
                raise RuntimeError("component has empty pts")
            if _component_is_strictly_negative_x(comp):
                skipped_negative_x += 1
                msg = f"[SKIP_NEGATIVE_X] file={dxf_path} layer={layer} out={path}"
                print(msg)
                _log_line(FAILED_OBJECTS_LOG, msg)
                continue
            with open(path, "w", encoding="utf-8") as f:
                json.dump(comp, f, indent=2, ensure_ascii=False)
            print(path)
            written += 1
            if rec.get("rescued", False):
                rescued_written += 1
        except Exception as e:
            _log_exception(
                FAILED_OBJECTS_LOG,
                header=f"[FAILED_OBJECT] file={dxf_path} comp_index={i} layer={layer} out={path}",
                exc=e,
            )
            print(f"[WARN] failed component {i} layer={layer} in {dxf_path}: {e}")

    print(
        f"[RESCUE_WRITTEN] file={dxf_path} rescued_written={rescued_written} "
        f"rescued_created={rescued_total}"
    )
    _log_line(
        FAILED_OBJECTS_LOG,
        f"[RESCUE_WRITTEN] file={dxf_path} rescued_written={rescued_written} "
        f"rescued_created={rescued_total}",
    )
    print(f"[SKIP_NEGATIVE_X_SUMMARY] file={dxf_path} skipped={skipped_negative_x}")
    _log_line(
        FAILED_OBJECTS_LOG,
        f"[SKIP_NEGATIVE_X_SUMMARY] file={dxf_path} skipped={skipped_negative_x}",
    )

    if written == 0:
        _log_line(FAILED_FILES_LOG, f"[FAILED_FILE] no components written file={dxf_path}")


def main():
    global ROOT_DIR, DXF_GLOB, GAP_TOL, FAILED_FILES_LOG, FAILED_OBJECTS_LOG

    args = _parse_args()
    ROOT_DIR = os.path.abspath(args.root_dir)
    DXF_GLOB = args.dxf_glob
    GAP_TOL = float(args.gap_tol)
    FAILED_FILES_LOG = os.path.join(ROOT_DIR, "failed_files.log")
    FAILED_OBJECTS_LOG = os.path.join(ROOT_DIR, "failed_objects.log")

    _reset_log(FAILED_FILES_LOG)
    _reset_log(FAILED_OBJECTS_LOG)

    if args.dxf_file:
        dxf_paths = sorted(os.path.abspath(p) for p in args.dxf_file)
    else:
        pattern = os.path.join(ROOT_DIR, DXF_GLOB)
        dxf_paths = sorted(glob.glob(pattern, recursive=True))

    if args.path_regex:
        rx = re.compile(args.path_regex)
        dxf_paths = [p for p in dxf_paths if rx.search(p)]

    if not dxf_paths:
        raise SystemExit(
            f"No DXF files matched. root_dir={ROOT_DIR} dxf_glob={DXF_GLOB} "
            f"path_regex={args.path_regex!r} dxf_file_count={len(args.dxf_file)}"
        )

    print(
        f"[run] files={len(dxf_paths)} root_dir={ROOT_DIR} "
        f"dxf_glob={DXF_GLOB!r} gap_tol={GAP_TOL} path_regex={args.path_regex!r}"
    )

    for p in dxf_paths:
        process_one_file(p)


if __name__ == "__main__":
    main()
