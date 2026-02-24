#!/usr/bin/env python3
import os, json, math, glob
import ezdxf
from ezdxf import edgesmith, edgeminer
from ezdxf.math import Vec3

# Point this at the directory that contains multiple slice_* subdirectories
ROOT_DIR = "/work/geometry/createGeometryFromCAD/DXF_slices_parts"
DXF_GLOB = "**/*.dxf"   # recursive
GAP_TOL  = 1e-6


def _vxy(p, nd=12):
    p = Vec3(p)
    return (round(p.x, nd), round(p.y, nd))

def _cross2(ax, ay, bx, by):
    return ax * by - ay * bx

def _arc_direction_ccw(x1, y1, x2, y2, cx, cy):
    return _cross2(x1 - cx, y1 - cy, x2 - cx, y2 - cy) > 0.0

def _knot_multiplicities(full_knots, tol=1e-12):
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

def _ellipse_params_from_dxf_ellipse(e):
    c = Vec3(e.dxf.center)
    maj = Vec3(e.dxf.major_axis)
    a = maj.magnitude
    b = a * float(e.dxf.ratio)
    phi_deg = math.degrees(math.atan2(maj.y, maj.x))
    return (float(c.x), float(c.y), float(a), float(b), float(phi_deg))

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

def _dist2_xy(a, b):
    return (a.x - b.x) ** 2 + (a.y - b.y) ** 2

def _reverse_bspline_payload_inplace(payload_dict):
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
            direction = _arc_direction_ccw(x1, y1, x2, y2, cx, cy)
            pts.append(["arc", x1, y1, cx, cy, direction])

        elif t == "ELLIPSE":
            cx, cy, a, b, phi_deg = _ellipse_params_from_dxf_ellipse(payload)
            inversed = not _arc_direction_ccw(x1, y1, x2, y2, cx, cy)
            pts.append(["ellipse", x1, y1, cx, cy, a, b, phi_deg, inversed])

        elif t == "SPLINE":
            degree = int(payload.dxf.degree)
            uniq_knots, mults = _knot_multiplicities(list(payload.knots))
            poles = [(float(Vec3(p).x), float(Vec3(p).y)) for p in payload.control_points]
            weights = [float(w) for w in getattr(payload, "weights", [])]
            flags = int(payload.dxf.flags)
            periodic = bool(flags & 2)
            if periodic:
                raise NotImplementedError("Periodic SPLINE not supported.")

            payload_dict = {
                "degree": degree,
                "poles": poles,
                "weights": weights,
                "knots": uniq_knots,
                "mults": mults,
                "periodic": periodic,
                "flags": flags,
            }

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

def dxf_loops_to_components(msp, gap_tol=GAP_TOL):
    entities = []
    for e in msp:
        if e.dxftype() == "INSERT":
            entities.extend(list(e.virtual_entities()))
        else:
            entities.append(e)

    edges = list(edgesmith.edges_from_entities_2d(entities))
    deposit = edgeminer.Deposit(edges, gap_tol=gap_tol)
    loops = list(edgeminer.find_all_loops(deposit))
    if not loops:
        loop = edgeminer.find_loop(deposit)
        if loop is None:
            return []
        loops = [loop]

    components = []
    for loop_edges in loops:
        chain = order_and_orient_closed_loop(loop_edges, gap_tol=gap_tol)
        components.append(_chain_to_component(chain))
    return components

def _indexed_dump_paths(base_json_path, count):
    base_json_path = os.path.abspath(base_json_path)
    root, ext = os.path.splitext(base_json_path)
    if ext.lower() != ".json":
        ext = ".json"
    if count <= 1:
        return [root + ext]
    return [f"{root}_{i}{ext}" for i in range(count)]

def _dxf_to_base_json_path(dxf_path):
    dxf_path = os.path.abspath(dxf_path)
    d = os.path.dirname(dxf_path)
    stem = os.path.splitext(os.path.basename(dxf_path))[0]
    return os.path.join(d, stem + ".json")

def process_one_file(dxf_path):
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()
    components = dxf_loops_to_components(msp, gap_tol=GAP_TOL)

    base_json = _dxf_to_base_json_path(dxf_path)
    os.makedirs(os.path.dirname(base_json), exist_ok=True)

    if not components:
        print(f"[SKIP] no loops: {dxf_path}")
        return

    paths = _indexed_dump_paths(base_json, len(components))
    for comp, path in zip(components, paths):
        with open(path, "w", encoding="utf-8") as f:
            json.dump(comp, f, indent=2, ensure_ascii=False)
        print(path)

def main():
    pattern = os.path.join(os.path.abspath(ROOT_DIR), DXF_GLOB)
    dxf_paths = sorted(glob.glob(pattern, recursive=True))
    if not dxf_paths:
        raise SystemExit(f"No DXF files found under: {ROOT_DIR}")

    for p in dxf_paths:
        try:
            process_one_file(p)
        except Exception as e:
            print(f"[ERROR] {p}: {e}")

if __name__ == "__main__":
    main()