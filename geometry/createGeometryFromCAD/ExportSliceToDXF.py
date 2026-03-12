import os
import re
import math
import argparse
import json
import csv

import FreeCAD
import Import
import Part

import ezdxf

import multiprocessing as mp


DEFAULT_STEP_PATH = "/work/geometry/createGeometryFromCAD/CAD_files/XENT-TPC_20250428.STEP"
DEFAULT_OUT_DIR = "/work/geometry/createGeometryFromCAD/DXF_slices_parts"
DEFAULT_WORKERS = 20
DEFAULT_TOL = 0.1


# -----------------------
# STEP load + object walk
# -----------------------

def load_step(step_path: str, doc_name: str = "imported_step"):
    doc = FreeCAD.newDocument(doc_name)
    Import.insert(step_path, doc.Name)
    doc.recompute()
    return doc

def find_step_root(doc):
    root = None
    for o in getattr(doc, "Objects", []) or []:
        if not getattr(o, "InList", []):
            root = o
            break
    if root is None:
        raise RuntimeError("Could not determine STEP root object")
    return root

def collect_leaf_geometry_entries(root):
    out = []

    def has_nonnull_shape(o) -> bool:
        try:
            sh = getattr(o, "Shape", None)
            return sh is not None and hasattr(sh, "isNull") and not sh.isNull()
        except Exception:
            return False

    def rec(parent, name_toks, label_toks, stack_ids):
        # Guard only against true recursion cycles. Do NOT globally deduplicate
        # objects by id: repeated assembly occurrences can legitimately reference
        # the same underlying object and still need separate exports.
        pid = id(parent)
        if pid in stack_ids:
            return
        stack_ids.add(pid)

        try:
            children = list(getattr(parent, "Group", []) or [])
        except Exception:
            children = []

        if not children:
            if parent is not root and has_nonnull_shape(parent):
                name = str(getattr(parent, "Name", "") or "")
                label = str(getattr(parent, "Label", "") or "")
                if name:
                    out.append(
                        {
                            "subpath": ".".join(name_toks + [name]) + ".",
                            "label_path": "/".join(label_toks + [label or name]),
                        }
                    )
            stack_ids.remove(pid)
            return

        for ch in children:
            name = str(getattr(ch, "Name", "") or "")
            label = str(getattr(ch, "Label", "") or "")
            if not name:
                continue

            rec(ch, name_toks + [name], label_toks + [label or name], stack_ids)

        stack_ids.remove(pid)

    root_label = str(getattr(root, "Label", "") or getattr(root, "Name", "") or "root")
    rec(root, [], [root_label], set())
    return out

def collect_leaf_geometry_subpaths(root):
    return [e["subpath"] for e in collect_leaf_geometry_entries(root)]

def _normalize_component_token(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (text or "").lower())

def _filter_leaf_subpaths(
    subpaths,
    *,
    component_match: str | None = None,
    component_regex: str | None = None,
    require_single_match: bool = False,
    search_aliases: dict | None = None,
):
    subpaths = list(subpaths or [])
    if not component_match and not component_regex:
        return subpaths

    rx = re.compile(component_regex) if component_regex else None
    token = _normalize_component_token(component_match) if component_match else None

    matched = []
    for sp in subpaths:
        alias = ""
        if search_aliases:
            alias = str(search_aliases.get(sp, "") or "")
        search_blob = f"{sp} {alias}".strip()
        sp_norm = _normalize_component_token(search_blob)
        if token and token not in sp_norm:
            continue
        if rx and not (rx.search(sp) or (alias and rx.search(alias))):
            continue
        matched.append(sp)

    if not matched:
        raise RuntimeError(
            f"No leaf components matched component_match={component_match!r} "
            f"component_regex={component_regex!r}"
        )

    if require_single_match and len(matched) != 1:
        preview = "\n".join(f"  - {s}" for s in matched[:25])
        suffix = "\n  - ..." if len(matched) > 25 else ""
        raise RuntimeError(
            f"Expected exactly one leaf component, but matched {len(matched)}.\n"
            f"Refine --component / --component-regex. Matches:\n{preview}{suffix}"
        )

    return matched

def _default_angles_with_offset(angle_offset_deg: float):
    return [15 * (i + 1) + float(angle_offset_deg) for i in range(0, 24)]

def _parse_cli_args():
    parser = argparse.ArgumentParser(
        description="Slice STEP geometry and export per-component DXF files."
    )
    parser.add_argument("--step-path", default=DEFAULT_STEP_PATH)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--tol",
        type=float,
        default=DEFAULT_TOL,
        help="Slicing/grouping tolerance in model units (default: 0.1).",
    )
    parser.add_argument(
        "--no-trim-u-ge-0",
        action="store_true",
        help="Disable half-plane trim and export full slice intersections.",
    )
    parser.add_argument(
        "--debug-artifacts-dir",
        default=None,
        help=(
            "Optional directory for per-component debug artifacts "
            "(DXF/SVG/FCStd/info) for unsplit and split-x slices."
        ),
    )
    parser.add_argument(
        "--debug-component",
        action="append",
        default=[],
        help=(
            "Component token for debug artifact export; may be passed multiple times. "
            "Matches normalized leaf subpath substrings."
        ),
    )
    parser.add_argument(
        "--component-audit-csv",
        default=None,
        help=(
            "Optional CSV path for per-component stage audit "
            "(selected, no-edges, exported mode, etc.)."
        ),
    )
    parser.add_argument(
        "--slice-fcstd-dir",
        default=None,
        help=(
            "Optional directory to dump final per-component sliced geometry as FCStd "
            "(includes final exported edges and trimmed/full candidates)."
        ),
    )

    parser.add_argument(
        "--angle-offset-deg",
        type=float,
        default=math.degrees(math.atan2(164.59, 1250.21)),
        help="Offset used by the default 24-angle schedule: 15*(i+1)+offset.",
    )
    parser.add_argument(
        "--single-angle-deg",
        type=float,
        default=None,
        help="Override to export only one explicit angle in degrees.",
    )
    parser.add_argument(
        "--first-angle-only",
        action="store_true",
        help="Override to export only the first scheduled angle (15 + offset).",
    )

    parser.add_argument(
        "--component",
        default=None,
        help=(
            "Normalized substring filter on leaf component subpaths. "
            "Can match one or many components."
        ),
    )
    parser.add_argument(
        "--component-regex",
        default=None,
        help=(
            "Regex filter on leaf component subpaths. "
            "Can be combined with --component."
        ),
    )
    parser.add_argument(
        "--require-single-component-match",
        action="store_true",
        help=(
            "Require exactly one filtered leaf component. "
            "Use this for strict single-component debug runs."
        ),
    )

    return parser.parse_args()

# -----------------------
# Geometry helpers
# -----------------------

def _slice_plane_normal(angle_deg: float, axis: str = "Y") -> FreeCAD.Vector:
    a = math.radians(float(angle_deg))
    ax = axis.upper().strip()
    if ax == "Y":
        n = FreeCAD.Vector(math.sin(a), 0.0, -math.cos(a))
    elif ax == "Z":
        n = FreeCAD.Vector(-math.sin(a), math.cos(a), 0.0)
    elif ax == "X":
        n = FreeCAD.Vector(0.0, -math.sin(a), math.cos(a))
    else:
        raise ValueError(f"axis must be 'X','Y','Z' (got {axis!r})")
    n.normalize()
    return n

def _plane_frame(angle_deg: float, axis: str = "Y"):
    a = math.radians(float(angle_deg))
    ax = axis.upper().strip()

    if ax == "X":
        A = FreeCAD.Vector(1, 0, 0)
        e1 = FreeCAD.Vector(0, 1, 0)
        e2 = FreeCAD.Vector(0, 0, 1)
    elif ax == "Y":
        A = FreeCAD.Vector(0, 1, 0)
        e1 = FreeCAD.Vector(1, 0, 0)
        e2 = FreeCAD.Vector(0, 0, 1)
    elif ax == "Z":
        A = FreeCAD.Vector(0, 0, 1)
        e1 = FreeCAD.Vector(1, 0, 0)
        e2 = FreeCAD.Vector(0, 1, 0)
    else:
        raise ValueError(f"axis must be 'X','Y','Z' (got {axis!r})")

    R = (math.cos(a) * e1) + (math.sin(a) * e2)
    try:
        A.normalize()
        R.normalize()
    except Exception:
        pass
    return A, R

def _to_2d(p: FreeCAD.Vector, origin: FreeCAD.Vector, A: FreeCAD.Vector, R: FreeCAD.Vector):
    d = p.sub(origin)
    return (float(d.dot(R)), float(d.dot(A)))

def _compute_global_bbox_diag(root):
    subpaths = collect_leaf_geometry_subpaths(root)
    bb = None
    for sp in subpaths:
        try:
            sh = root.getSubObject(sp)
        except Exception:
            continue
        if sh is None or not hasattr(sh, "isNull") or sh.isNull():
            continue
        try:
            b = sh.BoundBox
        except Exception:
            continue
        if bb is None:
            bb = FreeCAD.BoundBox(b.XMin, b.YMin, b.ZMin, b.XMax, b.YMax, b.ZMax)
        else:
            bb.add(b)
    if bb is None:
        return 1.0
    return math.sqrt(bb.XLength * bb.XLength + bb.YLength * bb.YLength + bb.ZLength * bb.ZLength) or 1.0

def make_keep_slab_u_ge_0(*, angle_deg: float, axis: str, origin, extent: float, thickness: float):
    origin3 = FreeCAD.Vector(*origin)
    A, R = _plane_frame(angle_deg=angle_deg, axis=axis)
    N = _slice_plane_normal(angle_deg=angle_deg, axis=axis)

    U = float(extent)
    V = float(extent)

    p1 = origin3 + R * 0.0 + A * (-V)
    p2 = origin3 + R * U   + A * (-V)
    p3 = origin3 + R * U   + A * ( V)
    p4 = origin3 + R * 0.0 + A * ( V)

    wire = Part.makePolygon([p1, p2, p3, p4, p1])
    face = Part.Face(wire)

    slab = face.extrude(N * thickness)
    slab.translate(N * (-0.5 * thickness))
    return slab

def make_keep_slab_x_ge_0(*, extent: float):
    e = float(extent)
    s = 2.0 * e
    return Part.makeBox(s, s, s, FreeCAD.Vector(0.0, -e, -e))

def make_keep_slab_x_lt_0(*, extent: float):
    e = float(extent)
    s = 2.0 * e
    return Part.makeBox(s, s, s, FreeCAD.Vector(-2.0 * e, -e, -e))

def _parse_debug_component_tokens(tokens):
    out = []
    seen = set()
    for raw in tokens or []:
        if not raw:
            continue
        for part in re.split(r"[,;\n]+", str(raw)):
            tok = _normalize_component_token(part.strip())
            if not tok or tok in seen:
                continue
            seen.add(tok)
            out.append(tok)
    return out

def _matches_debug_component(subpath: str, debug_tokens) -> bool:
    if not debug_tokens:
        return False
    s = _normalize_component_token(subpath)
    return any(t in s for t in debug_tokens)

def _shape_info(shape):
    info = {
        "solids": 0,
        "faces": 0,
        "shells": 0,
        "compounds": 0,
        "edges": 0,
    }
    if shape is None:
        return info
    for key in ("Solids", "Faces", "Shells", "Compounds", "Edges"):
        try:
            info[key.lower()] = len(list(getattr(shape, key)))
        except Exception:
            pass
    try:
        bb = shape.BoundBox
        info["bbox"] = {
            "xmin": float(bb.XMin),
            "xmax": float(bb.XMax),
            "ymin": float(bb.YMin),
            "ymax": float(bb.YMax),
            "zmin": float(bb.ZMin),
            "zmax": float(bb.ZMax),
            "xlen": float(bb.XLength),
            "ylen": float(bb.YLength),
            "zlen": float(bb.ZLength),
        }
    except Exception:
        pass
    return info

def _edge_polyline_2d_for_svg(edge, origin3, A, R):
    try:
        pts3 = list(edge.discretize(Number=80))
    except Exception:
        pts3 = []
    if not pts3:
        a, b = _edge_start_end_points(edge)
        if a is not None:
            pts3.append(a)
        if b is not None and (a is None or not _points_close(a, b, 1e-12)):
            pts3.append(b)
    return [_to_2d(p, origin3, A, R) for p in pts3]

def _write_edges_svg(edges, svg_path: str, *, angle_deg: float, axis: str, origin):
    os.makedirs(os.path.dirname(os.path.abspath(svg_path)), exist_ok=True)

    origin3 = FreeCAD.Vector(*origin)
    A, R = _plane_frame(angle_deg=angle_deg, axis=axis)
    polylines = []
    all_pts = []
    for e in edges or []:
        pts2 = _edge_polyline_2d_for_svg(e, origin3, A, R)
        if len(pts2) < 2:
            continue
        polylines.append(pts2)
        all_pts.extend(pts2)

    if not all_pts:
        with open(svg_path, "w", encoding="utf-8") as f:
            f.write(
                "<svg xmlns='http://www.w3.org/2000/svg' width='256' height='256' "
                "viewBox='0 0 256 256'>\n"
                "<rect x='0' y='0' width='256' height='256' fill='white'/>\n"
                "<text x='12' y='24' font-size='14' fill='black'>no edges</text>\n"
                "</svg>\n"
            )
        return

    xs = [p[0] for p in all_pts]
    ys = [p[1] for p in all_pts]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    span = max(xmax - xmin, ymax - ymin, 1.0)
    margin = 0.02 * span
    width = (xmax - xmin) + 2.0 * margin
    height = (ymax - ymin) + 2.0 * margin
    stroke_w = max(0.5, 0.002 * span)

    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(
            f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width:.9f} {height:.9f}'>\n"
        )
        f.write(f"<rect x='0' y='0' width='{width:.9f}' height='{height:.9f}' fill='white'/>\n")
        for poly in polylines:
            pts = []
            for x, y in poly:
                sx = (x - xmin) + margin
                sy = (ymax - y) + margin
                pts.append(f"{sx:.9f},{sy:.9f}")
            f.write(
                f"<polyline points='{' '.join(pts)}' fill='none' stroke='black' "
                f"stroke-width='{stroke_w:.9f}'/>\n"
            )
        f.write("</svg>\n")

def _save_component_slice_fcstd(
    fcstd_path: str,
    *,
    shape,
    unsplit_edges,
    split_x_ge_edges,
    split_x_lt_edges,
):
    doc_name = f"slice_debug_{abs(hash(fcstd_path)) % 1000000000}"
    try:
        if FreeCAD.getDocument(doc_name) is not None:
            FreeCAD.closeDocument(doc_name)
    except Exception:
        pass

    doc = FreeCAD.newDocument(doc_name)
    try:
        try:
            src = doc.addObject("Part::Feature", "SourceShape")
            src.Shape = shape
        except Exception:
            pass

        def add_edges(name, edges):
            if not edges:
                return
            try:
                obj = doc.addObject("Part::Feature", name)
                obj.Shape = Part.makeCompound(list(edges))
            except Exception:
                pass

        add_edges("SliceUnsplit", unsplit_edges)
        add_edges("SliceSplitXGe0", split_x_ge_edges)
        add_edges("SliceSplitXLt0", split_x_lt_edges)

        doc.recompute()
        os.makedirs(os.path.dirname(os.path.abspath(fcstd_path)), exist_ok=True)
        doc.saveAs(fcstd_path)
    finally:
        try:
            FreeCAD.closeDocument(doc.Name)
        except Exception:
            pass

def _save_export_slice_fcstd(
    fcstd_path: str,
    *,
    shape,
    final_edges,
    trimmed_edges=None,
    full_edges=None,
):
    """
    Save exact slice state used by exporter for GUI inspection.
    """
    doc_name = f"slice_export_{abs(hash(fcstd_path)) % 1000000000}"
    try:
        if FreeCAD.getDocument(doc_name) is not None:
            FreeCAD.closeDocument(doc_name)
    except Exception:
        pass

    doc = FreeCAD.newDocument(doc_name)
    try:
        try:
            src = doc.addObject("Part::Feature", "SourceShape")
            src.Shape = shape
        except Exception:
            pass

        def add_edges(name, edges):
            if not edges:
                return
            try:
                obj = doc.addObject("Part::Feature", name)
                obj.Shape = Part.makeCompound(list(edges))
            except Exception:
                pass

        add_edges("SliceFinal", final_edges)
        add_edges("SliceTrimmedCandidate", trimmed_edges)
        add_edges("SliceFullCandidate", full_edges)

        doc.recompute()
        os.makedirs(os.path.dirname(os.path.abspath(fcstd_path)), exist_ok=True)
        doc.saveAs(fcstd_path)
    finally:
        try:
            FreeCAD.closeDocument(doc.Name)
        except Exception:
            pass

def write_component_debug_artifacts(
    *,
    subpath: str,
    shape,
    angle_deg: float,
    axis: str,
    origin,
    plane_face,
    tol: float,
    per_solid_clean: bool,
    extent: float,
    debug_artifacts_dir: str,
):
    safe = _sanitize_layer(subpath, max_len=180)
    comp_dir = os.path.join(
        os.path.abspath(debug_artifacts_dir),
        f"slice_{float(angle_deg):06.2f}deg",
        safe,
    )
    os.makedirs(comp_dir, exist_ok=True)

    slab_x_ge = make_keep_slab_x_ge_0(extent=extent)
    slab_x_lt = make_keep_slab_x_lt_0(extent=extent)

    unsplit_edges = slice_component_edges(
        shape, plane_face=plane_face, keep_slab=None, tol=tol, per_solid_clean=per_solid_clean
    )
    split_x_ge_edges = slice_component_edges(
        shape, plane_face=plane_face, keep_slab=slab_x_ge, tol=tol, per_solid_clean=per_solid_clean
    )
    split_x_lt_edges = slice_component_edges(
        shape, plane_face=plane_face, keep_slab=slab_x_lt, tol=tol, per_solid_clean=per_solid_clean
    )

    cases = [
        ("unsplit", unsplit_edges),
        ("split_x_ge0", split_x_ge_edges),
        ("split_x_lt0", split_x_lt_edges),
    ]

    case_info = {}
    for label, edges in cases:
        dxf_path = os.path.join(comp_dir, f"{label}.dxf")
        svg_path = os.path.join(comp_dir, f"{label}.svg")
        mode, prim, gsum = export_component_wiregrouped_or_raw(
            subpath=f"{subpath}__{label}",
            edges=edges,
            dxf_path=dxf_path,
            angle_deg=angle_deg,
            axis=axis,
            origin=origin,
            tol=tol,
        )
        _write_edges_svg(edges, svg_path, angle_deg=angle_deg, axis=axis, origin=origin)
        case_info[label] = {
            "mode": mode,
            "edge_count": len(edges),
            "group_summary": gsum,
            "primitive_counts": prim,
            "dxf_path": dxf_path,
            "svg_path": svg_path,
        }

    fcstd_path = os.path.join(comp_dir, "slice_debug.FCStd")
    _save_component_slice_fcstd(
        fcstd_path,
        shape=shape,
        unsplit_edges=unsplit_edges,
        split_x_ge_edges=split_x_ge_edges,
        split_x_lt_edges=split_x_lt_edges,
    )

    info = {
        "subpath": subpath,
        "component_safe_name": safe,
        "angle_deg": float(angle_deg),
        "axis": axis,
        "origin": [float(origin[0]), float(origin[1]), float(origin[2])],
        "tolerance": float(tol),
        "shape_info": _shape_info(shape),
        "slice_cases": case_info,
        "fcstd_path": fcstd_path,
    }
    info_path = os.path.join(comp_dir, "info.json")
    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    return comp_dir


# -----------------------
# Slice per component (NO merging that can drop duplicates)
# -----------------------

def slice_component_edges(
    shape,
    *,
    plane_face,
    keep_slab=None,
    tol: float = DEFAULT_TOL,
    per_solid_clean: bool = True,
):
    tol_eff = float(max(tol, 1e-6))

    try:
        solids = list(shape.Solids)
    except Exception:
        solids = []
    # Some STEP leaves are shells/faces (no solids). Those must still be sliced.
    slicing_units = solids if solids else [shape]

    out_edges = []
    for unit in slicing_units:
        unit_trim = unit
        if keep_slab is not None:
            try:
                unit_trim = unit.common(keep_slab)
            except Exception:
                continue
            if unit_trim is None or unit_trim.isNull():
                continue

        try:
            sec = unit_trim.section(plane_face)
        except Exception:
            continue
        if sec is None or sec.isNull():
            continue

        if per_solid_clean:
            try:
                sec = sec.cleanTolerance(tol_eff)
            except Exception:
                pass

        try:
            out_edges.extend(list(sec.Edges))
        except Exception:
            pass

    return out_edges

# -----------------------
# Wire grouping diagnostics
# -----------------------

def wire_group_edges(edges, tol: float):
    """
    Returns:
      groups: list[list[Edge]] (from Part.sortEdges; fallback single group)
      status: list[str] per group: CLOSED / OPEN / WIRE_FAIL
      ok: bool -> True iff no group is WIRE_FAIL (i.e., wire-grouping succeeded)
    """
    tol_eff = float(max(tol, 1e-6))
    edges = list(edges) if edges else []
    if not edges:
        return [], [], True

    try:
        groups = [list(g) for g in Part.sortEdges(edges)]
    except Exception:
        groups = [edges]

    status = []
    ok = True
    for g in groups:
        try:
            w = Part.Wire(g)
        except Exception:
            status.append("WIRE_FAIL")
            ok = False
            continue

        try:
            w.fixTolerance(tol_eff)
        except Exception:
            pass

        try:
            status.append("CLOSED" if w.isClosed() else "OPEN")
        except Exception:
            status.append("OPEN")

    return groups, status, ok


# -----------------------
# DXF export (exact primitives; SPLINE fallback)
# -----------------------

def _sanitize_layer(name: str, max_len: int = 200) -> str:
    name = name.strip()
    name = re.sub(r'[<>/\\":;?*\|\=]', "_", name)
    name = name.replace("\n", "_").replace("\r", "_").replace("\t", "_")
    return (name or "UNNAMED")[:max_len]

def _edge_is_full_circle(edge, tol=1e-7):
    try:
        u0, u1 = edge.ParameterRange
        span = abs(float(u1 - u0))
    except Exception:
        return False

    try:
        vs = getattr(edge, "Vertexes", None)
        if not vs:
            return False
        if len(vs) == 1:
            closed = True
        else:
            p0 = vs[0].Point
            p1 = vs[-1].Point
            closed = (p0.distanceToPoint(p1) <= tol)
        if not closed:
            return False
    except Exception:
        return False

    two_pi = 2.0 * math.pi
    span_mod = span % two_pi
    eps = 1e-6
    return (span_mod <= eps) or (abs(span_mod - two_pi) <= eps)

def _angle_deg_in_plane(pt2, center2):
    return math.degrees(math.atan2(pt2[1] - center2[1], pt2[0] - center2[0]))

def _arc_is_ccw_in_export_coords(edge, circle_curve, center2, origin3, A, R):
    try:
        u0, u1 = edge.ParameterRange
        span = float(u1 - u0)
        if span == 0.0:
            return True
        du = span * 1e-6
        if du == 0.0:
            du = 1e-9
        try:
            p0 = circle_curve.value(u0)
            p1s = circle_curve.value(u0 + du)
        except Exception:
            p0 = circle_curve.value(u1)
            p1s = circle_curve.value(u1 - du)

        p0_2 = _to_2d(p0, origin3, A, R)
        p1_2 = _to_2d(p1s, origin3, A, R)

        x0, y0 = (p0_2[0] - center2[0], p0_2[1] - center2[1])
        x1, y1 = (p1_2[0] - center2[0], p1_2[1] - center2[1])
        return (x0 * y1 - y0 * x1) > 0.0
    except Exception:
        return True

def _add_spline_from_edge(msp, edge, layer, origin3, A, R, flipped: bool = False):
    bs = None
    last_err = None

    try:
        c = edge.Curve
        if hasattr(c, "toBSpline"):
            bs = c.toBSpline()
    except Exception as ex:
        last_err = ex
        bs = None

    if bs is None:
        try:
            nurbs_shape = edge.toNurbs()
            if nurbs_shape and hasattr(nurbs_shape, "Edges") and nurbs_shape.Edges:
                bs = nurbs_shape.Edges[0].Curve
        except Exception as ex:
            last_err = ex
            bs = None

    if bs is None:
        try:
            tmp = Part.makeCompound([edge])
            nurbs_tmp = tmp.toNurbs()
            if nurbs_tmp and hasattr(nurbs_tmp, "Edges") and nurbs_tmp.Edges:
                bs = nurbs_tmp.Edges[0].Curve
        except Exception as ex:
            last_err = ex
            bs = None

    if bs is None:
        raise RuntimeError(f"Could not convert edge to NURBS for DXF SPLINE export: {last_err}")

    try:
        u0, u1 = edge.ParameterRange
        bs_seg = bs.copy()
        bs_seg.segment(u0, u1)
        bs = bs_seg
    except Exception:
        pass

    degree = int(bs.Degree)
    poles3 = list(bs.getPoles())

    # Preserve topological semantics for downstream contour construction.
    # Do not infer closure from the base curve: an open edge can lie on a closed/periodic curve.
    try:
        edge_closed = bool(edge.isClosed())
    except Exception:
        edge_closed = False
    edge_periodic = False
    if edge_closed:
        try:
            c0 = edge.Curve
            if hasattr(c0, "isPeriodic"):
                edge_periodic = bool(c0.isPeriodic())
        except Exception:
            pass
        try:
            edge_periodic = edge_periodic or bool(bs.isPeriodic())
        except Exception:
            pass

    # If an open edge converted to a closed-loop pole set, rebuild as an open spline
    # from sampled edge points to preserve contour topology.
    p0, p1 = _edge_start_end_points(edge)
    edge_open_by_vertices = (
        (p0 is not None) and (p1 is not None) and (not _points_close(p0, p1, 1e-7))
    )
    if edge_open_by_vertices and len(poles3) >= 2 and _points_close(poles3[0], poles3[-1], 1e-7):
        samples = list(edge.discretize(Number=max(20, degree * 4)))
        if len(samples) >= 2:
            if _points_close(samples[0], samples[-1], 1e-7):
                samples = samples[:-1]
            if len(samples) >= 2:
                try:
                    bs_open = Part.BSplineCurve()
                    bs_open.interpolate(samples)
                    bs = bs_open
                    degree = int(bs.Degree)
                    poles3 = list(bs.getPoles())
                except Exception:
                    pass

    poles2 = [_to_2d(p, origin3, A, R) for p in poles3]

    knots = list(bs.getKnots())
    mults = list(bs.getMultiplicities())
    knot_vector = []
    for kv, m in zip(knots, mults):
        knot_vector.extend([float(kv)] * int(m))

    weights = None
    try:
        w = list(bs.getWeights())
        if any(abs(float(x) - 1.0) > 1e-12 for x in w):
            weights = [float(x) for x in w]
    except Exception:
        weights = None

    if flipped:
        poles2 = list(reversed(poles2))
        if weights is not None:
            weights = list(reversed(weights))
        if knot_vector:
            kmin = float(knot_vector[0])
            kmax = float(knot_vector[-1])
            knot_vector = [float((kmin + kmax) - k) for k in reversed(knot_vector)]

    spline = msp.add_spline(dxfattribs={"layer": layer})
    spline.dxf.degree = degree

    ok = False
    for cp_attr in ("control_points", "controlpoints", "points"):
        try:
            setattr(spline, cp_attr, poles2)
            ok = True
            break
        except Exception:
            pass
    if not ok:
        spline.set_control_points(poles2)

    ok = False
    for k_attr in ("knots", "knot_values"):
        try:
            setattr(spline, k_attr, knot_vector)
            ok = True
            break
        except Exception:
            pass
    if not ok:
        spline.set_knots(knot_vector)

    if weights is not None:
        try:
            spline.weights = weights
        except Exception:
            spline.set_weights(weights)

    flags = 0
    if edge_closed:
        flags |= 1  # closed
    if edge_periodic:
        flags |= 2  # periodic
    if weights is not None:
        flags |= 4  # rational
    flags |= 8      # planar in exported 2D slice
    try:
        spline.dxf.flags = flags
    except Exception:
        pass

def _add_ellipse_from_edge(msp, edge, layer, origin3, A, R, flipped: bool = False):
    c = edge.Curve
    required = ("Center", "MajorRadius", "MinorRadius", "MajorAxis")
    if not all(hasattr(c, k) for k in required):
        return False

    center3 = c.Center
    major_axis_dir3 = c.MajorAxis
    try:
        major_axis_dir3.normalize()
    except Exception:
        pass

    major_len = float(c.MajorRadius)
    minor_len = float(c.MinorRadius)
    if major_len <= 0 or minor_len <= 0:
        return False

    center2 = _to_2d(center3, origin3, A, R)

    maj_end3 = center3.add(major_axis_dir3.multiply(major_len))
    maj_end2 = _to_2d(maj_end3, origin3, A, R)
    major_vec2 = (maj_end2[0] - center2[0], maj_end2[1] - center2[1])

    ratio = float(minor_len / major_len)

    u0, u1 = edge.ParameterRange
    start_param = float(u0)
    end_param = float(u1)
    if flipped:
        start_param, end_param = end_param, start_param

    msp.add_ellipse(
        center=center2,
        major_axis=major_vec2,
        ratio=ratio,
        start_param=start_param,
        end_param=end_param,
        dxfattribs={"layer": layer},
    )
    return True

def _points_close(p0, p1, tol: float) -> bool:
    try:
        return p0.distanceToPoint(p1) <= tol
    except Exception:
        return False

def _point_distance(p0, p1) -> float:
    try:
        return p0.distanceToPoint(p1)
    except Exception:
        return float("inf")

def _unpack_edge_item(item):
    """
    Normalize edge item encoding.
    Supports raw Edge and nested (item, flipped_bool) tuples.
    Returns (edge_or_none, flipped_bool).
    """
    flipped = False
    e = item
    depth = 0
    while isinstance(e, tuple) and len(e) == 2 and isinstance(e[1], bool) and depth < 8:
        e, f = e
        flipped = (flipped != bool(f))
        depth += 1
    return e, flipped

def _edge_start_end_points(edge):
    vs = list(getattr(edge, "Vertexes", []) or [])
    if not vs:
        return None, None
    if len(vs) == 1:
        return vs[0].Point, vs[0].Point
    return vs[0].Point, vs[-1].Point

def _reverse_edge(edge):
    try:
        return edge.reversed()
    except Exception:
        try:
            out = edge.copy()
            out.reverse()
            return out
        except Exception:
            return edge

def _order_closed_group_edges(edges, tol: float):
    """
    Orders one closed group into a head-to-tail chain with explicit orientation.
    Returns list of (edge, flipped), where flipped=True means traverse end->start.
    """
    edges = list(edges or [])
    if not edges:
        return []
    if len(edges) == 1:
        return [(edges[0], False)]

    def _try_with_first_flip(first_flip: bool):
        remaining = list(edges)
        first = remaining.pop(0)
        a0, b0 = _edge_start_end_points(first)
        if a0 is None or b0 is None:
            return None

        if first_flip:
            start_start = b0
            cur_end = a0
        else:
            start_start = a0
            cur_end = b0

        ordered = [(first, first_flip)]

        while remaining:
            found = False
            for i, cand in enumerate(remaining):
                a, b = _edge_start_end_points(cand)
                if a is None or b is None:
                    continue
                if _points_close(cur_end, a, tol):
                    ordered.append((cand, False))
                    cur_end = b
                    remaining.pop(i)
                    found = True
                    break
                if _points_close(cur_end, b, tol):
                    ordered.append((cand, True))
                    cur_end = a
                    remaining.pop(i)
                    found = True
                    break
            if not found:
                return None

        if not _points_close(cur_end, start_start, max(tol, 1e-6)):
            return None
        return ordered

    chain = _try_with_first_flip(False)
    if chain is None:
        chain = _try_with_first_flip(True)
    if chain is None:
        return [(e, False) for e in edges]
    return chain

def _order_open_group_edges(edges, tol: float):
    """
    Greedily order an open group into a head-to-tail chain.
    Returns list of (edge, flipped).
    """
    raw_items = list(edges or [])
    if not raw_items:
        return []

    remaining = []
    for it in raw_items:
        e, f = _unpack_edge_item(it)
        if e is None:
            continue
        remaining.append((e, f))
    if not remaining:
        return []

    first_e, first_f = remaining.pop(0)
    ordered = [(first_e, first_f)]

    a0, b0 = _edge_start_end_points(first_e)
    if a0 is None or b0 is None:
        ordered.extend((e, f) for (e, f) in remaining)
        return ordered
    if first_f:
        a0, b0 = b0, a0
    cur_end = b0

    while remaining:
        best_i = None
        best_flip = False
        best_next_end = None
        best_d = None

        for i, (cand, cand_flip) in enumerate(remaining):
            a, b = _edge_start_end_points(cand)
            if a is None or b is None:
                continue
            if cand_flip:
                a, b = b, a

            da = _point_distance(cur_end, a)
            db = _point_distance(cur_end, b)

            if best_d is None or da < best_d:
                best_d = da
                best_i = i
                best_flip = cand_flip
                best_next_end = b

            if best_d is None or db < best_d:
                best_d = db
                best_i = i
                best_flip = (not cand_flip)
                best_next_end = a

        if best_i is None:
            break

        e, _f = remaining.pop(best_i)
        ordered.append((e, best_flip))
        cur_end = best_next_end

    # If ordering could not consume all edges, append the rest in original orientation.
    ordered.extend((e, f) for (e, f) in remaining)
    return ordered

def _edge_polyline_2d(edge, origin3, A, R, min_samples=20, flipped: bool = False):
    """
    Returns a sampled 2D polyline for one edge in export coordinates.
    """
    pts3 = []
    try:
        n = int(max(8, min_samples))
        pts3 = list(edge.discretize(Number=n))
    except Exception:
        pts3 = []

    if not pts3:
        a, b = _edge_start_end_points(edge)
        if a is not None:
            pts3.append(a)
        if b is not None and (a is None or not _points_close(a, b, 1e-12)):
            pts3.append(b)

    pts2 = [_to_2d(p, origin3, A, R) for p in pts3]
    if flipped:
        pts2 = list(reversed(pts2))
    return pts2

def _signed_area_xy(poly):
    """
    Shoelace signed area. Positive => CCW.
    """
    if len(poly) < 3:
        return 0.0
    s = 0.0
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        s += x1 * y2 - x2 * y1
    return 0.5 * s

def _point_in_poly_xy(pt, poly):
    """
    Ray-casting point-in-polygon test (boundary-inclusive best effort).
    """
    if len(poly) < 3:
        return False
    x, y = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            x_cross = x1 + (y - y1) * (x2 - x1) / ((y2 - y1) if (y2 - y1) != 0 else 1e-30)
            if x <= x_cross:
                inside = not inside
    return inside

def _poly_centroid_xy(poly):
    """
    Area-weighted centroid; falls back to arithmetic mean if near-degenerate.
    """
    if len(poly) < 3:
        if not poly:
            return (0.0, 0.0)
        sx = sum(p[0] for p in poly)
        sy = sum(p[1] for p in poly)
        return (sx / len(poly), sy / len(poly))

    a2 = 0.0
    cx = 0.0
    cy = 0.0
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        cr = x0 * y1 - x1 * y0
        a2 += cr
        cx += (x0 + x1) * cr
        cy += (y0 + y1) * cr

    if abs(a2) < 1e-20:
        sx = sum(p[0] for p in poly)
        sy = sum(p[1] for p in poly)
        return (sx / len(poly), sy / len(poly))

    inv = 1.0 / (3.0 * a2)
    return (cx * inv, cy * inv)

def _ring_probe_points_xy(ring):
    """
    Build a few stable probe points that should lie on/just inside the ring.
    """
    if not ring:
        return []

    n = len(ring)
    idxs = sorted(set([0, n // 3, (2 * n) // 3]))
    c = _poly_centroid_xy(ring)
    probes = []
    for idx in idxs:
        x, y = ring[idx]
        # tiny move toward centroid; helps avoid exact-boundary ambiguity
        px = x + 1e-9 * (c[0] - x)
        py = y + 1e-9 * (c[1] - y)
        probes.append((px, py))
    return probes

def _ring_inside_ring(inner_ring, outer_ring, inner_area: float, outer_area: float) -> bool:
    """
    Heuristic containment check robust to concentric shells:
    use area ordering + probe points from the candidate inner ring.
    """
    if len(inner_ring) < 3 or len(outer_ring) < 3:
        return False
    if abs(inner_area) >= abs(outer_area):
        return False

    for p in _ring_probe_points_xy(inner_ring):
        if _point_in_poly_xy(p, outer_ring):
            return True
    return False

def _build_group_ring_2d(edges, angle_deg: float, axis: str, origin):
    """
    Sample a closed edge chain into one 2D ring in export coordinates.
    """
    origin3 = FreeCAD.Vector(*origin)
    A, R = _plane_frame(angle_deg=angle_deg, axis=axis)

    ring = []
    for i, item in enumerate(edges):
        if isinstance(item, tuple) and len(item) == 2:
            e, flipped = item
        else:
            e, flipped = item, False
        seg = _edge_polyline_2d(e, origin3, A, R, min_samples=20, flipped=flipped)
        if not seg:
            continue
        if i > 0 and ring:
            if abs(seg[0][0] - ring[-1][0]) < 1e-9 and abs(seg[0][1] - ring[-1][1]) < 1e-9:
                seg = seg[1:]
        ring.extend(seg)

    # remove adjacent duplicates
    dedup = []
    for p in ring:
        if not dedup:
            dedup.append(p)
            continue
        if abs(p[0] - dedup[-1][0]) < 1e-9 and abs(p[1] - dedup[-1][1]) < 1e-9:
            continue
        dedup.append(p)

    if len(dedup) >= 2:
        if abs(dedup[0][0] - dedup[-1][0]) < 1e-9 and abs(dedup[0][1] - dedup[-1][1]) < 1e-9:
            dedup.pop()
    return dedup

def _orient_closed_groups_for_cutouts(groups, status, angle_deg: float, axis: str, origin, tol: float):
    """
    Enforce nested winding:
      - even containment depth (outer shells): CCW
      - odd containment depth (holes): CW
    Works on closed groups only; open/wire-fail groups are passed through.
    """
    oriented = [[(e, False) for e in g] for g in groups]

    closed_infos = []
    for idx, (g, st) in enumerate(zip(groups, status)):
        if st != "CLOSED":
            continue
        chain = _order_closed_group_edges(g, tol=max(tol, 1e-6))
        ring = _build_group_ring_2d(chain, angle_deg=angle_deg, axis=axis, origin=origin)
        if len(ring) < 3:
            oriented[idx] = chain
            continue
        area = _signed_area_xy(ring)
        rep = _poly_centroid_xy(ring)
        closed_infos.append({
            "idx": idx,
            "edges": chain,
            "ring": ring,
            "area": area,
            "rep": rep,
        })

    if not closed_infos:
        return oriented

    for i, info in enumerate(closed_infos):
        depth = 0
        for j, other in enumerate(closed_infos):
            if i == j:
                continue
            if _ring_inside_ring(
                info["ring"],
                other["ring"],
                info["area"],
                other["area"],
            ):
                depth += 1
        desired_ccw = (depth % 2 == 0)
        current_ccw = info["area"] > 0.0
        edges = info["edges"]
        if current_ccw != desired_ccw:
            edges = [(e, not flipped) for (e, flipped) in reversed(edges)]
        oriented[info["idx"]] = edges

    return oriented

def _export_edges_to_msp(
    *,
    msp,
    layer: str,
    edges,
    angle_deg: float,
    axis: str,
    origin,
    tol: float,
    bridge_open_gaps: bool = False,
    force_close_open: bool = False,
    flatten_open_splines: bool = False,
):
    origin3 = FreeCAD.Vector(*origin)
    A, R = _plane_frame(angle_deg=angle_deg, axis=axis)
    tol_eff = float(max(tol, 1e-6))

    prim_counts = {"LINE": 0, "ARC": 0, "CIRCLE": 0, "ELLIPSE": 0, "SPLINE": 0, "FAIL": 0}
    first_start2 = None
    prev_end2 = None

    for item in edges:
        e, flipped = _unpack_edge_item(item)
        if e is None or not hasattr(e, "Curve"):
            prim_counts["FAIL"] += 1
            continue

        start3, end3 = _edge_start_end_points(e)
        if start3 is not None and end3 is not None and flipped:
            start3, end3 = end3, start3
        start2 = _to_2d(start3, origin3, A, R) if start3 is not None else None
        end2 = _to_2d(end3, origin3, A, R) if end3 is not None else None

        if first_start2 is None and start2 is not None:
            first_start2 = start2

        if bridge_open_gaps and prev_end2 is not None and start2 is not None:
            if math.hypot(prev_end2[0] - start2[0], prev_end2[1] - start2[1]) > tol_eff:
                msp.add_line(prev_end2, start2, dxfattribs={"layer": layer})
                prim_counts["LINE"] += 1

        c = e.Curve
        cname = c.__class__.__name__
        try:
            if cname in ("Line", "GeomLine"):
                if start2 is None or end2 is None:
                    p0 = e.Vertexes[0].Point
                    p1 = e.Vertexes[-1].Point
                    if flipped:
                        p0, p1 = p1, p0
                    start2 = _to_2d(p0, origin3, A, R)
                    end2 = _to_2d(p1, origin3, A, R)
                msp.add_line(start2, end2, dxfattribs={"layer": layer})
                prim_counts["LINE"] += 1
                prev_end2 = end2
                continue

            if hasattr(c, "Center") and hasattr(c, "Radius"):
                center2 = _to_2d(c.Center, origin3, A, R)
                r = float(c.Radius)

                if _edge_is_full_circle(e, tol=max(tol, 1e-6)):
                    msp.add_circle(center2, r, dxfattribs={"layer": layer})
                    prim_counts["CIRCLE"] += 1
                    continue

                u0, u1 = e.ParameterRange
                p0 = c.value(u0)
                p1 = c.value(u1)
                if flipped:
                    p0, p1 = p1, p0
                p0_2 = start2 if start2 is not None else _to_2d(p0, origin3, A, R)
                p1_2 = end2 if end2 is not None else _to_2d(p1, origin3, A, R)

                a1 = _angle_deg_in_plane(p0_2, center2)
                a2 = _angle_deg_in_plane(p1_2, center2)

                ccw = _arc_is_ccw_in_export_coords(e, c, center2, origin3, A, R)
                if flipped:
                    ccw = not ccw
                if not ccw:
                    a1, a2 = a2, a1

                msp.add_arc(center2, r, a1, a2, dxfattribs={"layer": layer})
                prim_counts["ARC"] += 1
                prev_end2 = p1_2
                continue

            if cname in ("Ellipse", "GeomEllipse") or (
                hasattr(c, "MajorRadius") and hasattr(c, "MinorRadius") and hasattr(c, "MajorAxis")
            ):
                if _add_ellipse_from_edge(msp, e, layer, origin3, A, R, flipped=flipped):
                    prim_counts["ELLIPSE"] += 1
                    prev_end2 = end2
                    continue

            if flatten_open_splines:
                poly2 = _edge_polyline_2d(e, origin3, A, R, min_samples=80, flipped=flipped)
                if len(poly2) >= 2:
                    for a2, b2 in zip(poly2, poly2[1:]):
                        msp.add_line(a2, b2, dxfattribs={"layer": layer})
                    prim_counts["LINE"] += max(0, len(poly2) - 1)
                    prev_end2 = poly2[-1]
                    continue

            _add_spline_from_edge(msp, e, layer, origin3, A, R, flipped=flipped)
            prim_counts["SPLINE"] += 1
            prev_end2 = end2

        except Exception:
            prim_counts["FAIL"] += 1
            prev_end2 = end2

    if force_close_open and first_start2 is not None and prev_end2 is not None:
        if math.hypot(prev_end2[0] - first_start2[0], prev_end2[1] - first_start2[1]) > tol_eff:
            msp.add_line(prev_end2, first_start2, dxfattribs={"layer": layer})
            prim_counts["LINE"] += 1

    return prim_counts


def export_component_wiregrouped_or_raw(
    *,
    subpath: str,
    edges,
    dxf_path: str,
    angle_deg: float,
    axis: str,
    origin,
    tol: float,
):
    """
    If wire grouping succeeds (no WIRE_FAIL groups), export grouped on multiple layers (normal directory).
    Else export all raw edges on a single layer (raw directory).
    Returns (mode, prim_counts, group_summary)
    """
    tol_eff = float(max(tol, 1e-6))
    safe = _sanitize_layer(subpath, max_len=180)

    groups, status, ok = wire_group_edges(edges, tol=tol_eff)

    dxf = ezdxf.new("R2010")
    msp = dxf.modelspace()

    group_summary = {"groups": len(groups), "closed": 0, "open": 0, "wire_fail": 0}
    for st in status:
        if st == "CLOSED":
            group_summary["closed"] += 1
        elif st == "OPEN":
            group_summary["open"] += 1
        else:
            group_summary["wire_fail"] += 1

    if ok:
        groups = _orient_closed_groups_for_cutouts(
            groups=groups,
            status=status,
            angle_deg=angle_deg,
            axis=axis,
            origin=origin,
            tol=tol_eff,
        )

        # wire grouped: one layer per group (keeps visibility of where the break is)
        mode = "WIRE_GROUPED"
        prim_total = {"LINE": 0, "ARC": 0, "CIRCLE": 0, "ELLIPSE": 0, "SPLINE": 0, "FAIL": 0}

        for i, (g, st) in enumerate(zip(groups, status)):
            layer = _sanitize_layer(f"{safe}__G{i:03d}_{st}", max_len=200)
            if layer not in dxf.layers:
                dxf.layers.new(layer)

            export_edges = g
            bridge_open_gaps = False
            force_close_open = False
            flatten_open_splines = False
            if st == "OPEN":
                export_edges = _order_open_group_edges(g, tol=tol_eff)
                # Keep OPEN groups exactly as open chains (no synthetic bridges/closures).
                bridge_open_gaps = False
                force_close_open = False
                flatten_open_splines = False

            prim = _export_edges_to_msp(
                msp=msp,
                layer=layer,
                edges=export_edges,
                angle_deg=angle_deg,
                axis=axis,
                origin=origin,
                tol=tol,
                bridge_open_gaps=bridge_open_gaps,
                force_close_open=force_close_open,
                flatten_open_splines=flatten_open_splines,
            )
            for k in prim_total:
                prim_total[k] += prim.get(k, 0)

        dxf.saveas(dxf_path)
        return mode, prim_total, group_summary

    # raw: put everything on one layer for maximum forensic visibility
    mode = "RAW"
    layer = safe
    if layer not in dxf.layers:
        dxf.layers.new(layer)

    prim = _export_edges_to_msp(
        msp=msp,
        layer=layer,
        edges=edges,
        angle_deg=angle_deg,
        axis=axis,
        origin=origin,
        tol=tol,
    )
    dxf.saveas(dxf_path)
    return mode, prim, group_summary


# -----------------------
# Export: angle subdir + one DXF per component
#   - normal dir: wire-grouped outputs (if wire grouping succeeded)
#   - out_dir_raw: raw outputs ONLY for components that fail wire-grouping
# -----------------------

def export_slices_per_component(
    *,
    root,
    out_dir: str,
    angles_deg,
    axis: str = "Y",
    origin=(0.0, 0.0, 0.0),
    tol: float = DEFAULT_TOL,
    trim_u_ge_0: bool = True,
    per_solid_clean: bool = True,
    component_match: str | None = None,
    component_regex: str | None = None,
    require_single_component_match: bool = False,
    debug_artifacts_dir: str | None = None,
    debug_component_tokens=None,
    component_audit_csv: str | None = None,
    slice_fcstd_dir: str | None = None,
):
    out_dir = os.path.abspath(out_dir)
    out_dir_raw = out_dir + "_raw"  # requested: _raw at highest-order directory

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(out_dir_raw, exist_ok=True)

    entries = collect_leaf_geometry_entries(root)
    subpaths_all = [e["subpath"] for e in entries]
    alias_by_subpath = {e["subpath"]: e.get("label_path", "") for e in entries}

    subpaths = _filter_leaf_subpaths(
        subpaths_all,
        component_match=component_match,
        component_regex=component_regex,
        require_single_match=require_single_component_match,
        search_aliases=alias_by_subpath,
    )
    if component_match or component_regex:
        print(
            f"[filter] selected_leaf_components={len(subpaths)} "
            f"component_match={component_match!r} component_regex={component_regex!r}"
        )
    debug_tokens = _parse_debug_component_tokens(debug_component_tokens)
    if debug_artifacts_dir and debug_tokens:
        print(
            f"[debug] artifact_dir={os.path.abspath(debug_artifacts_dir)} "
            f"debug_component_tokens={debug_tokens}"
        )
    elif debug_artifacts_dir:
        print(
            f"[debug] artifact_dir={os.path.abspath(debug_artifacts_dir)} "
            "debug_component_tokens=<all selected components>"
        )
    if slice_fcstd_dir:
        print(f"[debug] slice_fcstd_dir={os.path.abspath(slice_fcstd_dir)}")

    diag = _compute_global_bbox_diag(root)
    extent = diag * 3.0
    thickness = diag * 0.25

    total_with_edges = 0
    failed_wiregroup = 0
    audit_rows = []

    for ang in angles_deg:
        ang = float(ang)

        angle_dir = os.path.join(out_dir, f"slice_{ang:06.2f}deg")
        angle_dir_raw = os.path.join(out_dir_raw, f"slice_{ang:06.2f}deg")
        os.makedirs(angle_dir, exist_ok=True)
        os.makedirs(angle_dir_raw, exist_ok=True)

        N = _slice_plane_normal(ang, axis=axis)
        ox, oy, oz = origin
        plane_face = Part.Plane(FreeCAD.Vector(ox, oy, oz), N).toShape()

        keep_slab = None
        if trim_u_ge_0:
            keep_slab = make_keep_slab_u_ge_0(
                angle_deg=ang, axis=axis, origin=origin, extent=extent, thickness=thickness
            )

        print(f"[angle] {ang:.2f}° → {angle_dir}  (raw failures → {angle_dir_raw})")

        for sp in subpaths:
            label_path = alias_by_subpath.get(sp, "")
            try:
                sh = root.getSubObject(sp)
            except Exception:
                audit_rows.append(
                    {
                        "angle_deg": f"{ang:.8f}",
                        "subpath": sp,
                        "label_path": label_path,
                        "status": "GET_SUBOBJECT_FAILED",
                        "edge_source": "",
                        "edge_count": "0",
                        "mode": "",
                        "dxf_path": "",
                    }
                )
                continue
            if sh is None or not hasattr(sh, "isNull") or sh.isNull():
                audit_rows.append(
                    {
                        "angle_deg": f"{ang:.8f}",
                        "subpath": sp,
                        "label_path": label_path,
                        "status": "NULL_SHAPE",
                        "edge_source": "",
                        "edge_count": "0",
                        "mode": "",
                        "dxf_path": "",
                    }
                )
                continue

            debug_selected = (not debug_tokens) or _matches_debug_component(sp, debug_tokens)
            if debug_artifacts_dir and debug_selected:
                try:
                    dbg_dir = write_component_debug_artifacts(
                        subpath=sp,
                        shape=sh,
                        angle_deg=ang,
                        axis=axis,
                        origin=origin,
                        plane_face=plane_face,
                        tol=tol,
                        per_solid_clean=per_solid_clean,
                        extent=extent,
                        debug_artifacts_dir=debug_artifacts_dir,
                    )
                    print(f"  [debug] {_sanitize_layer(sp)} -> {dbg_dir}")
                except Exception as ex:
                    print(f"  [debug] FAILED {sp}: {type(ex).__name__}: {ex}")

            edge_source = "trimmed" if keep_slab is not None else "full"
            edges = slice_component_edges(
                sh,
                plane_face=plane_face,
                keep_slab=keep_slab,
                tol=tol,
                per_solid_clean=per_solid_clean,
            )
            trimmed_candidate_edges = list(edges) if keep_slab is not None else None
            full_candidate_edges = None
            if not edges:
                audit_rows.append(
                    {
                        "angle_deg": f"{ang:.8f}",
                        "subpath": sp,
                        "label_path": label_path,
                        "status": "NO_EDGES",
                        "edge_source": edge_source,
                        "edge_count": "0",
                        "mode": "",
                        "dxf_path": "",
                    }
                )
                continue

            total_with_edges += 1
            safe = _sanitize_layer(sp, max_len=180)

            if slice_fcstd_dir and debug_selected:
                try:
                    fcstd_path = os.path.join(
                        os.path.abspath(slice_fcstd_dir),
                        f"slice_{ang:06.2f}deg",
                        safe,
                        "final_slice.FCStd",
                    )
                    _save_export_slice_fcstd(
                        fcstd_path,
                        shape=sh,
                        final_edges=edges,
                        trimmed_edges=trimmed_candidate_edges,
                        full_edges=full_candidate_edges,
                    )
                    print(f"  [slice-fcstd] {safe} -> {fcstd_path}")
                except Exception as ex:
                    print(f"  [slice-fcstd] FAILED {safe}: {type(ex).__name__}: {ex}")

            # decide output destination based on wire grouping success
            # wire-grouped -> normal dir; raw only if grouping fails -> raw dir
            dxf_path_normal = os.path.join(angle_dir, f"{safe}.dxf")
            dxf_path_raw = os.path.join(angle_dir_raw, f"{safe}.dxf")

            # Try wire-grouped first; if it fails, write raw into _raw dir
            mode, prim, gsum = export_component_wiregrouped_or_raw(
                subpath=sp,
                edges=edges,
                dxf_path=dxf_path_normal,  # will be used if WIRE_GROUPED
                angle_deg=ang,
                axis=axis,
                origin=origin,
                tol=tol,
            )

            if mode == "WIRE_GROUPED":
                audit_rows.append(
                    {
                        "angle_deg": f"{ang:.8f}",
                        "subpath": sp,
                        "label_path": label_path,
                        "status": "EXPORTED",
                        "edge_source": edge_source,
                        "edge_count": str(len(edges)),
                        "mode": mode,
                        "dxf_path": dxf_path_normal,
                    }
                )
                print(
                    f"  [comp] {safe} mode=WIRE source={edge_source} "
                    f"groups={gsum} edges={len(edges)} prim={prim}"
                )
            else:
                # re-export to raw location (so raw is only written on failures)
                failed_wiregroup += 1
                mode2, prim2, gsum2 = export_component_wiregrouped_or_raw(
                    subpath=sp,
                    edges=edges,
                    dxf_path=dxf_path_raw,
                    angle_deg=ang,
                    axis=axis,
                    origin=origin,
                    tol=tol,
                )
                audit_rows.append(
                    {
                        "angle_deg": f"{ang:.8f}",
                        "subpath": sp,
                        "label_path": label_path,
                        "status": "EXPORTED",
                        "edge_source": edge_source,
                        "edge_count": str(len(edges)),
                        "mode": mode2,
                        "dxf_path": dxf_path_raw,
                    }
                )
                # mode2 will still be RAW; keep the RAW print
                print(
                    f"  [comp] {safe} mode=RAW source={edge_source} "
                    f"groups={gsum2} edges={len(edges)} prim={prim2}"
                )

    frac = (failed_wiregroup / total_with_edges) if total_with_edges else 0.0
    print(
        f"[summary] components_with_edges={total_with_edges} "
        f"wiregroup_failed={failed_wiregroup} "
        f"fraction_failed={frac:.6f}"
    )

    if component_audit_csv:
        apath = os.path.abspath(component_audit_csv)
        os.makedirs(os.path.dirname(apath), exist_ok=True)
        fields = [
            "angle_deg",
            "subpath",
            "label_path",
            "status",
            "edge_source",
            "edge_count",
            "mode",
            "dxf_path",
        ]
        with open(apath, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in audit_rows:
                w.writerow(r)
        print(f"[audit] wrote {len(audit_rows)} rows -> {apath}")

def _export_one_angle_worker(args):
    """
    Runs in a separate process.
    """
    (
        step_path,
        out_dir,
        ang,
        axis,
        origin,
        tol,
        trim_u_ge_0,
        per_solid_clean,
        component_match,
        component_regex,
        require_single_component_match,
        debug_artifacts_dir,
        debug_component_tokens,
        component_audit_csv,
        slice_fcstd_dir,
    ) = args

    # Load in this process
    doc = load_step(step_path, doc_name=f"nb_step_{ang:.2f}")
    root = find_step_root(doc)

    # Export only this single angle
    export_slices_per_component(
        root=root,
        out_dir=out_dir,
        angles_deg=[ang],
        axis=axis,
        origin=origin,
        tol=tol,
        trim_u_ge_0=trim_u_ge_0,
        per_solid_clean=per_solid_clean,
        component_match=component_match,
        component_regex=component_regex,
        require_single_component_match=require_single_component_match,
        debug_artifacts_dir=debug_artifacts_dir,
        debug_component_tokens=debug_component_tokens,
        component_audit_csv=component_audit_csv,
        slice_fcstd_dir=slice_fcstd_dir,
    )
    return ang

def export_slices_per_component_parallel(
    *,
    step_path: str,
    out_dir: str,
    angles_deg,
    axis: str = "Y",
    origin=(0.0, 0.0, 0.0),
    tol: float = DEFAULT_TOL,
    trim_u_ge_0: bool = True,
    per_solid_clean: bool = True,
    workers: int | None = None,
    component_match: str | None = None,
    component_regex: str | None = None,
    require_single_component_match: bool = False,
    debug_artifacts_dir: str | None = None,
    debug_component_tokens=None,
    component_audit_csv: str | None = None,
    slice_fcstd_dir: str | None = None,
):
    angles = [float(a) for a in angles_deg]
    def _audit_path_for_angle(ang: float):
        if not component_audit_csv:
            return None
        if len(angles) == 1:
            return component_audit_csv
        root, ext = os.path.splitext(component_audit_csv)
        if not ext:
            ext = ".csv"
        return f"{root}_angle_{ang:06.2f}{ext}"

    args = [
        (
            step_path,
            out_dir,
            ang,
            axis,
            origin,
            tol,
            trim_u_ge_0,
            per_solid_clean,
            component_match,
            component_regex,
            require_single_component_match,
            debug_artifacts_dir,
            debug_component_tokens,
            _audit_path_for_angle(ang),
            slice_fcstd_dir,
        )
        for ang in angles
    ]

    n = int(workers) if workers is not None else max(1, (os.cpu_count() or 2) - 1)
    if n < 1:
        raise ValueError(f"workers must be >= 1 (got {n})")

    # Debug path: force true serial execution in current process.
    if n == 1:
        print("[run] serial mode (workers=1)")
        for a in args:
            ang = _export_one_angle_worker(a)
            print(f"[done] angle={ang:.2f}")
        return

    # Use spawn for safety with FreeCAD
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n) as pool:
        for ang in pool.imap_unordered(_export_one_angle_worker, args):
            print(f"[done] angle={ang:.2f}")


if __name__ == "__main__":
    args = _parse_cli_args()
    default_angles = _default_angles_with_offset(args.angle_offset_deg)
    if args.single_angle_deg is not None:
        angles = [float(args.single_angle_deg)]
    elif args.first_angle_only:
        angles = [float(default_angles[0])]
    else:
        angles = default_angles

    print(
        f"[run] angles={', '.join(f'{a:.5f}' for a in angles)} "
        f"workers={args.workers} tol={args.tol} trim_u_ge_0={not args.no_trim_u_ge_0} "
        f"component={args.component!r} component_regex={args.component_regex!r} "
        f"debug_artifacts_dir={args.debug_artifacts_dir!r} debug_component={args.debug_component!r} "
        f"component_audit_csv={args.component_audit_csv!r} "
        f"slice_fcstd_dir={args.slice_fcstd_dir!r}"
    )

    export_slices_per_component_parallel(
        step_path=args.step_path,
        out_dir=args.out_dir,
        angles_deg=angles,
        axis="Y",
        origin=(0.0, 0.0, 0.0),
        tol=args.tol,
        trim_u_ge_0=(not args.no_trim_u_ge_0),
        per_solid_clean=True,
        workers=args.workers,
        component_match=args.component,
        component_regex=args.component_regex,
        require_single_component_match=args.require_single_component_match,
        debug_artifacts_dir=args.debug_artifacts_dir,
        debug_component_tokens=args.debug_component,
        component_audit_csv=args.component_audit_csv,
        slice_fcstd_dir=args.slice_fcstd_dir,
    )
