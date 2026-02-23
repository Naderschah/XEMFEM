import os
import re
import math

import FreeCAD
import Import
import Part

import ezdxf


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

def collect_leaf_geometry_subpaths(root):
    out = []
    seen = set()

    def has_nonnull_shape(o) -> bool:
        try:
            sh = getattr(o, "Shape", None)
            return sh is not None and hasattr(sh, "isNull") and not sh.isNull()
        except Exception:
            return False

    def rec(parent, toks):
        try:
            children = list(getattr(parent, "Group", []) or [])
        except Exception:
            children = []

        if not children:
            if parent is not root and has_nonnull_shape(parent):
                name = str(getattr(parent, "Name", "") or "")
                if name:
                    out.append(".".join(toks + [name]) + ".")
            return

        for ch in children:
            oid = id(ch)
            if oid in seen:
                continue
            seen.add(oid)

            name = str(getattr(ch, "Name", "") or "")
            if not name:
                continue

            rec(ch, toks + [name])

    rec(root, [])
    return out


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


# -----------------------
# Slice per component (NO merging that can drop duplicates)
# -----------------------

def slice_component_edges(
    shape,
    *,
    plane_face,
    keep_slab=None,
    tol: float = 1e-7,
    per_solid_clean: bool = True,
):
    """
    Returns a flat list of edges for THIS component only.
    Key behavior: preserve duplicates by NOT doing compound removeSplitter/cleanTolerance.
    """
    tol_eff = float(max(tol, 1e-6))

    try:
        solids = list(shape.Solids)
    except Exception:
        solids = []
    if not solids:
        return []

    out_edges = []
    for sol in solids:
        sol_trim = sol
        if keep_slab is not None:
            try:
                sol_trim = sol.common(keep_slab)
            except Exception:
                continue
            if sol_trim is None or sol_trim.isNull() or len(sol_trim.Solids) == 0:
                continue

        try:
            sec = sol_trim.section(plane_face)
        except Exception:
            continue
        if sec is None or sec.isNull():
            continue

        # minimal cleaning, per-solid only (does not merge across solids)
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

def _add_spline_from_edge(msp, edge, layer, origin3, A, R):
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

def _add_ellipse_from_edge(msp, edge, layer, origin3, A, R):
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

    msp.add_ellipse(
        center=center2,
        major_axis=major_vec2,
        ratio=ratio,
        start_param=start_param,
        end_param=end_param,
        dxfattribs={"layer": layer},
    )
    return True


def export_component_edges_to_dxf(
    *,
    subpath: str,
    edges,
    dxf_path: str,
    angle_deg: float,
    axis: str,
    origin,
    tol: float,
):
    origin3 = FreeCAD.Vector(*origin)
    A, R = _plane_frame(angle_deg=angle_deg, axis=axis)

    dxf = ezdxf.new("R2010")
    msp = dxf.modelspace()

    layer = _sanitize_layer(subpath)
    if layer not in dxf.layers:
        dxf.layers.new(layer)

    prim_counts = {"LINE": 0, "ARC": 0, "CIRCLE": 0, "ELLIPSE": 0, "SPLINE": 0, "FAIL": 0}

    for e in edges:
        c = e.Curve
        cname = c.__class__.__name__
        try:
            if cname in ("Line", "GeomLine"):
                p0 = e.Vertexes[0].Point
                p1 = e.Vertexes[-1].Point
                msp.add_line(_to_2d(p0, origin3, A, R), _to_2d(p1, origin3, A, R), dxfattribs={"layer": layer})
                prim_counts["LINE"] += 1
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
                p0_2 = _to_2d(p0, origin3, A, R)
                p1_2 = _to_2d(p1, origin3, A, R)

                a1 = _angle_deg_in_plane(p0_2, center2)
                a2 = _angle_deg_in_plane(p1_2, center2)

                ccw = _arc_is_ccw_in_export_coords(e, c, center2, origin3, A, R)
                if not ccw:
                    a1, a2 = a2, a1

                msp.add_arc(center2, r, a1, a2, dxfattribs={"layer": layer})
                prim_counts["ARC"] += 1
                continue

            if cname in ("Ellipse", "GeomEllipse") or (
                hasattr(c, "MajorRadius") and hasattr(c, "MinorRadius") and hasattr(c, "MajorAxis")
            ):
                if _add_ellipse_from_edge(msp, e, layer, origin3, A, R):
                    prim_counts["ELLIPSE"] += 1
                    continue

            _add_spline_from_edge(msp, e, layer, origin3, A, R)
            prim_counts["SPLINE"] += 1

        except Exception:
            prim_counts["FAIL"] += 1

    dxf.saveas(dxf_path)
    return prim_counts


# -----------------------
# Export: angle subdir + one DXF per component
# -----------------------

def export_slices_per_component(
    *,
    root,
    out_dir: str,
    angles_deg,
    axis: str = "Y",
    origin=(0.0, 0.0, 0.0),
    tol: float = 1e-7,
    trim_u_ge_0: bool = True,
    per_solid_clean: bool = True,
):
    os.makedirs(out_dir, exist_ok=True)
    subpaths = collect_leaf_geometry_subpaths(root)

    diag = _compute_global_bbox_diag(root)
    extent = diag * 3.0
    thickness = diag * 0.25

    for ang in angles_deg:
        ang = float(ang)
        angle_dir = os.path.join(out_dir, f"slice_{ang:06.2f}deg")
        os.makedirs(angle_dir, exist_ok=True)

        N = _slice_plane_normal(ang, axis=axis)
        ox, oy, oz = origin
        plane_face = Part.Plane(FreeCAD.Vector(ox, oy, oz), N).toShape()

        keep_slab = None
        if trim_u_ge_0:
            keep_slab = make_keep_slab_u_ge_0(
                angle_deg=ang, axis=axis, origin=origin, extent=extent, thickness=thickness
            )

        print(f"[angle] {ang:.2f}° → {angle_dir}")

        for sp in subpaths:
            try:
                sh = root.getSubObject(sp)
            except Exception:
                continue
            if sh is None or not hasattr(sh, "isNull") or sh.isNull():
                continue

            edges = slice_component_edges(
                sh,
                plane_face=plane_face,
                keep_slab=keep_slab,
                tol=tol,
                per_solid_clean=per_solid_clean,
            )
            if not edges:
                continue

            # file name per component
            safe = _sanitize_layer(sp, max_len=180)
            dxf_path = os.path.join(angle_dir, f"{safe}.dxf")

            prim = export_component_edges_to_dxf(
                subpath=sp,
                edges=edges,
                dxf_path=dxf_path,
                angle_deg=ang,
                axis=axis,
                origin=origin,
                tol=tol,
            )

            # Debug: confirms "every edge we got was attempted"
            print(f"  [comp] {safe} edges={len(edges)} prim={prim}")


# -----------------------
# Example run
# -----------------------

STEP_PATH = "/work/geometry/createGeometryFromCAD/CAD_files/XENT-TPC_20250428.STEP"
doc = load_step(STEP_PATH, doc_name="nb_step")
root = find_step_root(doc)

export_slices_per_component(
    root=root,
    out_dir="/work/geometry/createGeometryFromCAD/DXF_slices_parts",
    angles_deg=[15 * (i + 1) for i in range(0, 24)],
    axis="Y",
    origin=(0.0, 0.0, 0.0),
    tol=1e-7,
    trim_u_ge_0=True,
    per_solid_clean=True,
)