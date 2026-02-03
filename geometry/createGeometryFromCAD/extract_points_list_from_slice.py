#!/usr/bin/env python3
# needs pip install matplotlib

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union, Set

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import FreeCAD  # type: ignore
import Import   # type: ignore
import Part     # type: ignore

Pt2 = Tuple[float, float]
Instr = Union[List[float], List[object]]  # [x,y] or ["line",x,y] or ["arc",x,y,cx,cy,cw] etc.


# -------------------- Data model --------------------

@dataclass(frozen=True)
class ComponentSpec:
    parent: str
    subcomponent: str
    material: str
    mesh: bool
    bc: str
    comment: str


# -------------------- Small helpers --------------------

def _safe(obj: Any, attr: str, default: Any = None) -> Any:
    try:
        return getattr(obj, attr)
    except Exception:
        return default

def _truthy_mesh(value: str) -> bool:
    v = (value or "").strip().lower()
    if v in {"1", "true", "yes", "y", "on"}:
        return True
    if v in {"0", "false", "no", "n", "off", ""}:
        return False
    try:
        return int(v) != 0
    except Exception:
        return False

def _sanitize_filename(s: str) -> str:
    s = (s or "").replace("*", "")
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("._-")
    return s or "unnamed"

def _angles_deg_default() -> List[float]:
    return [30.0 + 15.0 * k for k in range(12)]

def _canon(obj: Any, nd: int = 10) -> Any:
    if isinstance(obj, float):
        return round(obj, nd)
    if isinstance(obj, (int, str, bool)) or obj is None:
        return obj
    if isinstance(obj, list):
        return [_canon(x, nd=nd) for x in obj]
    if isinstance(obj, tuple):
        return [_canon(x, nd=nd) for x in obj]
    if isinstance(obj, dict):
        return {k: _canon(obj[k], nd=nd) for k in sorted(obj.keys())}
    return str(obj)

def _signature_instances(per_obj_outputs: List[Dict[str, Any]], nd: int = 10) -> str:
    per_obj_outputs_sorted = sorted(
        per_obj_outputs,
        key=lambda d: (d.get("object_label", ""), d.get("object_name", ""), d.get("subpath", ""))
    )
    payload = {"instances": _canon(per_obj_outputs_sorted, nd=nd)}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))

def _build_spec_index(targets: List[ComponentSpec]) -> Dict[str, ComponentSpec]:
    return {t.subcomponent: t for t in targets}


# -------------------- CSV loading --------------------

def load_component_csv(csv_path: str = "./naming_convention.csv") -> Tuple[List[ComponentSpec], Dict[str, List[ComponentSpec]]]:
    specs: List[ComponentSpec] = []
    by_bc_material: Dict[str, List[ComponentSpec]] = {}

    with open(csv_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        required = {"Parent", "Subcomponent", "Material", "Mesh", "BC", "Comment"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}. Found: {reader.fieldnames}")

        for row in reader:
            parent = (row.get("Parent") or "").strip()
            sub = (row.get("Subcomponent") or "").strip()
            material = (row.get("Material") or "").strip() or "NA"
            mesh = _truthy_mesh(row.get("Mesh") or "")
            bc = (row.get("BC") or "").strip() or "NA"
            comment = (row.get("Comment") or "").strip()

            if not parent and not sub:
                continue

            spec = ComponentSpec(
                parent=parent,
                subcomponent=sub,
                material=material,
                mesh=mesh,
                bc=bc,
                comment=comment,
            )
            specs.append(spec)

            key = f"{bc}|{material}"
            by_bc_material.setdefault(key, []).append(spec)

    return specs, by_bc_material

def mesh_targets_from_specs(specs: List[ComponentSpec]) -> List[ComponentSpec]:
    return [s for s in specs if s.mesh]

def load_component_specs_next_to_step(step_path: str, csv_name: Optional[str] = None) -> Tuple[str, List[ComponentSpec], Dict[str, List[ComponentSpec]]]:
    # You currently hardcode naming_convention.csv; keep that behavior.
    csv_path = "./naming_convention.csv"
    specs, by_group = load_component_csv(csv_path)
    return csv_path, specs, by_group


# -------------------- STEP loading --------------------

def load_step(step_path: str, doc_name: str) -> Any:
    doc = FreeCAD.newDocument(doc_name)
    Import.insert(step_path, doc.Name)
    doc.recompute()
    return doc

def find_step_root(doc: Any) -> Any:
    # Choose a root that is not in anyone’s InList (same as your approach)
    root = None
    for o in getattr(doc, "Objects", []) or []:
        if not getattr(o, "InList", []):
            root = o
            break
    if root is None:
        raise RuntimeError("Could not determine STEP root object")
    return root


# -------------------- Indexing: labels & paths --------------------

def build_label_object_index(doc: Any) -> Dict[str, List[Any]]:
    """
    Label -> [objects]
    Includes App::Part containers and leaf objects with non-null Shape.
    """
    idx: Dict[str, List[Any]] = {}
    seen: set[int] = set()

    def add(obj: Any):
        if obj is None:
            return
        oid = id(obj)
        if oid in seen:
            return
        seen.add(oid)

        lab = str(getattr(obj, "Label", "") or "").strip()
        if lab:
            idx.setdefault(lab, []).append(obj)

        try:
            grp = list(getattr(obj, "Group", []))
        except Exception:
            grp = []
        for ch in grp:
            add(ch)

    for o in getattr(doc, "Objects", []) or []:
        try:
            type_id = str(getattr(o, "TypeId", "") or "")
        except Exception:
            type_id = ""

        is_part = (type_id == "App::Part")
        has_shape = False
        if not is_part:
            try:
                sh = getattr(o, "Shape", None)
                has_shape = (sh is not None and not sh.isNull())
            except Exception:
                has_shape = False

        if is_part or has_shape:
            add(o)
        else:
            try:
                grp = list(getattr(o, "Group", []))
            except Exception:
                grp = []
            for ch in grp:
                add(ch)

    return idx

def build_object_subpath_index_from_root(root: Any) -> Dict[Any, str]:
    """
    Build object -> getSubObject()-compatible path.

    We walk root.Group recursively and build a token chain using .Name.
    Returned path format is: "TokenA.TokenB.TokenC." (trailing dot required for containers/subobjects).

    IMPORTANT:
      - This assumes the STEP import hierarchy is represented via Group containers.
      - If an object appears multiple times, the first discovered path is used.
    """
    idx: Dict[Any, str] = {}
    seen: set[int] = set()

    def children(o: Any) -> List[Any]:
        try:
            return list(getattr(o, "Group", [])) or []
        except Exception:
            return []

    def rec(o: Any, chain_tokens: List[str]):
        oid = id(o)
        if oid in seen:
            return
        seen.add(oid)

        name = str(getattr(o, "Name", "") or "").strip()
        if not name:
            # fall back to Label (less stable)
            name = str(getattr(o, "Label", "") or "").strip() or "<?>"

        new_chain = chain_tokens + [name]
        path = ".".join(new_chain) + "."
        idx[o] = path

        for ch in children(o):
            rec(ch, new_chain)

    rec(root, [])
    return idx

def summarize_targets(
    targets: List[ComponentSpec],
    resolved: Dict[str, List[str]],
    unresolved: List[str],
) -> None:
    print(f"Mesh targets (CSV): {len(targets)}")
    print(f"Resolved targets:   {sum(1 for k in resolved if resolved[k])} labels")
    print(f"Unresolved targets: {len(unresolved)} labels")
    if unresolved:
        print("Unresolved labels:")
        for s in unresolved:
            print(f"  - {s}")


# -------------------- Resolution: CSV -> paths --------------------

def resolve_targets_from_csv_paths(
    doc: Any,
    root: Any,
    step_path: str,
    csv_name_or_path: str | None,
    label_index: Dict[str, List[Any]],
    obj2path: Dict[Any, str],
) -> Tuple[List[ComponentSpec], Dict[str, List[str]], List[str], List[ComponentSpec], Any]:
    """
    Resolve CSV mesh targets to *SubObject path strings* (root.getSubObject(path)).

    Returns:
      targets, resolved, unresolved, specs, by_group

    Where:
      resolved[sub_label] = [ "RootToken.ChildToken.LeafToken.", ... ]
    """
    csv_path, specs, by_group = load_component_specs_next_to_step(step_path, csv_name_or_path)
    targets = mesh_targets_from_specs(specs)

    resolved: Dict[str, List[str]] = {}
    unresolved: List[str] = []

    def objs_for_label(sub: str) -> List[Any]:
        # exact label
        if sub in label_index:
            return list(label_index[sub])
        # wildcard suffix "***" => prefix match
        if sub.endswith("***"):
            prefix = sub[:-3]
            out: List[Any] = []
            for lab, lst in label_index.items():
                if lab.startswith(prefix):
                    out.extend(lst)
            return out
        return []

    for t in targets:
        sub = (t.subcomponent or "").strip()
        if not sub:
            resolved[sub] = []
            unresolved.append(sub)
            continue

        objs = objs_for_label(sub)
        paths: List[str] = []
        for o in objs:
            p = obj2path.get(o)
            if p:
                paths.append(p)

        # de-dup while preserving order
        seenp: set[str] = set()
        paths_unique: List[str] = []
        for p in paths:
            if p not in seenp:
                seenp.add(p)
                paths_unique.append(p)

        resolved[sub] = paths_unique
        if not paths_unique:
            unresolved.append(sub)

    print(f"STEP: {os.path.abspath(step_path)}")
    print(f"CSV:  {os.path.abspath(csv_path)}")
    summarize_targets(targets, resolved, unresolved)

    return targets, resolved, unresolved, specs, by_group


# -------------------- Shape retrieval via root.getSubObject --------------------

def _is_toposhape(x: Any) -> bool:
    return (
        x is not None
        and hasattr(x, "isNull") and callable(getattr(x, "isNull", None))
        and hasattr(x, "section") and callable(getattr(x, "section", None))
        and hasattr(x, "BoundBox")
    )

def resolve_shape_via_root_subpath(root: Any, subpath: str) -> Part.TopoShape:
    """
    Resolve a subpath from a root object, returning a Part.TopoShape in the
    *instance/global* frame by applying the accumulated placement from the path.

    Key idea:
      - root.getSubObject(path) gives you the referenced thing
      - root.getSubObject(path, retType=3) (when supported) gives the accumulated placement
        along the path (critical for Links / nested containers)

    If retType=3 is not supported in this FreeCAD build, we fall back to returning the
    shape as-is (but then Links may be wrong).
    """
    s = str(subpath or "").strip()
    if not s.endswith("."):
        s += "."
    pyObject = -1
    mat = FreeCAD.Matrix()
    sub = root.getSubObject(s,  matrix=mat, transform = True)
    print(type(sub))
    if sub is None:
        raise RuntimeError(f"getSubObject returned None for subpath={s!r}")
    print(type(sub))
    # --- obtain a TopoShape from `sub` ---
    sh = None

    # Case 1: already a TopoShape
    if _is_toposhape(sub) and not sub.isNull():
        sh = sub

    # Case 2: document object with Shape
    elif hasattr(sub, "Shape"):
        try:
            cand = sub.Shape
        except Exception:
            cand = None
        if _is_toposhape(cand) and not cand.isNull():
            sh = cand

    # Case 3: subshape with toShape()
    if sh is None:
        to_shape = getattr(sub, "toShape", None)
        if callable(to_shape):
            try:
                cand = sub.toShape()
            except Exception:
                cand = None
            if _is_toposhape(cand) and not cand.isNull():
                sh = cand

    if sh is None:
        raise RuntimeError(
            f"resolve_shape_via_root_subpath: could not obtain TopoShape for subpath={s!r} "
            f"(type={type(sub).__name__})"
        )

    # --- get accumulated placement and apply it ---
    pl = None
    try:
        # This is the feature you said is "always provided" for the working case.
        pl = root.getSubObject(s, retType=3)
    except Exception:
        pl = None

    # If pl is a Placement-like object, apply its matrix to the shape.
    if pl is not None and hasattr(pl, "toMatrix"):
        try:
            m = pl.toMatrix()
            # transformGeometry returns a new transformed shape
            sh2 = sh.copy()
            try:
                sh2 = sh2.transformGeometry(m)
            except Exception:
                # some builds use transformShape instead
                sh2.transformShape(m)
            return sh2
        except Exception:
            # If something goes wrong, fall back to the untransformed shape
            return sh

    return sh

# -------------------- Slicing --------------------

def slice_paths_with_vertical_plane(
    root: Any,
    subpaths: List[str],
    angle_deg: float,
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    axis: str = "Y",
    tol: float = 1e-7,
) -> List[Part.Edge]:
    ox, oy, oz = origin
    a = math.radians(angle_deg)
    ax = axis.upper().strip()

    if ax == "Y":
        nx, ny, nz = math.sin(a), 0.0, -math.cos(a)
    elif ax == "Z":
        nx, ny, nz = -math.sin(a), math.cos(a), 0.0
    elif ax == "X":
        nx, ny, nz = 0.0, -math.sin(a), math.cos(a)
    else:
        raise ValueError(f"axis must be 'X','Y','Z' (got {axis!r})")

    plane = Part.Plane(FreeCAD.Vector(ox, oy, oz), FreeCAD.Vector(nx, ny, nz))
    plane_face = plane.toShape()

    section_edges: List[Part.Edge] = []

    for sp in subpaths:
        try:
            shape = resolve_shape_via_root_subpath(root, sp)
        except Exception as e:
            print(e)
            continue

        if shape is None or shape.isNull():
            continue

        try:
            sec = shape.section(plane_face)
        except Exception:
            continue

        if sec is None or sec.isNull():
            continue

        try:
            sec = sec.removeSplitter()
        except Exception:
            pass
        try:
            sec = sec.cleanTolerance(tol)
        except Exception:
            pass

        try:
            section_edges.extend(list(sec.Edges))
        except Exception:
            pass

    return section_edges


# -------------------- Projection to 2D (axis-aware) --------------------

def _v3_to_2d_in_vertical_plane(
    v: FreeCAD.Vector,
    angle_deg: float,
    axis: str = "Y",
) -> Pt2:
    a = math.radians(angle_deg)
    ax = axis.upper().strip()

    if ax == "X":
        A  = (1.0, 0.0, 0.0)
        e1 = (0.0, 1.0, 0.0)
        e2 = (0.0, 0.0, 1.0)
    elif ax == "Y":
        A  = (0.0, 1.0, 0.0)
        e1 = (1.0, 0.0, 0.0)
        e2 = (0.0, 0.0, 1.0)
    elif ax == "Z":
        A  = (0.0, 0.0, 1.0)
        e1 = (1.0, 0.0, 0.0)
        e2 = (0.0, 1.0, 0.0)
    else:
        raise ValueError(f"axis must be 'X','Y','Z' (got {axis!r})")

    rx = math.cos(a) * e1[0] + math.sin(a) * e2[0]
    ry = math.cos(a) * e1[1] + math.sin(a) * e2[1]
    rz = math.cos(a) * e1[2] + math.sin(a) * e2[2]

    vx, vy, vz = float(v.x), float(v.y), float(v.z)
    axial = vx * A[0] + vy * A[1] + vz * A[2]
    u = vx * rx + vy * ry + vz * rz

    return (float(u), float(axial))

def _close_enough(p: Pt2, q: Pt2, tol: float) -> bool:
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 <= tol * tol

def _pt_key(p: Pt2, tol: float) -> Tuple[int, int]:
    return (int(round(p[0] / tol)), int(round(p[1] / tol)))

def _edge_endpoints_2d_oriented(
    e: Part.Edge,
    fwd: bool,
    angle_deg: float,
    axis: str = "Y",
) -> Tuple[Pt2, Pt2]:
    u0, u1 = e.ParameterRange
    pA3 = e.valueAt(u0 if fwd else u1)
    pB3 = e.valueAt(u1 if fwd else u0)
    return (
        _v3_to_2d_in_vertical_plane(pA3, angle_deg, axis=axis),
        _v3_to_2d_in_vertical_plane(pB3, angle_deg, axis=axis),
    )

def _circle_center_2d(e: Part.Edge, angle_deg: float, axis: str = "Y") -> Pt2:
    c = e.Curve
    if not isinstance(c, Part.Circle):
        raise TypeError(f"Edge curve is not a Circle: {type(c).__name__}")
    return _v3_to_2d_in_vertical_plane(c.Center, angle_deg, axis=axis)

def _arc_cw_flag(p0: Pt2, p1: Pt2, c: Pt2) -> bool:
    ax, ay = p0[0] - c[0], p0[1] - c[1]
    bx, by = p1[0] - c[0], p1[1] - c[1]
    cross = ax * by - ay * bx
    return cross < 0.0


# -------------------- Curve instruction builders (axis-aware) --------------------

def make_bspline_instruction(
    c: Part.BSplineCurve,
    pA: Pt2,
    angle_deg: float,
    axis: str = "Y",
) -> Instr:
    degree = int(getattr(c, "Degree"))

    poles3 = list(c.getPoles())
    poles2 = [list(_v3_to_2d_in_vertical_plane(v, angle_deg, axis=axis)) for v in poles3]

    knots = [float(k) for k in c.getKnots()]
    mults = [int(m) for m in c.getMultiplicities()]

    rational = bool(getattr(c, "isRational")())
    periodic = bool(getattr(c, "isPeriodic")())

    payload = {
        "degree": degree,
        "poles": poles2,
        "knots": knots,
        "mults": mults,
        "periodic": periodic,
        "rational": rational,
    }

    if rational:
        w = getattr(c, "getWeights", None)
        if callable(w):
            payload["weights"] = [float(x) for x in c.getWeights()]
        else:
            raise RuntimeError("BSplineCurve reports rational=True but has no getWeights()")

    return ["bspline", pA[0], pA[1], payload]

def _ellipse_params_2d(c: Part.Ellipse, angle_deg: float, axis: str = "Y"):
    center2d = _v3_to_2d_in_vertical_plane(c.Center, angle_deg, axis=axis)
    a = float(c.MajorRadius)
    b = float(c.MinorRadius)

    phi_deg = 0.0
    try:
        p_major = c.value(0.0)
        p_major2d = _v3_to_2d_in_vertical_plane(p_major, angle_deg, axis=axis)
        du = float(p_major2d[0] - center2d[0])
        dv = float(p_major2d[1] - center2d[1])
        phi_deg = float(math.degrees(math.atan2(dv, du)))
    except Exception:
        try:
            p_major = c.valueAt(0.0)
            p_major2d = _v3_to_2d_in_vertical_plane(p_major, angle_deg, axis=axis)
            du = float(p_major2d[0] - center2d[0])
            dv = float(p_major2d[1] - center2d[1])
            phi_deg = float(math.degrees(math.atan2(dv, du)))
        except Exception:
            phi_deg = 0.0

    return center2d, a, b, phi_deg

def make_ellipse_instruction(
    c: Part.Ellipse,
    pA: Pt2,
    pB: Pt2,
    angle_deg: float,
    axis: str = "Y",
) -> Instr:
    center, a, b, phi_deg = _ellipse_params_2d(c, angle_deg, axis=axis)
    cw = _arc_cw_flag(pA, pB, center)
    return ["ellipse", pA[0], pA[1], center[0], center[1], a, b, phi_deg, cw]

def _hyperbola_params_2d(c: Part.Hyperbola, angle_deg: float, axis: str = "Y"):
    center2d = _v3_to_2d_in_vertical_plane(c.Center, angle_deg, axis=axis)
    a = float(c.MajorRadius)
    b = float(c.MinorRadius)

    phi_deg = 0.0
    try:
        t = c.tangent(0.0)  # FreeCAD.Vector direction
        o2d = _v3_to_2d_in_vertical_plane(FreeCAD.Vector(0, 0, 0), angle_deg, axis=axis)
        t2d = _v3_to_2d_in_vertical_plane(t, angle_deg, axis=axis)
        du = float(t2d[0] - o2d[0])
        dv = float(t2d[1] - o2d[1])
        phi_deg = float(math.degrees(math.atan2(dv, du)))
    except Exception:
        try:
            p0 = c.value(0.0)
            p02d = _v3_to_2d_in_vertical_plane(p0, angle_deg, axis=axis)
            du = float(p02d[0] - center2d[0])
            dv = float(p02d[1] - center2d[1])
            phi_deg = float(math.degrees(math.atan2(dv, du)))
        except Exception:
            phi_deg = 0.0

    return center2d, a, b, phi_deg

def make_hyperbola_instruction(
    c: Part.Hyperbola,
    pA: Pt2,
    pB: Pt2,
    angle_deg: float,
    axis: str = "Y",
) -> Instr:
    center, a, b, phi_deg = _hyperbola_params_2d(c, angle_deg, axis=axis)
    cw = _arc_cw_flag(pA, pB, center)
    return ["hyperbola", pA[0], pA[1], center[0], center[1], a, b, phi_deg, cw]


# -------------------- Connectivity + ordering + transcription (axis-aware) --------------------

def sort_edges_into_wires(edges: List[Part.Edge]) -> List[List[Part.Edge]]:
    if not edges:
        return []
    try:
        groups = Part.sortEdges(edges)
    except Exception as ex:
        raise RuntimeError(f"Part.sortEdges failed: {ex!r}")
    return [list(g) for g in groups if g]

def split_and_order_edges_into_chains(
    wire_edges: List[Part.Edge],
    angle_deg: float,
    tol: float,
    axis: str = "Y",
) -> List[List[Tuple[Part.Edge, bool]]]:
    if not wire_edges:
        return []

    endpoints: List[Tuple[Part.Edge, Tuple[int, int], Tuple[int, int]]] = []
    for e in wire_edges:
        u0, u1 = e.ParameterRange
        p0 = _v3_to_2d_in_vertical_plane(e.valueAt(u0), angle_deg, axis=axis)
        p1 = _v3_to_2d_in_vertical_plane(e.valueAt(u1), angle_deg, axis=axis)
        k0 = _pt_key(p0, tol)
        k1 = _pt_key(p1, tol)
        endpoints.append((e, k0, k1))

    v2e: Dict[Tuple[int, int], List[int]] = {}
    for i, (_e, k0, k1) in enumerate(endpoints):
        v2e.setdefault(k0, []).append(i)
        v2e.setdefault(k1, []).append(i)

    e_adj: Dict[int, Set[int]] = {i: set() for i in range(len(endpoints))}
    for eis in v2e.values():
        for i in eis:
            for j in eis:
                if i != j:
                    e_adj[i].add(j)

    comps: List[List[int]] = []
    unvisited = set(range(len(endpoints)))
    while unvisited:
        seed = next(iter(unvisited))
        stack = [seed]
        comp: List[int] = []
        unvisited.remove(seed)
        while stack:
            i = stack.pop()
            comp.append(i)
            for j in e_adj[i]:
                if j in unvisited:
                    unvisited.remove(j)
                    stack.append(j)
        comps.append(comp)

    chains: List[List[Tuple[Part.Edge, bool]]] = []
    for comp in comps:
        adj: Dict[Tuple[int, int], List[Tuple[int, bool]]] = {}
        for idx in comp:
            _e, k0, k1 = endpoints[idx]
            adj.setdefault(k0, []).append((idx, True))
            adj.setdefault(k1, []).append((idx, False))

        branched = [v for v, inc in adj.items() if len(inc) > 2]
        if branched:
            raise RuntimeError(
                f"Branched section graph (deg>2) under tol={tol}; cannot serialize as chains."
            )

        deg1 = [v for v, inc in adj.items() if len(inc) == 1]
        cur_v = deg1[0] if deg1 else next(iter(adj.keys()))

        used: Set[int] = set()
        ordered: List[Tuple[Part.Edge, bool]] = []

        while True:
            candidates = [t for t in adj.get(cur_v, []) if t[0] not in used]
            if not candidates:
                break
            edge_idx, end_is_start = candidates[0]
            used.add(edge_idx)

            e, k0, k1 = endpoints[edge_idx]
            if end_is_start:
                ordered.append((e, True))
                cur_v = k1
            else:
                ordered.append((e, False))
                cur_v = k0

        if len(used) != len(comp):
            raise RuntimeError(f"Failed to order all edges in a component: used {len(used)}/{len(comp)} (tol={tol})")

        chains.append(ordered)

    return chains

def wire_to_pts_instructions_strict(
    ordered_edges: List[Tuple[Part.Edge, bool]],
    angle_deg: float,
    tol: float,
    axis: str = "Y",
) -> List[Instr]:
    if not ordered_edges:
        return []

    pts: List[Instr] = []

    e0, fwd0 = ordered_edges[0]
    u0, u1 = e0.ParameterRange
    p0_3 = e0.valueAt(u0 if fwd0 else u1)
    start0 = _v3_to_2d_in_vertical_plane(p0_3, angle_deg, axis=axis)
    current = start0

    for (e, fwd) in ordered_edges:
        pA, pB = _edge_endpoints_2d_oriented(e, fwd, angle_deg, axis=axis)

        if pts and not _close_enough(pA, current, tol):
            raise RuntimeError("Wire transcription lost continuity (start mismatch).")

        c = e.Curve
        if isinstance(c, Part.Line):
            pts.append(["line", pA[0], pA[1]])

        elif isinstance(c, Part.Circle):
            if _close_enough(pA, pB, tol):
                raise ValueError("Encountered full circle edge; cannot represent as single arc.")
            center = _circle_center_2d(e, angle_deg, axis=axis)
            cw = _arc_cw_flag(pA, pB, center)
            pts.append(["arc", pA[0], pA[1], center[0], center[1], cw])

        elif isinstance(c, Part.Ellipse):
            if _close_enough(pA, pB, tol):
                raise ValueError("Encountered full ellipse edge; cannot represent as single ellipse-arc.")
            pts.append(make_ellipse_instruction(c, pA, pB, angle_deg, axis=axis))

        elif isinstance(c, Part.Hyperbola):
            if _close_enough(pA, pB, tol):
                raise ValueError("Encountered closed hyperbola trim; cannot represent.")
            pts.append(make_hyperbola_instruction(c, pA, pB, angle_deg, axis=axis))

        elif isinstance(c, Part.BSplineCurve):
            pts.append(make_bspline_instruction(c, pA, angle_deg, axis=axis))

        else:
            raise TypeError(
                f"Unsupported curve type: {type(c).__name__} "
                f"(allowed: Line, Circle, Ellipse, Hyperbola, BSplineCurve)."
            )

        current = pB

    return pts

def edges_to_pts_chains_strict(
    edges: List[Part.Edge],
    angle_deg: float,
    tol: float = 1e-5,
    axis: str = "Y",
) -> Tuple[List[List[Instr]], List[List[Instr]], Dict[str, Any] | None]:
    try:
        wire_groups = sort_edges_into_wires(edges)
        all_pts_lists: List[List[Instr]] = []

        for g in wire_groups:
            chains = split_and_order_edges_into_chains(g, angle_deg=angle_deg, tol=tol, axis=axis)
            for ordered in chains:
                pts = wire_to_pts_instructions_strict(ordered, angle_deg=angle_deg, tol=tol, axis=axis)
                if pts:
                    all_pts_lists.append(pts)

        return all_pts_lists, [], None

    except Exception as ex:
        print("Exception in edges to points ", ex)
        return [], [], {
            "type": type(ex).__name__,
            "message": str(ex),
            "traceback": traceback.format_exc(limit=10),
        }


# -------------------- Plotting --------------------

def _instr_xy(instr: Instr) -> Tuple[float, float]:
    if isinstance(instr, list) and instr and isinstance(instr[0], str):
        return float(instr[1]), float(instr[2])
    return float(instr[0]), float(instr[1])

def save_slice_plot_instances(
    per_obj_outputs: List[Dict[str, Any]],
    png_path: str,
    title: str | None = None,
) -> None:
    plt.figure()
    if title:
        plt.title(title)

    any_drawn = False
    for obj_out in per_obj_outputs:
        pts_lists = obj_out.get("pts_lists", []) or []
        if not pts_lists:
            continue
        any_drawn = True

        for pts in pts_lists:
            if not pts:
                continue
            xy = [_instr_xy(i) for i in pts]
            if not xy:
                continue

            xy_closed = xy + [xy[0]]
            xs = [p[0] for p in xy_closed]
            ys = [p[1] for p in xy_closed]

            plt.plot(xs, ys, linewidth=0.8)
            plt.scatter([p[0] for p in xy], [p[1] for p in xy], s=8)

    if not any_drawn:
        plt.text(0.5, 0.5, "no intersection", ha="center", va="center", transform=plt.gca().transAxes)

    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()


# -------------------- Output writing (now uses subpaths) --------------------

def write_pts_outputs_for_all_components(
    root: Any,
    resolved: Dict[str, List[str]],          # sub_label -> [subpath, ...]
    targets: List[ComponentSpec],
    out_base_dir: str,
    angles_deg: List[float] | None = None,
    tol: float = 1e-5,
    axis: str = "Y",
) -> None:
    if angles_deg is None:
        angles_deg = _angles_deg_default()

    pts_root = os.path.join(out_base_dir, "pts_out")
    os.makedirs(pts_root, exist_ok=True)

    spec_index = _build_spec_index(targets)

    for sub_label, subpaths in resolved.items():
        if not subpaths:
            continue

        spec = spec_index.get(sub_label, None)
        bc = (spec.bc if spec else "NA") or "NA"

        bc_dir = os.path.join(pts_root, _sanitize_filename(bc))
        comp_dir = os.path.join(bc_dir, _sanitize_filename(sub_label))
        os.makedirs(comp_dir, exist_ok=True)

        print(f"\nComponent: {sub_label}  (BC={bc})")
        print(f"  Instances(paths): {len(subpaths)}")

        buckets: Dict[str, Dict[str, Any]] = {}

        for ang in angles_deg:
            per_inst_outputs: List[Dict[str, Any]] = []

            # Slice each instance/path independently (keep outputs separate)
            for sp in subpaths:
                edges = slice_paths_with_vertical_plane(
                    root=root,
                    subpaths=[sp],
                    angle_deg=float(ang),
                    origin=(0, 0, 0),
                    axis=axis,
                    tol=tol,
                )

                pts_lists, holes, err = edges_to_pts_chains_strict(
                    edges,
                    angle_deg=float(ang),
                    tol=tol,
                    axis=axis,
                )

                # For inspection, also capture the resolved object if available
                obj_label = ""
                obj_name = ""
                obj_type = ""
                try:
                    o = root.getSubObject(sp if sp.endswith(".") else sp + ".")
                    obj_label = str(_safe(o, "Label", ""))
                    obj_name = str(_safe(o, "Name", ""))
                    obj_type = str(_safe(o, "TypeId", ""))
                except Exception:
                    pass

                per_inst_outputs.append({
                    "subpath": str(sp),
                    "object_name": obj_name,
                    "object_label": obj_label,
                    "object_type": obj_type,
                    "pts_lists": pts_lists,
                    "holes": holes,
                    "error": err,
                })

                del edges

            plot_path = os.path.join(comp_dir, f"slice_{int(round(float(ang))):03d}.png")
            save_slice_plot_instances(per_inst_outputs, plot_path, title=f"{sub_label} @ {float(ang):.1f}°  axis={axis}")

            sig = _signature_instances(per_inst_outputs, nd=10)
            if sig not in buckets:
                buckets[sig] = {"angles": [], "instances": per_inst_outputs}
            buckets[sig]["angles"].append(float(ang))

        for entry in buckets.values():
            angs_sorted = sorted(entry["angles"])
            name = "_".join(f"{int(round(a)):03d}" for a in angs_sorted) + ".json"
            out_path = os.path.join(comp_dir, name)

            payload = {
                "component": sub_label,
                "bc": bc,
                "axis": axis,
                "angles_deg": angs_sorted,
                "instances": entry["instances"],
            }

            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, ensure_ascii=False)

        print(f"  Wrote {len(buckets)} unique slice file(s) -> {comp_dir}")


def default_out_base_dir_from_script(argv0: str) -> str:
    script_dir = os.path.dirname(os.path.abspath(argv0))
    return os.path.join(script_dir, "out")


# -------------------- Main --------------------

def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("step_path", help="Path to STEP file")
    ap.add_argument("--csv", dest="csv_name_or_path", default=None,
                    help="CSV filename (relative to STEP dir) or absolute path")
    ap.add_argument("--json", dest="json_out", default=None,
                    help="Write resolved mapping to JSON (labels->paths)")
    ap.add_argument("--doc-name", default="imported_step",
                    help="FreeCAD document name")
    ap.add_argument("--axis", default="Y", help="Slice axis: X, Y, or Z")
    ap.add_argument("--tol", type=float, default=1e-5, help="Geometric tolerance")
    args = ap.parse_args(argv)

    doc = load_step(args.step_path, args.doc_name)
    root = find_step_root(doc)

    label_index = build_label_object_index(doc)
    obj2path = build_object_subpath_index_from_root(root)

    targets, resolved, unresolved, specs, by_group = resolve_targets_from_csv_paths(
        doc=doc,
        root=root,
        step_path=args.step_path,
        csv_name_or_path=args.csv_name_or_path,
        label_index=label_index,
        obj2path=obj2path,
    )

    if args.json_out:
        out_map: Dict[str, Any] = {}
        for sub_label, paths in resolved.items():
            out_map[sub_label] = [{"subpath": p} for p in (paths or [])]
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(out_map, f, indent=2, ensure_ascii=False)

    out_base = default_out_base_dir_from_script(sys.argv[0])

    write_pts_outputs_for_all_components(
        root=root,
        resolved=resolved,
        targets=targets,
        out_base_dir=out_base,
        angles_deg=None,
        tol=float(args.tol),
        axis=str(args.axis),
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
