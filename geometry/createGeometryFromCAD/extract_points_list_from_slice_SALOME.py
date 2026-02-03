#!/usr/bin/env python3
# SALOME SHAPER TUI script (no GEOM usage)

from __future__ import annotations

import csv
import json
import math
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Set

from salome.shaper import model

from GeomAPI import GeomAPI_Pnt, GeomAPI_Dir, GeomAPI_ShapeExplorer
from GeomAlgoAPI import GeomAlgoAPI_FaceBuilder, GeomAlgoAPI_ShapeTools


# -------------------- Hardcoded IO --------------------

STEP_PATH = "/work/geometry/createGeometryFromCAD/mount/XENT-TPC_20250428.STEP"
CSV_PATH  = "/work/geometry/createGeometryFromCAD/naming_convention.csv"

# Output base directory (mirrors your FreeCAD behavior: "out" next to script)
OUT_BASE_DIR = "/work/geometry/createGeometryFromCAD/out"

# Slice settings (first 6 angles: 30, 45, 60, 75, 90, 105)
DEFAULT_ANGLES_DEG = [30.0 + 15.0 * k for k in range(1)]

# Slicing axis (must support X/Y/Z)
DEFAULT_AXIS = "Y"

# Plane square size used for intersection face (large enough to cut the whole model)
PLANE_FACE_SIZE = 1.0e7


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

def _ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def _write_json(path: str, payload: Any) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)

def _wname(obj: Any) -> str:
    """
    Best-effort: object browser name. In Shaper, many objects expose:
      - obj.data().name()  (Model_Data::name)  :contentReference[oaicite:0]{index=0}
    """
    try:
        d = obj.data()
        n = d.name()
        # name() is a std::wstring; Python usually returns str already
        return str(n)
    except Exception:
        pass
    try:
        return str(obj.name())
    except Exception:
        return ""


# -------------------- CSV loading --------------------

def load_component_csv(csv_path: str) -> List[ComponentSpec]:
    specs: List[ComponentSpec] = []
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

            specs.append(ComponentSpec(
                parent=parent,
                subcomponent=sub,
                material=material,
                mesh=mesh,
                bc=bc,
                comment=comment,
            ))
    return specs

def mesh_targets(specs: List[ComponentSpec]) -> List[ComponentSpec]:
    return [s for s in specs if s.mesh]


# -------------------- STEP import (Shaper) --------------------

def import_step_into_new_part(step_path: str) -> Tuple[Any, Any, Any]:
    """
    Uses Shaper TUI command:
      model.addImportSTEP(Part_doc, FileNameString, scalInterUnits, materials, colors) :contentReference[oaicite:1]{index=1}
    """
    model.begin()
    part_set = model.moduleDocument()
    part = model.addPart(part_set)
    part_doc = part.document()

    # scaleToUnits=True, materials=True, colors=True (matches typical usage)
    imp = model.addImportSTEP(part_doc, step_path, False, False, False)  # :contentReference[oaicite:2]{index=2}
    model.do()
    model.end()
    return part_set, part_doc, imp


# -------------------- Object resolution (names from CSV) --------------------

def index_named_results(part_doc: Any) -> Dict[str, List[Any]]:
    """
    Scans all objects in the document (may include hidden/history objects).
    ModelAPI_Document::allObjects() exists. :contentReference[oaicite:3]{index=3}
    """
    idx: Dict[str, List[Any]] = {}
    for obj in part_doc.allObjects():  # :contentReference[oaicite:4]{index=4}
        name = _wname(obj).strip()
        if name:
            idx.setdefault(name, []).append(obj)
    return idx

def resolve_targets_to_objects_or_error(
    targets: List[ComponentSpec],
    name_index: Dict[str, List[Any]],
) -> Dict[str, List[Any]]:
    """
    CSV Subcomponent label:
      - exact match
      - suffix '***' means prefix match (same as your FreeCAD script)
    If nothing matches => error (as requested).
    """
    resolved: Dict[str, List[Any]] = {}
    missing: List[str] = []

    for t in targets:
        key = (t.subcomponent or "").strip()
        if not key:
            missing.append(key)
            resolved[key] = []
            continue

        matches: List[Any] = []
        if key.endswith("***"):
            prefix = key[:-3]
            for nm, objs in name_index.items():
                if nm.startswith(prefix):
                    matches.extend(objs)
        else:
            matches.extend(name_index.get(key, []))

        # de-dup by object identity
        seen: Set[int] = set()
        uniq: List[Any] = []
        for o in matches:
            oid = id(o)
            if oid not in seen:
                seen.add(oid)
                uniq.append(o)

        resolved[key] = uniq
        if not uniq:
            missing.append(key)

    if missing:
        missing_sorted = sorted(set(missing))
        raise RuntimeError(
            "Some CSV mesh targets were not found by name in the imported STEP results.\n"
            "Missing targets:\n  - " + "\n  - ".join(missing_sorted)
        )

    return resolved


# -------------------- Slicing planes (axis-aware) --------------------

def compute_plane_normal(angle_deg: float, axis: str) -> Tuple[float, float, float]:
    """
    Matches your FreeCAD convention:
      axis=Y => (sin(a), 0, -cos(a))
      axis=Z => (-sin(a), cos(a), 0)
      axis=X => (0, -sin(a), cos(a))
    """
    a = math.radians(float(angle_deg))
    ax = axis.upper().strip()
    if ax == "Y":
        return (math.sin(a), 0.0, -math.cos(a))
    if ax == "Z":
        return (-math.sin(a), math.cos(a), 0.0)
    if ax == "X":
        return (0.0, -math.sin(a), math.cos(a))
    raise ValueError(f"axis must be 'X','Y','Z' (got {axis!r})")

def make_square_plane_face(
    angle_deg: float,
    axis: str,
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    size: float = PLANE_FACE_SIZE,
):
    """
    Creates a square planar face by center + normal + size:
      GeomAlgoAPI_FaceBuilder::squareFace(center, normal, size) :contentReference[oaicite:5]{index=5}
    Direction constructed by GeomAPI_Dir(x,y,z). :contentReference[oaicite:6]{index=6}
    """
    ox, oy, oz = origin
    nx, ny, nz = compute_plane_normal(angle_deg, axis)

    center = GeomAPI_Pnt(float(ox), float(oy), float(oz))
    normal = GeomAPI_Dir(float(nx), float(ny), float(nz))  # :contentReference[oaicite:7]{index=7}
    face = GeomAlgoAPI_FaceBuilder.squareFace(center, normal, float(size))  # :contentReference[oaicite:8]{index=8}
    return face


# -------------------- Section + connectivity (no tolerances) --------------------

def _shape_of_result(obj: Any):
    """
    Many Shaper results expose .shape() returning GeomAPI_Shape.
    We keep this best-effort because the exact Python wrapper type depends on what allObjects() yields.
    """
    for attr in ("shape", "result", "defaultResult"):
        try:
            m = getattr(obj, attr, None)
            if callable(m):
                candidate = m()
                if candidate is not None and hasattr(candidate, "shape") and callable(candidate.shape):
                    return candidate.shape()
                if candidate is not None and hasattr(candidate, "isNull") and callable(candidate.isNull):
                    return candidate
        except Exception:
            pass
    try:
        if hasattr(obj, "isNull") and callable(obj.isNull):
            return obj
    except Exception:
        pass
    return None

def intersect_shape_with_face(shape: Any, face: Any):
    """
    Uses GeomAPI_Shape::intersect(otherShape) (method is on GeomAPI_Shape; listed among shape methods). :contentReference[oaicite:9]{index=9}
    """
    try:
        return shape.intersect(face)  # :contentReference[oaicite:10]{index=10}
    except Exception:
        return None

def explore_edges(section_shape: Any) -> List[Any]:
    """
    Iterate edges using GeomAPI_ShapeExplorer. (User-provided acceptable API.)
    """
    edges: List[Any] = []
    try:
        exp = GeomAPI_ShapeExplorer(section_shape, GeomAPI_ShapeExplorer.EDGE)
        while exp.more():
            edges.append(exp.current())
            exp.next()
    except Exception:
        # Fallback: older wrappers may use different enum access
        exp = GeomAPI_ShapeExplorer(section_shape, "EDGE")
        while exp.more():
            edges.append(exp.current())
            exp.next()
    return edges

def vertex_point_xyz(vtx: Any) -> Tuple[float, float, float]:
    """
    Vertex -> point coordinates.
    GeomAPI_Vertex provides point(). :contentReference[oaicite:11]{index=11}
    GeomAPI_XYZ provides x(), y(), z(). :contentReference[oaicite:12]{index=12}
    """
    p = vtx.point()  # :contentReference[oaicite:13]{index=13}
    xyz = p.xyz()
    return (float(xyz.x()), float(xyz.y()), float(xyz.z()))  # :contentReference[oaicite:14]{index=14}

def edge_end_vertices(edge: Any) -> Tuple[Any, Any]:
    """
    Get start/end vertices of an edge using:
      GeomAlgoAPI_ShapeTools.findBounds(shape, v1, v2) :contentReference[oaicite:15]{index=15}
    """
    v1 = None
    v2 = None
    GeomAlgoAPI_ShapeTools.findBounds(edge, v1, v2)  # :contentReference[oaicite:16]{index=16}
    return v1, v2

def build_vertex_edge_adjacency(edges: List[Any]) -> Tuple[Dict[int, Any], Dict[int, List[int]], List[Tuple[int, int]]]:
    """
    Returns:
      - vertices_by_id: id(vertex) -> vertex
      - incident_edges: id(vertex) -> list of edge indices
      - edge_ends: list[(id(v_start), id(v_end))] aligned with edges list
    """
    vertices_by_id: Dict[int, Any] = {}
    incident_edges: Dict[int, List[int]] = {}
    edge_ends: List[Tuple[int, int]] = []

    for ei, e in enumerate(edges):
        v1, v2 = edge_end_vertices(e)
        if v1 is None or v2 is None:
            continue
        i1 = id(v1)
        i2 = id(v2)
        vertices_by_id[i1] = v1
        vertices_by_id[i2] = v2
        incident_edges.setdefault(i1, []).append(ei)
        incident_edges.setdefault(i2, []).append(ei)
        edge_ends.append((i1, i2))

    return vertices_by_id, incident_edges, edge_ends

def chains_from_edges(edges: List[Any]) -> List[List[Dict[str, float]]]:
    """
    Produces vertex-chains by walking edge connectivity via shared vertices.
    Each chain is returned as a list of vertex coordinate dicts: {"x":..,"y":..,"z":..}
    """
    if not edges:
        return []

    vertices_by_id, incident, ends = build_vertex_edge_adjacency(edges)
    if not ends:
        return []

    used_edges: Set[int] = set()
    chains: List[List[Dict[str, float]]] = []

    # Build degree map (topological)
    degree: Dict[int, int] = {vid: len(eis) for vid, eis in incident.items()}

    def other_end(ei: int, vid: int) -> Optional[int]:
        a, b = ends[ei]
        if vid == a:
            return b
        if vid == b:
            return a
        return None

    # Start vertices: degree 1 preferred (open chains), else any (loops)
    start_vertices = [vid for vid, deg in degree.items() if deg == 1]
    if not start_vertices:
        start_vertices = list(vertices_by_id.keys())

    for sv in start_vertices:
        # If all incident edges already used, skip
        if all(ei in used_edges for ei in incident.get(sv, [])):
            continue

        chain_vids: List[int] = [sv]
        cur_v = sv

        while True:
            # pick an unused incident edge
            cand = [ei for ei in incident.get(cur_v, []) if ei not in used_edges]
            if not cand:
                break

            ei = cand[0]
            used_edges.add(ei)

            nxt = other_end(ei, cur_v)
            if nxt is None:
                break

            chain_vids.append(nxt)
            cur_v = nxt

            # stop on open end
            if degree.get(cur_v, 0) == 1:
                break

            # stop if looped back
            if cur_v == sv:
                break

        # Convert to coordinates
        chain_pts: List[Dict[str, float]] = []
        for vid in chain_vids:
            vtx = vertices_by_id.get(vid)
            if vtx is None:
                continue
            x, y, z = vertex_point_xyz(vtx)
            chain_pts.append({"x": x, "y": y, "z": z})

        if len(chain_pts) >= 2:
            chains.append(chain_pts)

    # Any remaining loop edges not covered by start-vertex pass
    for ei in range(len(ends)):
        if ei in used_edges:
            continue
        # start from one end and walk until closure
        a, _b = ends[ei]
        sv = a
        chain_vids = [sv]
        cur_v = sv
        while True:
            cand = [cei for cei in incident.get(cur_v, []) if cei not in used_edges]
            if not cand:
                break
            cei = cand[0]
            used_edges.add(cei)
            nxt = other_end(cei, cur_v)
            if nxt is None:
                break
            chain_vids.append(nxt)
            cur_v = nxt
            if cur_v == sv:
                break

        chain_pts: List[Dict[str, float]] = []
        for vid in chain_vids:
            vtx = vertices_by_id.get(vid)
            if vtx is None:
                continue
            x, y, z = vertex_point_xyz(vtx)
            chain_pts.append({"x": x, "y": y, "z": z})
        if len(chain_pts) >= 2:
            chains.append(chain_pts)

    return chains


# -------------------- Main workflow --------------------

def run(
    step_path: str = STEP_PATH,
    csv_path: str = CSV_PATH,
    out_base_dir: str = OUT_BASE_DIR,
    angles_deg: Optional[List[float]] = None,
    axis: str = DEFAULT_AXIS,
) -> None:
    if angles_deg is None:
        angles_deg = list(DEFAULT_ANGLES_DEG)

    # Load CSV
    print("Loading CSV")
    specs = load_component_csv(csv_path)
    targets = mesh_targets(specs)

    # Import STEP
    print("Importing Step")
    _part_set, part_doc, _imp = import_step_into_new_part(step_path)

    # Index names and resolve targets (error if missing)
    print("Indexing Names")
    idx = index_named_results(part_doc)
    resolved = resolve_targets_to_objects_or_error(targets, idx)

    # Output folder structure: out/pts_out/<BC>/<Subcomponent>/
    pts_root = os.path.join(out_base_dir, "pts_out")
    _ensure_dir(pts_root)

    print("Building Lookup")
    # Build a quick lookup subcomponent->BC
    bc_by_sub: Dict[str, str] = {}
    for t in targets:
        bc_by_sub[t.subcomponent] = t.bc or "NA"

    # Process each CSV mesh label
    print("Processing Each CSV label")
    for sub_label, objs in resolved.items():
        print(sub_label)
        bc = bc_by_sub.get(sub_label, "NA") or "NA"

        comp_dir = os.path.join(pts_root, _sanitize_filename(bc), _sanitize_filename(sub_label))
        _ensure_dir(comp_dir)

        # Slice each angle and write one JSON per angle (simple, explicit IO)
        print("Angle Loop")
        for ang in angles_deg:
            plane_face = make_square_plane_face(angle_deg=float(ang), axis=axis, origin=(0.0, 0.0, 0.0))

            instances_out: List[Dict[str, Any]] = []
            print("Iterating Objects")
            for obj in objs:
                obj_name = _wname(obj)
                shape = _shape_of_result(obj)
                if shape is None:
                    print("Shape Is None")
                    instances_out.append({
                        "object_name": obj_name,
                        "error": {"type": "NoShape", "message": "Could not retrieve GeomAPI_Shape from resolved object."},
                        "chains": [],
                    })
                    continue
                print("Computing Intersect")
                sec = intersect_shape_with_face(shape, plane_face)
                if sec is None:
                    print("Sec Is None")
                    instances_out.append({
                        "object_name": obj_name,
                        "error": {"type": "IntersectionFailed", "message": "Shape.intersect(plane_face) returned None/failed."},
                        "chains": [],
                    })
                    continue

                try:
                    print("Exploring Edges")
                    edges = explore_edges(sec)
                    print("Making chains from Edges")
                    chains = chains_from_edges(edges)
                    instances_out.append({
                        "object_name": obj_name,
                        "error": None,
                        "chains": chains,
                        "num_edges": len(edges),
                        "num_chains": len(chains),
                    })
                except Exception as ex:
                    print("Error on pts list creation")
                    instances_out.append({
                        "object_name": obj_name,
                        "error": {"type": type(ex).__name__, "message": str(ex)},
                        "chains": [],
                    })

            payload = {
                "step_path": step_path,
                "csv_path": csv_path,
                "component": sub_label,
                "bc": bc,
                "axis": axis,
                "angle_deg": float(ang),
                "instances": instances_out,
            }
            print("Writing Payload")
            out_path = os.path.join(comp_dir, f"slice_{int(round(float(ang))):03d}.json")
            _write_json(out_path, payload)

    print("Done.")
    print(f"Wrote outputs under: {os.path.join(out_base_dir, 'pts_out')}")

run()
