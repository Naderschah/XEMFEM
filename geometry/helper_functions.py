import re
import os
import sys

from salome.shaper import model
from SketchAPI import SketchAPI_Circle, SketchAPI_Point, SketchAPI_Arc, SketchAPI_Line
import ModelAPI
from salome.shaper import geom

# Post Partition identification
from GeomAlgoAPI import GeomAlgoAPI_ShapeTools
from GeomAPI import GeomAPI_ShapeExplorer 
import SHAPERSTUDY
from salome.smesh import smeshBuilder
import SMESH

# ------------------------------- Sketching Functions ----------------------------------------------
def make_anchors(Sketch, lines, pts, line_indices):
    """
    Anchors a geometry by attaching fixed construction points to line midpoints.
    Assumes straight lines and that each line connects two points in `pts`.

    Makes auxiliary (construction) points so they are *not* part of the solid geometry.
    """
    n_pts = len(pts)
    anchors = []
    for line_idx in line_indices:
        # indices of the endpoints of the line
        i0 = line_idx
        i1 = (line_idx + 1) % n_pts

        x_mid = 0.5 * (pts[i0][1] + pts[i1][1])
        y_mid = 0.5 * (pts[i0][2] + pts[i1][2])

        # Create anchor as construction geometry
        anchor = Sketch.addPoint(x_mid, y_mid)
        anchor.setAuxiliary(True)  # <-- Mark it as construction geometry

        # Coincident constraint pins it to the line
        Sketch.setCoincident(anchor, lines[line_idx].result())

        # Fix its position for solver stabilization
        Sketch.setFixed(anchor)

        anchors.append(anchor)

    return anchors

# ----------------------------- Selection Functions ---------------------------------
def get_all_sketch_features(Sketch):
    """Mainly usefull for fillets
    I have still not figured out how to select the produced edges
    """
    feat_list = Sketch.features()
    ids = set()
    for i in range(feat_list.size()):
        obj = feat_list.object(i)          # ModelAPI_Object
        ids.add(obj)
    return ids

def filter_for_fillet_edges(before_ids, after_ids):
    """Filter for fillet result edges"""
    new_feats = []
    for obj in after_ids:
        if obj in before_ids:
            continue
        feat = ModelAPI.objectToFeature(obj)
        if feat:
            new_feats.append(feat)

    fillet_arc_features = [feat for feat in new_feats if feat.getKind() == "SketchArc"]
    return fillet_arc_features


def face_feature_to_selections(src):
    """
    Turn any 'face-like' object into a list of FACE selections.

    src can be:
      - ModelHighAPI_Selection       (already a selection)
      - feature with .results()      (Cut, Group, etc., possibly multi-volume)
      - feature with .result()       (Face with single result)
      - ModelAPI_Result              (possibly multi-face via subResult)
    """

    # 1) Already a selection?
    if hasattr(src, "type") and callable(getattr(src, "type", None)):
        # Assume it's already a FACE selection; if you need, check src.type() == "FACE"
        return [src]

    # 2) Feature with multiple results (e.g. Cut)
    if hasattr(src, "results"):
        res_list = list(src.results())

    # 3) Feature with single result (e.g. addFace)
    elif hasattr(src, "result"):
        res_list = [src.result()]

    else:
        # 4) Assume it's already a Result
        res_list = [src]

    sels = []
    for r in res_list:
        # r is a Result; check for sub-faces
        if hasattr(r, "numberOfSubs") and r.numberOfSubs() > 0:
            for i in range(r.numberOfSubs()):
                sub_sel = r.subResult(i)  # This is already a ModelHighAPI_Selection
                sels.append(sub_sel)
        else:
            # Single face under this result
            sels.append(model.selection("FACE", r.name()))
    return sels

def safe_face_selections(feature, base_name):
    """Selects all faces associated with this object and names them"""
    res = face_feature_to_selections(feature) if feature is not None else []
    for idx, sub in enumerate(feature.results(), start=1):
            sub.setName(f"{base_name}_Face{idx}")
    return res

def safe_face_selections_for_dict(face_dict):
    """Do this for all elements in the dictionary"""
    res = {}
    for key, value in face_dict.items():
        res[key] = safe_face_selections(value, key)
    return res

def makeGroupBySelection(document, selection, baseName):
    # Works but finniky if the selection inputs are off
    
    faces = model.addGroup(document,selection)
    faces.setName(baseName+"Faces")
    faces.result().setName(baseName+"Faces")
    group = model.addGroupShape(document,[model.selection("COMPOUND", baseName+"Faces")])
    group.setName(baseName+"Group")
    group.result().setName(baseName+"Group")
    return group

def makeGroupByName(document, baseName, face_dict):
    # To select all faces created by sketch result <- Used when i didnt explicitly track each face
    #faces = model.addGroup(Cryostat_doc,[model.selection("FACE", "all-in-"+baseName) for baseName in face_dict.keys()])
    # To select by face name directly
    faces = model.addGroup(document,[model.selection("FACE", baseName) for baseName in face_dict.keys()])
    faces.setName(baseName+"Faces")
    faces.result().setName(baseName+"Faces")
    group = model.addGroupShape(document,[model.selection("COMPOUND", baseName+"Faces")])
    group.setName(baseName+"Group")
    group.result().setName(baseName+"Group")
    return group


# -------------------------------- Dict Parsing ------------------------------------
def check_circle_params(component):
    return "Radius" in component

def check_box_params(component):
    return "Height" in component and "Width" in component

def check_fillet_params(component):
    return "FilletRadius" in component

def check_repetition(component):
    return "HorizontalPitch" in component or "VerticalPitch" in component

def check_hull(component):
    return bool(component.get("hull", False))

def make_box(component):
    """Populates the pts field in component for boxes"""
    x0 = component["RadialPosition"]
    y0 = component["VerticalPosition"]
    W = component["Width"]
    H = component["Height"]

    pts = [
        [x0,     y0],
        [x0,     y0 + H],
        [x0 + W, y0 + H],
        [x0 + W, y0],
    ]

    component["pts"] = pts

def make_shape(Sketch, component, name="<component>"):
    """
    component["pts"]: list of N vertices around boundary.
      Each vertex entry:
        ["line", x, y]                -> segment from this vertex to next is a line
        ["arc",  x, y, cx, cy]        -> segment from this vertex to next is an arc with center (cx, cy)

    Returns list of created sketch objects (lines/arcs) forming a closed loop.
    """
    pts = component["pts"]
    n = len(pts)
    lines = []

    for i in range(n):
        x1, y1 = pts[i][1], pts[i][2]
        x2, y2 = pts[(i + 1) % n][1], pts[(i + 1) % n][2]

        if pts[i][0] == "line":
            l = Sketch.addLine(x1, y1, x2, y2)

        elif pts[i][0] == "arc":
            if len(pts[i]) == 5:
                cx, cy = pts[i][3], pts[i][4]
                direction = True
            elif len(pts[i]) == 6:
                cx, cy, direction = pts[i][3], pts[i][4], pts[i][5]
            else:
                raise Exception("pts for arc invalid length in component " + str(name))

            l = Sketch.addArc(
                cx, cy,
                x1, y1,
                x2, y2,
                direction
            )

        else:
            raise Exception("Argument " + str(pts[i][0]) + " not recognized for line creation in component " + str(name))

        lines.append(l)

    for i in range(n):
        Sketch.setCoincident(
            lines[i].endPoint(),
            lines[(i + 1) % n].startPoint()
        )

    return lines

def make_circle(Sketch, component):
    circle = Sketch.addCircle(
        component["RadialPosition"],     # x
        component["VerticalPosition"],   # y
        component["Radius"]
    )
    circle_vertex, circle_edges = circle.results()
    return [circle_edges]

def repeat_object(Sketch, component, objects):
    # Inclusive of original object
    count = component["Number"]

    HorizontalPitch = component.get("HorizontalPitch", 0.0)
    VerticalPitch = component.get("VerticalPitch", 0.0)

    # Total offset from first to last instance
    x_offset = HorizontalPitch * count
    y_offset = VerticalPitch * count

    x0 = component.get("RadialPosition", 0.0)
    y0 = component.get("VerticalPosition", 0.0)

    start_pt = Sketch.addPoint(x0, y0)
    end_pt = Sketch.addPoint(x0 + x_offset / count, y0 + y_offset / count)
    start_pt.setAuxiliary(True)
    end_pt.setAuxiliary(True)
    multi = Sketch.addTranslation(
        objects,
        start_pt.coordinates(),
        end_pt.coordinates(),
        count
    )
    
    return multi

def make_fillets(Sketch, component, lines):
    prior_ids = get_all_sketch_features(Sketch)

    if "FilletIndeces" in component:
        fillet_indices = component["FilletIndeces"]
    else:
        fillet_indices = list(range(len(lines)))

    n = len(lines)

    # During fillets objects like to move (anchor a few edges by midpoint)
    # NOTE: make_anchors expects the original point list, not len(lines).
    # If you don't have pts here, anchoring by line indices still works if make_anchors
    # is updated accordingly; otherwise remove anchors.
    if "pts" in component:
        _ = make_anchors(Sketch, lines, component["pts"], [0, 1, 2])

    for i in fillet_indices:
        vertex = lines[(i - 1) % n].endPoint()
        Sketch.setFilletWithRadius(vertex, component["FilletRadius"])

    after_ids = get_all_sketch_features(Sketch)
    fillet_edges = filter_for_fillet_edges(prior_ids, after_ids)
    return fillet_edges


def validate_pts(component):
    for idx, item in enumerate(component["pts"]):
        if not isinstance(item[0], str):
            component["pts"][idx].insert(0, "line")

def geom_sketcher(
    component,
    part_doc,
    name="<component>",
    sketch=None,
    shrink_below_y=None,
    shrinkage_factor=1.0,
    apply_shrinkage=False,
):
    """
    Optional shrinkage:
      If apply_shrinkage is True and shrink_below_y is not None:
        for any point/arc endpoint/arc center with y < shrink_below_y:
            y := y * shrinkage_factor
    Shrinkage is applied ONLY to the geometry being drawn (local copy), not stored back into component.
    """
    # Check what kind of structure is to be built
    has_box_params = check_box_params(component)
    has_circle_params = check_circle_params(component)
    has_fillet_params = check_fillet_params(component)
    has_repeat = check_repetition(component)
    is_hull = check_hull(component)

    if has_box_params and has_circle_params:
        raise Exception("Shape recognized as circle and box! " + str(name))
    if has_circle_params and has_fillet_params:
        raise Exception("Circle can not have fillet! " + str(name))

    # Make a sketch
    if sketch is None:
        sketch = model.addSketch(part_doc, model.defaultPlane("XOY"))
        sketch.setName(str(name))

    # TODO Just rename the things
    Sketch = sketch

    # Populate pts for boxes
    if has_box_params:
        make_box(component)

    # Helper: make a shrinked copy of pts (does not mutate component)
    def _maybe_shrink_pts(src_pts):
        if (not apply_shrinkage) or (shrink_below_y is None) or (shrinkage_factor == 1.0):
            return src_pts

        out = []
        for item in src_pts:
            # item is either ["line", x, y] or ["arc", x, y, cx, cy]
            if item[0] == "line":
                x, y = item[1], item[2]
                if y < shrink_below_y:
                    y = y * shrinkage_factor
                out.append(["line", x, y])

            elif item[0] == "arc":
                if len(item) != 5:
                    print(item)
                    raise Exception("pts for arc invalid length in component " + str(name))
                x, y, cx, cy = item[1], item[2], item[3], item[4]

                if y < shrink_below_y:
                    y = y * shrinkage_factor
                if cy < shrink_below_y:
                    cy = cy * shrinkage_factor

                out.append(["arc", x, y, cx, cy])
            else:
                raise Exception("Argument " + str(item[0]) + " not recognized for line creation in component " + str(name))
        return out

    # Draw geometry unless this node is a hull-only container
    lines = []
    if not is_hull:
        if has_circle_params:
            # Circle is defined by center (x,y) and radius. Shrinking only in y is ambiguous for a circle;
            # leave unchanged unless you later decide to convert to an ellipse explicitly.
            lines = make_circle(Sketch, component)
        else:
            validate_pts(component)
            local_component = component
            if apply_shrinkage and "pts" in component:
                local_component = dict(component)
                local_component["pts"] = _maybe_shrink_pts(component["pts"])
            lines = make_shape(Sketch, local_component, name=name)

        if has_fillet_params and lines:
            filleted = make_fillets(Sketch, component, lines)
            lines = list(lines) + list(filleted)

        if has_repeat and lines:
            repeat_object(Sketch, component, lines)

    # Recurse into sub-sketches
    sub = component.get("sub_sketches", {})
    if isinstance(sub, dict):
        for child_name, child_component in sub.items():
            Sketch = geom_sketcher(
                child_component,
                part_doc,
                name=child_name,
                sketch=Sketch,
                shrink_below_y=shrink_below_y,
                shrinkage_factor=shrinkage_factor,
                apply_shrinkage=apply_shrinkage,
            )
    elif isinstance(sub, list):
        for idx, child_component in enumerate(sub):
            Sketch = geom_sketcher(
                child_component,
                part_doc,
                name=f"{name}_sub{idx}",
                sketch=Sketch,
                shrink_below_y=shrink_below_y,
                shrinkage_factor=shrinkage_factor,
                apply_shrinkage=apply_shrinkage,
            )


    return Sketch


def sketch_from_dict(
    sketch_dict,
    part_doc,
    makeface = True,
    apply_shrinkage=False,
    shrink_below_y=None,
    shrinkage_factor=1.0,
):
    """
    sketch_dict: dict-of-dicts describing geometry
    Returns: dict mapping component name -> Sketch or Face
    """
    sketch_map = {}
    face_map = {}

    for name, component in sketch_dict.items():
        if name == 'shrinkage_factor':
            continue
        sk = geom_sketcher(
            component,
            part_doc,
            name=name,
            sketch=None,
            apply_shrinkage=apply_shrinkage,
            shrink_below_y=shrink_below_y,
            shrinkage_factor=shrinkage_factor,
        )
        if makeface:
            model.do()
            face = model.addFace(part_doc, [sk.result()])
            face.setName(name)
            for idx, face_i in enumerate(face.results()):
                name_sub = name+str(idx)
                face_map[name_sub] = face_i
                face_i.setName(name_sub) 
            sketch_map[name] = face
        else:
            sketch_map[name] = sk

    return sketch_map, face_map
# -------------------------------- Misc --------------------------------------
def flatten_dict(sel_dict):
    out = []
    for key in sel_dict.keys():
        out.append(sel_dict[key])
    return out



# ---------------------------- Post Partition Naming -------------------------
def get_area(FaceSelection, doc):
    name   = FaceSelection.name() 
    sel = model.selection("FACE", name)
    _, area, _ = model.getGeometryCalculation(doc, sel)
    return area 
def center_of_weight(part_doc, face_selection):
    # Somehow this creates the center of mass point - I do not get why it works but it does 
    point = model.addPoint(part_doc, face_selection)
    coords = model.getPointCoordinates(
        part_doc,
        model.selection("VERTEX", point.name())
    )
    return coords[0], coords[1] # coords is a list [x, y, z]
def match_and_rename_partition_faces(partition, pre_partition, post_partition,
                                     dist_tol=1e-4):
    """
    For each post-partition face, find the pre-partition face whose center of
    weight is closest (within dist_tol) and rename the post face to
    pre_name + "_part".

    partition: the Partition feature
    pre_partition: dict {pre_name: [area, cx, cy, cz]}
    post_partition: dict {post_name: [area, cx, cy, cz]}
    dist_tol: maximum allowed distance between centers (in model units)
    """
    # Map current post face names to their indices
    n_post = partition.result().numberOfSubs()
    post_index_by_name = {
        partition.result().subResult(i).name(): i
        for i in range(n_post)
    }
    renamed_count = 0
    names = []
    for post_name, (area_p, px, py) in post_partition.items():
        best_pre_name = None
        best_dist = None
        for pre_name, (area_pre, cx, cy) in pre_partition.items():
            dx = px - cx
            dy = py - cy
            d = math.sqrt(dx*dx + dy*dy)
            if best_dist is None or d < best_dist:
                best_dist = d
                best_pre_name = pre_name
        # If closest pre-face is within tolerance, rename post-face
        if best_pre_name is not None and best_dist is not None and best_dist <= dist_tol:
            idx = post_index_by_name[post_name]
            partition.result().subResult(idx).setName(best_pre_name + "_part")
            renamed_count += 1
            names.append(best_pre_name + "_part")

    print("Successfully renamed {} of {} post-partition faces".format(
        renamed_count, len(post_partition)
    ))
    return names, renamed_count

### BBox logic ---------------------------------------------
def _validate_rec(name, rec, debug=False, role=""):
    """
    rec must be a dict with keys: area, cx, cy, bbox
    bbox = (xmin, ymin, xmax, ymax)
    """
    if not debug:
        return True

    ok = True
    msgs = []

    if rec is None:
        msgs.append("rec is None")
        ok = False
    elif not isinstance(rec, dict):
        msgs.append(f"rec is not dict: {type(rec)}")
        ok = False
    else:
        for k in ("area", "cx", "cy", "bbox"):
            if k not in rec:
                msgs.append(f"missing key '{k}'")
                ok = False

        if "area" in rec:
            a = rec["area"]
            if not isinstance(a, (int, float)) or a <= 0:
                msgs.append(f"invalid area: {a}")
                ok = False

        if "cx" in rec and "cy" in rec:
            for k in ("cx", "cy"):
                v = rec[k]
                if not isinstance(v, (int, float)) or not math.isfinite(v):
                    msgs.append(f"invalid {k}: {v}")
                    ok = False

        if "bbox" in rec:
            bb = rec["bbox"]
            if (
                not isinstance(bb, (tuple, list)) or
                len(bb) != 4
            ):
                msgs.append(f"invalid bbox type/len: {bb}")
                ok = False
            else:
                xmin, ymin, xmax, ymax = bb
                if not all(isinstance(v, (int, float)) for v in bb):
                    msgs.append(f"bbox has non-numeric entries: {bb}")
                    ok = False
                if xmax <= xmin or ymax <= ymin:
                    msgs.append(f"degenerate bbox: {bb}")
                    ok = False

    if not ok:
        print(f"[DEBUG][INVALID {role}] {name}")
        for m in msgs:
            print("   -", m)

    return ok

def _validate_probe(x, y, bbox, debug=False, label=""):
    if not debug:
        return True

    xmin, ymin, xmax, ymax = bbox
    ok = True

    if not all(math.isfinite(v) for v in (x, y)):
        print(f"[DEBUG][BAD PROBE {label}] non-finite: ({x},{y})")
        ok = False

    if not (xmin <= x <= xmax and ymin <= y <= ymax):
        print(f"[DEBUG][BAD PROBE {label}] outside bbox: ({x},{y}) not in {bbox}")
        ok = False

    return ok


from GeomAlgoAPI import GeomAlgoAPI_ShapeTools
from GeomAPI import GeomAPI_ShapeExplorer, GeomAPI_Shape


def bbox_from_face_name_xy(part_doc, face_result_name):
    """
    Returns (xmin, ymin, xmax, ymax) computed numerically from the vertices of the face (XY only).
    Uses filters: BelongsTo -> select("Vertex")
    """

    face_sel = model.selection("FACE", face_result_name)

    flt = model.filters(part_doc, [
        model.addFilter(name="BelongsTo", args=[face_sel])
    ])
    v_sels = flt.select("Vertex")  # list of vertex selections :contentReference[oaicite:2]{index=2}

    if not v_sels:
        raise RuntimeError(f"No vertices selected for face {face_result_name}")

    # init from first vertex
    x0, y0, _ = model.getPointCoordinates(part_doc, v_sels[0])
    xmin = xmax = x0
    ymin = ymax = y0

    for vs in v_sels[1:]:
        x, y, _ = model.getPointCoordinates(part_doc, vs)
        if x < xmin: xmin = x
        if x > xmax: xmax = x
        if y < ymin: ymin = y
        if y > ymax: ymax = y

    return xmin, ymin, xmax, ymax
from math import hypot
import math
import random


def _bbox_area(b):
    xmin, ymin, xmax, ymax = b
    return (xmax - xmin) * (ymax - ymin)

def _contains(b_outer, b_inner, eps=0.0):
    ox0, oy0, ox1, oy1 = b_outer
    ix0, iy0, ix1, iy1 = b_inner
    return (ox0 <= ix0 + eps and oy0 <= iy0 + eps and
            ox1 >= ix1 - eps and oy1 >= iy1 - eps)
def _to_rec(vals):
    # vals: [area, cx, cy, (xmin,ymin,xmax,ymax)]
    return {"area": vals[0], "cx": vals[1], "cy": vals[2], "bbox": vals[3]}

def _bbox_area(b):
    xmin, ymin, xmax, ymax = b
    return (xmax - xmin) * (ymax - ymin)


DEBUG_FACE = True  # set to a post face name string to restrict debug, or None for all
DEBUG_PRE  = True  # set to a pre face name string to restrict debug, or None for all

def _dbg_enabled(post_name, pre_name):
    if DEBUG_FACE is not None and post_name != DEBUG_FACE:
        return False
    if DEBUG_PRE is not None and pre_name != DEBUG_PRE:
        return False
    return True

def faces_containing_point(part_doc, x, y, z=0.0):
    p = model.addPoint(part_doc, float(x), float(y), float(z))
    model.do()
    vsel = model.selection("VERTEX", "all-in-" + p.name())

    flt = model.filters(part_doc, [model.addFilter(name="BelongsTo", args=[vsel])])
    face_sels = flt.select("Face")

    out = []
    for fs in face_sels:
        ctx_res, subshape = fs.resultSubShapePair()
        out.append(ctx_res.name())
    return out  # keep order for readable prints


def deterministic_probe_points(bbox, margin=1e-9):
    xmin, ymin, xmax, ymax = bbox
    xmin += margin; ymin += margin; xmax -= margin; ymax -= margin
    mx, my = 0.5*(xmin+xmax), 0.5*(ymin+ymax)

    yield ("center", mx, my)
    yield ("BL", 0.5*(xmin+mx), 0.5*(ymin+my))
    yield ("BR", 0.5*(mx+xmax), 0.5*(ymin+my))
    yield ("TL", 0.5*(xmin+mx), 0.5*(my+ymax))
    yield ("TR", 0.5*(mx+xmax), 0.5*(my+ymax))


def random_probe_points(bbox, n, margin=1e-9):
    xmin, ymin, xmax, ymax = bbox
    xmin += margin; ymin += margin; xmax -= margin; ymax -= margin
    for k in range(n):
        yield (f"R{k}", random.uniform(xmin, xmax), random.uniform(ymin, ymax))


def passes_containment_sampling(part_doc, post_name, pre_name, bbox,
                                z=0.0, need_successes=4, rand_tries=200, margin=1e-9):
    """
    Debug prints focus on the condition checks:
      - what faces SHAPER reports at each probe point
      - membership tests for post_name and pre_name
    """
    successes = 0

    def test_point(tag, x, y):
        nonlocal successes

        if _dbg_enabled(post_name, pre_name):
            _validate_probe(tag, x, y, z, bbox, debug=True)

        faces = faces_containing_point(part_doc, x, y, z=z)

        in_post = post_name in faces
        in_pre  = pre_name in faces

        if _dbg_enabled(post_name, pre_name):
            print(f"[probe {tag}] point=({x:.6g},{y:.6g},{z:.6g})")
            print(f"  post='{post_name}' in_faces={in_post}")
            print(f"  pre ='{pre_name}' in_faces={in_pre}")
            print(f"  faces_returned({len(faces)}): {faces}")

        if in_post and in_pre:
            successes += 1
            if _dbg_enabled(post_name, pre_name):
                print(f"  SUCCESS {successes}/{need_successes}\n")
            return successes >= need_successes
        else:
            if _dbg_enabled(post_name, pre_name):
                print("  fail\n")
            return False

    if _dbg_enabled(post_name, pre_name):
        print(f"=== containment check: post='{post_name}' pre='{pre_name}' bbox={bbox} need={need_successes} z={z} ===")

    for tag, x, y in deterministic_probe_points(bbox, margin=margin):
        if test_point(tag, x, y):
            return True

    for tag, x, y in random_probe_points(bbox, rand_tries, margin=margin):
        if test_point(tag, x, y):
            return True

    if _dbg_enabled(post_name, pre_name):
        print(f"=== END: successes={successes}/{need_successes} (FAIL) ===\n")

    return False


def associate_by_bbox_smallest_container(pre_partition, post_partition, part_doc, eps=0.0,
                                        z=0.0, need_successes=4, rand_tries=200, margin=1e-9,
                                        debug=False):
    pre = {n: _to_rec(v) for n, v in pre_partition.items()}
    post = {n: _to_rec(v) for n, v in post_partition.items()}

    if debug:
        print("[DEBUG] validating pre_partition records")
        for n, r in pre.items():
            _validate_rec(n, r, debug=True, role="pre")

        print("[DEBUG] validating post_partition records")
        for n, r in post.items():
            _validate_rec(n, r, debug=True, role="post")

    pre_items = list(pre.items())

    mapping = {}
    unresolved = []

    for post_name, post_rec in post.items():
        pb = post_rec["bbox"]
        pcx, pcy = post_rec["cx"], post_rec["cy"]

        if debug and (DEBUG_FACE is None or post_name == DEBUG_FACE):
            print(f"\n### POST '{post_name}' bbox={pb} cog=({pcx},{pcy})")

        candidates = []
        for pre_name, pre_rec in pre_items:
            ok_bbox = _contains(pre_rec["bbox"], pb, eps=eps)
            if debug and _dbg_enabled(post_name, pre_name):
                print(f"  bbox? pre='{pre_name}' pre_bbox={pre_rec['bbox']} contains={ok_bbox}")
            if ok_bbox:
                candidates.append((pre_name, pre_rec))

        if not candidates:
            unresolved.append(post_name)
            if debug and (DEBUG_FACE is None or post_name == DEBUG_FACE):
                print(f"  -> no bbox candidates")
            continue

        filtered = []
        for pre_name, pre_rec in candidates:
            ok = passes_containment_sampling(
                part_doc, post_name, pre_name, pb,
                z=z, need_successes=need_successes, rand_tries=rand_tries, margin=margin
            )
            if ok:
                filtered.append((pre_name, pre_rec))

        if not filtered:
            unresolved.append(post_name)
            if debug and (DEBUG_FACE is None or post_name == DEBUG_FACE):
                print(f"  -> no candidates passed sampling containment")
            continue

        def key_fn(item):
            pre_name, pre_rec = item
            barea = _bbox_area(pre_rec["bbox"])
            dx = pre_rec["cx"] - pcx
            dy = pre_rec["cy"] - pcy
            dist = math.sqrt(dx*dx + dy*dy)
            return (barea, dist, pre_rec["area"], pre_name)

        best_pre_name, _ = min(filtered, key=key_fn)
        mapping[post_name] = best_pre_name

        if debug and (DEBUG_FACE is None or post_name == DEBUG_FACE):
            print(f"  -> PICK '{best_pre_name}'")

    return mapping, unresolved
## Manual Renaming
def rename_partition_subresults_manual(partition, idx_to_base):
    """
    Renames partition.result().subResult(idx) using:
      <baseName>_part_<n>_manual

    Ensures uniqueness:
      - n is counted per baseName within idx_to_base
      - also avoids collisions with existing names in the partition result set
    """
    res = partition.result()
    n_sub = res.numberOfSubs()  # if unavailable in your build, replace with a known count

    existing = set(res.subResult(i).name() for i in range(n_sub))

    per_base = {}  # baseName -> next n (starting at 1)
    new_names = []

    for idx, base in idx_to_base.items():
        idx = int(idx)
        sr = res.subResult(idx)

        if 'partition' not in sr.name().lower():
            raise Exception("Element with idx " + str(idx) + " already has name " + sr.name() + " interrupted renaming to " + base)


        n = per_base.get(base, 0) + 1
        per_base[base] = n

        # propose name and bump n until it is unique against existing + already assigned
        while True:
            candidate = f"{base}_part_{n}_manual"
            if candidate not in existing:
                break
            n += 1
            per_base[base] = n

        print(f"idx {idx}: '{sr.name()}' -> '{candidate}'")
        sr.setName(candidate)

        existing.add(candidate)
        new_names.append(candidate)

    model.do()
    return new_names


## Last Renaming strategy ----------------------------------------
# LXe and GXe
def rename_two_largest_partition_faces(partition, post_partition):
    """
    post_partition: dict {face_name: [area, cx, cy] or [area, cx, cy, ...]}
    Renames the largest-area face to 'LXe_part' and the second-largest to 'GXe_part'.
    """
    # sort by area descending
    top2 = sorted(post_partition.items(), key=lambda kv: kv[1][0], reverse=True)[:2]
    if len(top2) < 2:
        raise RuntimeError(f"Need at least 2 faces, got {len(top2)}")

    (name1, _), (name2, _) = top2

    pres = partition.result()
    idx_by_name = {pres.subResult(i).name(): i for i in range(pres.numberOfSubs())}

    pres.subResult(idx_by_name[name1]).setName("LXe_part")
    pres.subResult(idx_by_name[name2]).setName("GXe_part")

    return ["LXe_part", "GXe_part"], [name1, name2]  # new names, original names
def _assign_unique_name(base, used):
    """
    base like "GateInsulatingFrame0_part"
    If already used, increment trailing integer before "_part" (or append one) until unique.
    """
    if base not in used:
        used.add(base)
        return base

    if base.endswith("_part"):
        stem = base[:-5]
        # split trailing digits
        j = len(stem)
        while j > 0 and stem[j-1].isdigit():
            j -= 1
        prefix = stem[:j]
        num = stem[j:]
        start = int(num) if num else 0

        k = start + 1
        while True:
            cand = f"{prefix}{k}_part"
            if cand not in used:
                used.add(cand)
                return cand
            k += 1
    else:
        k = 1
        while True:
            cand = f"{base}{k}"
            if cand not in used:
                used.add(cand)
                return cand
            k += 1

