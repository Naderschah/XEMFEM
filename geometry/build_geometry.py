import os
import shutil
import sys
import yaml
import time 
import inspect
import salome
from pathlib import Path
import importlib.util
from ModelAPI import ModelAPI_Session, StringList

# ============== Paths for path unaware TUI =====================
base_path = os.environ.get("BASE_PATH", "/work/geometry/")
sys.path.insert(0, base_path)
from helper_functions import *

# ==============  Load Config File ===========================
with open(base_path + "config.yaml") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Resolve which geometry file to use
import glob

geom = config["mesh"]["geometry"]
geom_override = config["mesh"].get("geom_volt_path_overwrite")
is_tpc = config["mesh"].get("is_tpc", True)
mesh_path = config["mesh"]['path']
manual_only = config["mesh"].get("manual_naming_only", False)
mesh_path = os.path.normpath(
    os.path.expandvars(os.path.expanduser(mesh_path))
)
shrinkage_factor = config["mesh"].get('shrinkage_factor', 1-0.014)

# Optional external-boundary tagging (independent of is_tpc).
# Defaults preserve existing behavior.
mark_external_boundary = config["mesh"].get("mark_external_boundary", is_tpc)
external_boundary_name = config["mesh"].get("external_boundary_name", "BC_Cryostat")
split_axis_boundary = config["mesh"].get("split_axis_boundary", is_tpc)
axis_boundary_name = config["mesh"].get("axis_boundary_name", "BC_r0")
axis_boundary_tol = float(config["mesh"].get("axis_boundary_tol", 1e-9))

# If it looks like a file, strip filename
if '.' in mesh_path.split('/')[-1]:
    mesh_path = os.path.dirname(mesh_path)

if geom_override:
    # Allow relative paths, resolve against base_path
    geom_path = (
        geom_override
        if os.path.isabs(geom_override)
        else os.path.join(base_path, geom_override)
    )
else:
    geom_path = os.path.join(base_path, "geometries")
print("globbing: ", os.path.join(geom_path, f"*{geom}*.py"))
fname = glob.glob(
    os.path.join(geom_path, f"*{geom}*.py")
)

if len(fname) == 1:
    module_name = Path(fname[0]).stem # Strip .py
    spec = importlib.util.spec_from_file_location(module_name, fname[0])
    geometry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
else:
    raise Exception("Found " + str(len(fname)) + " choices for geometry source " + config['mesh']['geometry'])

slice_alignment = getattr(geometry, "_SLICE_ALIGNMENT", None)
component_store = getattr(geometry, "_COMPONENT_STORE", None)
if isinstance(slice_alignment, dict) and isinstance(component_store, dict):
    ref_name = slice_alignment.get("reference_component")
    if ref_name and ref_name not in component_store:
        print(
            "[slice_alignment] reference component missing for this slice; "
            "disabling runtime alignment:",
            ref_name,
        )
        geometry._SLICE_ALIGNMENT = {}

# ==============  Debug Options ===========================
makeface = True
stop_after_sketch = (
    config
    .get('debug', {})
    .get('StopAfterSketh', False)
)
save_after_sketch = (
    config
    .get('debug', {})
    .get('SaveAfterSketch', False)
)
save_on_partition_failure = (
    config
    .get('debug', {})
    .get('SaveOnPartitionFailure', save_after_sketch)
)


def _study_save_path(tag: str) -> str:
    base_name = "".join(c if (c.isalnum() or c in "._-") else "_" for c in str(geom))
    return os.path.join(mesh_path, f"{base_name}_{tag}.shaper")


def _save_current_study(tag: str):
    os.makedirs(mesh_path, exist_ok=True)
    study_path = os.path.abspath(_study_save_path(tag))
    try:
        shutil.rmtree(study_path, ignore_errors=True)
        os.makedirs(study_path, exist_ok=True)
        a_session = ModelAPI_Session.get()
        a_files = StringList()
        a_session.save(study_path, a_files)
        print(f"[study] saved SHAPER session {study_path}")
        return study_path
    except Exception as exc:
        print(f"[study] failed to save {study_path}: {exc}")
        return None


def _finalize_model_for_save():
    # SHAPER content is not persisted into the study until the transaction is
    # published and ended. Saving earlier produces an HDF that opens empty.
    try:
        model.publishToShaperStudy()
    except Exception as exc:
        print(f"[study] publish before save failed: {exc}")
    try:
        model.end()
    except Exception as exc:
        print(f"[study] model.end() before save failed: {exc}")

def _component_has_spline(component):
    if not isinstance(component, dict):
        return False

    pts = component.get("pts", [])
    if isinstance(pts, list):
        for item in pts:
            if (
                isinstance(item, (list, tuple))
                and len(item) > 0
                and isinstance(item[0], str)
                and item[0].lower() == "spline"
            ):
                return True

    sub = component.get("sub_sketches", {})
    if isinstance(sub, dict):
        for child in sub.values():
            if _component_has_spline(child):
                return True
    elif isinstance(sub, list):
        for child in sub:
            if _component_has_spline(child):
                return True

    return False

def _sketch_dict_has_spline(sketch_dict):
    for name, component in sketch_dict.items():
        if name == "shrinkage_factor":
            continue
        if _component_has_spline(component):
            return True
    return False

#def main():
# ----------------- Build sketch dicts (dict-of-dicts) -----------------
build_sketch_dicts = geometry.build_sketch_dicts
build_sig = inspect.signature(build_sketch_dicts)
if "liquid_level" in build_sig.parameters:
    sketch_build = build_sketch_dicts(
        shrinkage_factor,
        liquid_level=config["mesh"]["liquid_level"],
    )
else:
    sketch_build = build_sketch_dicts(shrinkage_factor)
if len(sketch_build) == 4:
    ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping = sketch_build
    sketch_repeat = True
elif len(sketch_build) == 5:
    ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping, sketch_repeat = sketch_build
else:
    raise Exception("build_sketch_dicts must return 4 or 5 values")

has_spline = (
    _sketch_dict_has_spline(ptfe_sketches)
    or _sketch_dict_has_spline(electrode_sketches)
    or _sketch_dict_has_spline(xenon_sketches)
)
use_sketch_repeat = bool(sketch_repeat) and (not has_spline)
print(
    "repeat mode:",
    "SKETCH_REPEAT" if use_sketch_repeat else "FACE_REPEAT",
    "(sketch_repeat=", bool(sketch_repeat), ", has_spline=", has_spline, ")"
)

# ----------------- SALOME Sketching -------------------------
model.begin()
partSet = model.moduleDocument()
Cryostat = model.addPart(partSet)
Cryostat_doc = Cryostat.document()

# --------------------- Draw sketches
ptfe_sketch_objs, ptfe_face_map = sketch_from_dict(ptfe_sketches, Cryostat_doc, makeface = makeface,
                                                apply_shrinkage = False, shrink_below_y = config['mesh']['liquid_level'],
                                                shrinkage_factor = shrinkage_factor,
                                                use_sketch_repeat = use_sketch_repeat)
electrode_sketch_objs, electrode_face_map = sketch_from_dict(electrode_sketches, Cryostat_doc, makeface = makeface,
                                                            use_sketch_repeat = use_sketch_repeat)
xenon_sketch_objs, xenon_face_map = sketch_from_dict(xenon_sketches, Cryostat_doc, makeface = makeface,
                                                    use_sketch_repeat = use_sketch_repeat)
if "GXe_0" not in xenon_face_map.keys(): xenon_face_map["GXe_0"] = None
if "LXe_0" not in xenon_face_map.keys(): xenon_face_map["LXe_0"] = None
model.do()

if stop_after_sketch:
    if save_after_sketch:
        _finalize_model_for_save()
        _save_current_study("after_sketch")
    raise Exception("Stopping after sketching as per debug config") if stop_after_sketch else None
# MAke selection Groups and Partition
try:
    PTFEGrp = None
    ElectrodeGrp = None
    GXeGrp = None
    LXeGrp = None
    if len(ptfe_face_map) > 0:      PTFEGrp = makeGroupByName(Cryostat_doc, "PTFE", ptfe_face_map) 
    if len(electrode_face_map) > 0: ElectrodeGrp = makeGroupByName(Cryostat_doc, "Electrodes", electrode_face_map) 
    if xenon_face_map["GXe_0"]:      GXeGrp = makeGroupByName(Cryostat_doc, "GXe", {"GXe_0":xenon_face_map["GXe_0"]})
    if xenon_face_map["LXe_0"]:      LXeGrp = makeGroupByName(Cryostat_doc, "LXe", {"LXe_0":xenon_face_map["LXe_0"]})

    to_partition = [] 
    if PTFEGrp is not None:       to_partition += [PTFEGrp.result()]
    if ElectrodeGrp is not None:  to_partition += [ElectrodeGrp.result()]
    if GXeGrp is not None:        to_partition += [GXeGrp.result()]
    if LXeGrp is not None:        to_partition += [LXeGrp.result()]

    partition = model.addPartition(Cryostat_doc,to_partition)
    partition.setName("partition_surfaces")
    partition.result().setName("partition_surfaces")
    # Name and Match
    # First by center of weight and area (for all the small components)
    pre_partition_faces = flatten_dict(ptfe_face_map) + flatten_dict(electrode_face_map) + flatten_dict(xenon_face_map)
    pre_partition_face_results = {i.name(): i for i in pre_partition_faces if i is not None}
    ptfe_face_names = set(ptfe_face_map.keys())
    electrode_names = list(electrode_sketches.keys())
    split_bases = {"FieldShapingRings", "FieldShapingGuard"}
    total_partition_subs = partition.result().numberOfSubs()
    pre_partition = {i.name(): [get_area(i, Cryostat_doc), *center_of_weight(Cryostat_doc,i)] for i in pre_partition_faces if i is not None}
    post_partition = {partition.result().subResult(i).name(): [get_area(partition.result().subResult(i), Cryostat_doc), *center_of_weight(Cryostat_doc,partition.result().subResult(i))] for i in range(total_partition_subs)}
    if not manual_only:
        names_in_partition, renamed_cnt = match_and_rename_partition_faces(partition, pre_partition, post_partition)
        print(f"[naming] center/area step named {renamed_cnt} objects ({len(names_in_partition)}/{total_partition_subs} total named)")
    else:
        names_in_partition, renamed_cnt = [], 0
        print(f"[naming] center/area step skipped by manual_only ({len(names_in_partition)}/{total_partition_subs} total named)")
    # Make a residuals list TODO Dont recompute
    post_partition_sub = {partition.result().subResult(i).name():[get_area(partition.result().subResult(i), Cryostat_doc), *center_of_weight(Cryostat_doc,partition.result().subResult(i))] for i in range(total_partition_subs) if ('Partition' in partition.result().subResult(i).name())}
    # Select LXe and GXe by area
    if xenon_face_map["GXe_0"] and xenon_face_map["LXe_0"] and len(post_partition_sub.keys()) > 0 and not manual_only:
        LXeGXeNames, _  = rename_two_largest_partition_faces(partition, post_partition_sub)
        names_in_partition += LXeGXeNames
        print(f"[naming] LXe/GXe area step named {len(LXeGXeNames)} objects ({len(names_in_partition)}/{total_partition_subs} total named)")
    else:
        print(f"[naming] LXe/GXe area step named 0 objects ({len(names_in_partition)}/{total_partition_subs} total named)")

    if not manual_only:
        containment_names, containment_unresolved = rename_partition_faces_by_containment(
            partition,
            pre_partition_face_results,
            ptfe_face_names,
            electrode_names,
            split_bases,
            sample_points=int(config["mesh"].get("partition_containment_sample_points", 2)),
        )
        names_in_partition += containment_names
        print(f"[naming] containment step named {len(containment_names)} objects ({len(names_in_partition)}/{total_partition_subs} total named)")
        for idx, part_name, reason in containment_unresolved:
            print(f"[containment unresolved] idx {idx} ({part_name}): {reason}")
    else:
        print(f"[naming] containment step skipped by manual_only ({len(names_in_partition)}/{total_partition_subs} total named)")

    # Rename by manual mapping
    if len(list(manual_mapping.keys())) > 0:
        if max(list(manual_mapping.keys())) > partition.result().numberOfSubs():
            raise Exception("Highest Manual Naming index larger than number of subs available")
        names = rename_partition_subresults_manual(partition, manual_mapping)
        names_in_partition += names
        print(f"[naming] manual step named {len(names)} objects ({len(names_in_partition)}/{total_partition_subs} total named)")
    else:
        print(f"[naming] manual step named 0 objects ({len(names_in_partition)}/{total_partition_subs} total named)")

    # Improve stability of renaming (Only req for NVidia GPU I think)
    time.sleep(1)
    model.do()
    time.sleep(1)

    # This is included here as selecting later causes recomputation which makes the idx reassignment mess up
    # Selection doesnt work since i name both the result and sub result objects the same thing (I think) so we jsut select the object here FIXME
    ptfe_wall_name_for_charge_buildup = "PTFEWall0_part"
    ptfe_wall_selection = None 
    ## Check everything is named 
    any_unnamed = False
    for i in range(partition.result().numberOfSubs()):
        if "partition" in partition.result().subResult(i).name().lower():
            print("Partition result " + str(i+1) + " is not named (idx:"+str(i)+")" )
            any_unnamed = True
        # FIXME Remove when either naming is fixed or partition naming is fully automated
        if partition.result().subResult(i).name() == ptfe_wall_name_for_charge_buildup: 
            ptfe_wall_selection = partition.result().subResult(i)
    if any_unnamed:
        raise Exception("Stopping: rename partition objects")

    ## Make Name Lists for submesh groups 
    GXe_post_part_names = [i for i in names_in_partition if "GXe" in i]
    LXe_post_part_names = [i for i in names_in_partition if "LXe" in i]
    ptfe_names = list(ptfe_sketches.keys())
    PTFE_post_part_names = []
    for n in names_in_partition:
        for base in ptfe_names:
            # ^Base + optional synthetic face suffix (_<idx>[ _<ridx> ]) + _part + optional auto/manual suffix
            pat = rf"^{re.escape(base)}(?:_\d+(?:_\d+)?)?_part(?:_\d+_(?:manual|auto))?$"
            if re.match(pat, n):
                PTFE_post_part_names.append(n)
                break
    #---------- Make meshing groups
    GXe_faces = [model.selection("FACE", name) for name in GXe_post_part_names]
    GXe_group = model.addGroup(Cryostat_doc, "FACE", GXe_faces)
    GXe_group.setName("GXe_Meshing")
    GXe_group.result().setName("GXe_Meshing")
    # LXE
    LXe_faces = [model.selection("FACE", name) for name in LXe_post_part_names]
    LXe_group = model.addGroup(Cryostat_doc, "FACE", LXe_faces)
    LXe_group.setName("LXe_Meshing")
    LXe_group.result().setName("LXe_Meshing")
    # PTFE
    PTFE_faces = [model.selection("FACE", name) for name in PTFE_post_part_names]
    PTFE_group = model.addGroup(Cryostat_doc, "FACE", PTFE_faces)
    PTFE_group.setName("PTFE_Meshing")
    PTFE_group.result().setName("PTFE_Meshing")
    # Supergroup (must include all meshing groups)
    super_group = model.addGroup(Cryostat_doc, "FACE", GXe_faces  + PTFE_faces + LXe_faces) 
    super_group.setName("Super_Meshing")
    super_group.result().setName("Super_Meshing")
    # -------  Boundary Condition Groups
    Electrode_post_part_by_base = {}
    for n in names_in_partition:
        matched = False
        for base in electrode_names:
            if base in split_bases: 
                # Branch for field shaping that need individual BCs 
                m = re.match(rf"^{re.escape(base)}(\d+)(?:_\d+(?:_\d+)?)?_part(?:_\d+_(?:manual|auto))?$", n)
                if m:
                    nr = m.group(1)
                    key = f"{base}{nr}"
                    Electrode_post_part_by_base.setdefault(key, []).append(n)
                    matched = True
                    break
            elif base != 'shrinkage_factor':
                # Normal rule: optional synthetic face suffix after the base name.
                if re.match(rf"^{re.escape(base)}(?:_\d+(?:_\d+)?)?_part(?:_\d+_(?:manual|auto))?$", n):
                    Electrode_post_part_by_base.setdefault(base, []).append(n)
                    matched = True
                    break
    electrode_grps = []
    for electrode in Electrode_post_part_by_base.keys():
        selec_names = Electrode_post_part_by_base[electrode]
        grp = model.addGroup(Cryostat_doc, "FACE", [model.selection("FACE",i) for i in selec_names])
        grp.setName("BC_"+electrode)
        grp.result().setName("BC_"+electrode)
        electrode_grps.append(grp)

    # ------   For wall charge up we also need a PTFE Wall group
    # We selected earlier swithc to model.selection when the fixme is gone
    # We need to find all Vertex points for each edge
    # There is a bunch of points on the rhs of the wall due to touching the field shaping rings
    if is_tpc:
        tol = 1e-6 # TODO Load config param 
        res, sub = ptfe_wall_selection.resultSubShapePair()
        exp = GeomAPI_ShapeExplorer(res.shape(), GeomAPI_Shape.EDGE)
        i = 0
        lines = {}
        while exp.more():
            cur = exp.current()
            ed = cur.edge()
            if ed is None:
                print(f"[{i}] cur.edge() is None -> skip")
                exp.next()
                i += 1
                continue

            p1 = ed.firstPoint()
            p2 = ed.lastPoint()
            x1, y1, z1 = p1.x(), p1.y(), p1.z()
            x2, y2, z2 = p2.x(), p2.y(), p2.z()

            dx = x2 - x1
            dy = y2 - y1

            if abs(dx) <= tol and abs(dy) > tol:
                key = ("V", round(0.5*(x1+x2), 12))
            elif abs(dy) <= tol and abs(dx) > tol:
                key = ("H", round(0.5*(y1+y2), 12))
            else:
                key = ("O", i)  # other / degenerate / angled

            lines.setdefault(key, []).append({
                "idx": i,
                "edge_shape": cur,      # GeomAPI_Shape for the edge (handle)
                "p1": (x1, y1, z1),
                "p2": (x2, y2, z2),
            })

            exp.next()
            i += 1

        # --- pick the left-most vertical column ---
        v_keys = [k for k in lines.keys() if k[0] == "V"]
        if not v_keys:
            raise Exception("No vertical lines found")
        left_key = min(v_keys, key=lambda k: k[1])   # smallest x
        left_edges = lines[left_key]
        # --- convert edge shapes -> selections (try both common bindings) ---
        edge_sels = []
        for item in left_edges:
            edge_shape = item["edge_shape"]   # GeomAPI_Shape (EDGE)
            edge_sels.append(ModelHighAPI_Selection(res, edge_shape))

        PTFE_Wall_group = model.addGroup(Cryostat_doc, "EDGE", edge_sels )
        PTFE_Wall_group.setName("PTFE_Wall_group")
        PTFE_Wall_group.result().setName("PTFE_Wall_Charge")

    time.sleep(1)
    model.do()
    time.sleep(1)

    # Publish and collect groups
    model.publishToShaperStudy()
    all_groups = SHAPERSTUDY.shape(model.featureStringId(partition))
    # Build a name to group dictionary
    named_groups = {g.GetName(): g for g in all_groups}
    partition_surfaces      = named_groups.get("partition_surfaces")
    GXeMesh   = named_groups.get("GXe_Meshing")
    LXeMesh   = named_groups.get("LXe_Meshing")
    PTFEMesh  = named_groups.get("PTFE_Meshing")
    SuperMesh = named_groups.get("Super_Meshing")
    #For surface charge
    used_names = (i.GetName() for i in (partition_surfaces, GXeMesh, LXeMesh, PTFEMesh, SuperMesh) if i is not None)
    BC_Groups = [g for g in all_groups if g.GetName() not in used_names]
    model.end()
except Exception:
    if save_on_partition_failure:
        _finalize_model_for_save()
        _save_current_study("partition_failure")
    raise

# Mesh : TODO Make yaml entry for this
smesh = smeshBuilder.New()
mesh = smesh.Mesh(SuperMesh)
NETGEN_1D_2D = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
# ---------- Configure Netgen Meshing rules -----------

params = config['mesh']['NetgenParams']
# https://docs.salome-platform.org/latest/gui/NETGENPLUGIN/netgen_2d_3d_hypo_page.html
NETGEN_2D_Params = NETGEN_1D_2D.Parameters()
# Get a preset

if params["debugCoarse"]:
    NETGEN_2D_Params.SetFineness(smeshBuilder.VeryCoarse)
else:
    NETGEN_2D_Params.SetFineness(smeshBuilder.Custom)  
    # --- General sizing ---
    # Min and max edge lengths
    NETGEN_2D_Params.SetMaxSize(params['maxSize']) # 1/100 of our largest part
    NETGEN_2D_Params.SetMinSize(params['minSize'])        # Dont set a minimum
    NETGEN_2D_Params.SetGrowthRate(params['growthRate'])  # Mesh element growth rate relative to one another 
    # --- Curvature-driven sizing ---
    # Lots of curves we want to optimize against
    NETGEN_2D_Params.SetUseSurfaceCurvature(params['UseSurfaceCurvature'])
    # Minimum segments per topological edge
    NETGEN_2D_Params.SetNbSegPerEdge(params['NBSegmentsPerEdge'])
    # For a circle the local target size is radius / N
    # Number of segments on a circle then 2piR / (R/N) = 2piN
    # So with N NbSegPerRadius we get ~2piN segments around the circle
    # So for 12 edges on a circle we want ~2=N
    NETGEN_2D_Params.SetNbSegPerRadius(params['NBSegmentsPerEdge']) # 8  
    # Alternatively we can use this with a maximum allowed deviation - not enabled
    NETGEN_2D_Params.SetChordalErrorEnabled(params['UseChordalError'])
    NETGEN_2D_Params.SetChordalError(params['ChordalErrorValue'])
    # --- Local Size Control ---
    # This can be used to set local sizes on individual shapes
    # A source file can also be used
    #NETGEN_2D_Params.SetLocalSizeOnShape
    # --- Advanced Quality Algorithms ---
    NETGEN_2D_Params.SetElemSizeWeight(params['AQ_ElemSizeWeight']) # Prioritize mesh quality over accuracy
    NETGEN_2D_Params.SetUseDelauney(params['AQ_UseDelauney']) # advancing front is a bit more fragile but sometimes nicer mesh
    NETGEN_2D_Params.SetCheckOverlapping(params['AQ_CheckOverlapping']) # Surface elements cant overlap
    NETGEN_2D_Params.SetCheckChartBoundary(params['AQ_CheckChartBoundary']) # Strict checking
    NETGEN_2D_Params.SetFuseEdges(params['AQ_FuseEdges'])        # Avoid Duplicate edges 
    # --- Element Types and Optimization ---
    NETGEN_2D_Params.SetQuadAllowed(params['Opt_QuadAllowed'])  # Dont allow quad meshes
    NETGEN_2D_Params.SetSecondOrder(params['Opt_SecondOrder'])      # Overkill for our purposes - also not sure MFEM supports this by default
    NETGEN_2D_Params.SetOptimize(params['Opt_Optimize'])      # Post meshing optimization of the mesh 

# ================================ Surface Groups ====================================
MSH_GXe_face_group  = mesh.GroupOnGeom(GXeMesh,  "GXe",  SMESH.FACE)
MSH_LXe_face_group  = mesh.GroupOnGeom(LXeMesh,  "LXe",  SMESH.FACE)
MSH_PTFE_face_group = mesh.GroupOnGeom(PTFEMesh, "PTFE", SMESH.FACE)
# ================================ Boundary Groups ====================================
BCs = [mesh.GroupOnGeom(msh_grp, msh_grp.GetName(), SMESH.EDGE) for msh_grp in BC_Groups]
# ------------------  Need to compute such that Borders exist ------------------
mesh.Compute()
# --------------------------------  Also mark the cryostat boundary --------------------------------
# Select only the free borders of the mesh (ie edges that dont connect to a triangle) 
if mark_external_boundary:
    free_border_edges = mesh.MakeGroup("FreeBordersAll",
                                        SMESH.EDGE,
                                        SMESH.FT_FreeBorders)

    all_bc_edges = None
    if len(BCs) > 0:
        # Subtract all already-defined BC edges from free borders.
        all_bc_edges = mesh.UnionListOfGroups(BCs, "AllBC_Edges")
        all_bc_edges = mesh.ConvertToStandalone(all_bc_edges)
        free_borders_clean = mesh.CutGroups(
            free_border_edges, 
            all_bc_edges,
            "FreeBorders_NoBC"
        )
        free_borders_clean = mesh.ConvertToStandalone(free_borders_clean)
    else:
        free_borders_clean = mesh.ConvertToStandalone(free_border_edges)

    edge_ids = free_borders_clean.GetIDs()
    keep_edges = list(edge_ids)
    r0_edges = []

    if split_axis_boundary:
        # Split out edges lying on x=0 axis to a separate boundary group.
        def on_axis(node_id):
            x, y, z = mesh.GetNodeXYZ(node_id)
            return abs(x) < axis_boundary_tol

        keep_edges = []
        r0_edges = []
        for eid in edge_ids:
            nodes = mesh.GetElemNodes(eid)
            if not (on_axis(nodes[0]) and on_axis(nodes[1])):
                keep_edges.append(eid)
            else:
                r0_edges.append(eid)

    cleaned = mesh.CreateEmptyGroup(SMESH.EDGE, external_boundary_name)
    if len(keep_edges) > 0:
        cleaned.Add(keep_edges)
    cleaned = mesh.ConvertToStandalone(cleaned)

    if split_axis_boundary:
        cleaned2 = mesh.CreateEmptyGroup(SMESH.EDGE, axis_boundary_name)
        if len(r0_edges) > 0:
            cleaned2.Add(r0_edges)
        cleaned2 = mesh.ConvertToStandalone(cleaned2)

    # Delete helper groups
    mesh.RemoveGroup(free_borders_clean)
    if all_bc_edges is not None:
        mesh.RemoveGroup(all_bc_edges)
    mesh.RemoveGroup(free_border_edges)


# --------------------------------- And save the mesh --------------------------------- 
if not os.path.isdir(mesh_path):
    os.mkdir(mesh_path)
mesh.ExportMED(
    os.path.join(mesh_path, "mesh.med"),
    auto_groups=False,  
    minor=42,            # MED 4.2
)
#return 

#if __name__ == '__main__':
#    main()
