import os
import sys
import yaml
import time 
from pathlib import Path
import importlib.util

# ============== Paths for path unaware TUI =====================
base_path = '/work/geometry/'
sys.path.insert(0, base_path)
from helper_functions import *

# ==============  Load Config File ===========================
with open(base_path + "config.yaml") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# Resolve which geometry file to use
fname = [base_path + i for i in os.listdir(base_path) if (config['mesh']['geometry'] in i) and (".py" in i)]
if len(fname) == 1:
    module_name = Path(fname[0]).stem # Strip .py
    spec = importlib.util.spec_from_file_location(module_name, fname[0])
    geometry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
else:
    raise Exception("Found " + len(fname) + " choices for geometry source " + config['mesh']['geometry'])

# ==============  Debug Options ===========================
makeface = True
stop_after_sketch = False

#def main():
shrinkage_factor = 1.0 - 0.014

# ----------------- Build sketch dicts (dict-of-dicts) -----------------
ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping = geometry.build_sketch_dicts(shrinkage_factor)

# ----------------- SALOME Sketching -------------------------
model.begin()
partSet = model.moduleDocument()
Cryostat = model.addPart(partSet)
Cryostat_doc = Cryostat.document()

# --------------------- Draw sketches
ptfe_sketch_objs, ptfe_face_map = sketch_from_dict(ptfe_sketches, Cryostat_doc, makeface = makeface,
                                                apply_shrinkage = False, shrink_below_y = geometry.liquid_level,
                                                shrinkage_factor = shrinkage_factor)
electrode_sketch_objs, electrode_face_map = sketch_from_dict(electrode_sketches, Cryostat_doc, makeface = makeface)
xenon_sketch_objs, xenon_face_map = sketch_from_dict(xenon_sketches, Cryostat_doc, makeface = makeface)

model.do()

#if stop_after_sketch: return

# MAke selection Groups and Partition
PTFEGrp = makeGroupByName(Cryostat_doc, "PTFE", ptfe_face_map) 
ElectrodeGrp = makeGroupByName(Cryostat_doc, "Electrodes", electrode_face_map) 
GXeGrp = makeGroupByName(Cryostat_doc, "GXe", {"GXe0":xenon_face_map["GXe0"]}) 
LXeGrp = makeGroupByName(Cryostat_doc, "LXe", {"LXe0":xenon_face_map["LXe0"]}) 

partition = model.addPartition(Cryostat_doc,[GXeGrp.result(), LXeGrp.result(), PTFEGrp.result(), ElectrodeGrp.result()])
partition.setName("partition_surfaces")
partition.result().setName("partition_surfaces")

# Name and Match
# First by center of weight and area (for all the small components)
pre_partition_faces = flatten_dict(ptfe_face_map) + flatten_dict(electrode_face_map) + flatten_dict(xenon_face_map)
pre_partition = {i.name(): [get_area(i, Cryostat_doc), *center_of_weight(Cryostat_doc,i)] for i in pre_partition_faces}
post_partition = {partition.result().subResult(i).name(): [get_area(partition.result().subResult(i), Cryostat_doc), *center_of_weight(Cryostat_doc,partition.result().subResult(i))] for i in range(partition.result().numberOfSubs())}
names_in_partition, renamed_cnt = match_and_rename_partition_faces(partition, pre_partition, post_partition)

# Select LXe and GXe by area
LXeGXeNames, _  = rename_two_largest_partition_faces(partition, post_partition)
names_in_partition += LXeGXeNames
print("Named LXe and GXe by Volume")

"""
Selection by containment prooved very difficult, 
there is a function ported to python via GeomAlgoAPI_ShapeTools
https://docs.salome-platform.org/latest/tui/SHAPER/classGeomAlgoAPI__ShapeTools.html#a3f2a4dbb9d3479ef975593de9db7a329
However I could not figure out how to use it reliably
Points on edges are not guaranteed to return true for containment so a few points inside the 
shape must be sampled, for convenience 
https://docs.salome-platform.org/latest/tui/SHAPER/classGeomAlgoAPI__PointCloudOnFace.html#details
may be used, but I didn't try it yet.

TODO Implement this
"""

# Rename by manual mapping
if max(list(manual_mapping.keys())) > partition.result().numberOfSubs():
    raise Exception("Highest Manual Naming index larger than number of subs available")
names = rename_partition_subresults_manual(partition, manual_mapping)
names_in_partition += names

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
        # ^Base + optional digits + _part + optional _<nr>_manual
        pat = rf"^{re.escape(base)}(\d+)?_part(_\d+_manual)?$"
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
electrode_names = list(electrode_sketches.keys())
split_bases = {"FieldShapingRings", "FieldShapingGuard"}
Electrode_post_part_by_base = {}
for n in names_in_partition:
    matched = False
    for base in electrode_names:
        if base in split_bases: 
            # Branch for field shaping that need individual BCs 
            m = re.match(rf"^{re.escape(base)}(\d+)_part(_\d+_manual)?$", n)
            if m:
                nr = m.group(1)
                key = f"{base}{nr}"
                Electrode_post_part_by_base.setdefault(key, []).append(n)
                matched = True
                break
        elif base != 'shrinkage_factor':
            # Normal rule: optional number, key by base (e.g. Gate, Anode, Cathode, etc.)
            if re.match(rf"^{re.escape(base)}(\d+)?_part(_\d+_manual)?$", n):
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
used_names = (i.GetName() for i in (partition_surfaces, GXeMesh, LXeMesh, PTFEMesh, SuperMesh))
BC_Groups = [g for g in all_groups if g.GetName() not in used_names]
model.end()

# Mesh : TODO Make yaml entry for this
smesh = smeshBuilder.New()
mesh = smesh.Mesh(SuperMesh)
NETGEN_1D_2D = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
# ---------- Configure Netgen Meshing rules -----------

params = config['mesh']['NetgenParams']
# https://docs.salome-platform.org/latest/gui/NETGENPLUGIN/netgen_2d_3d_hypo_page.html
NETGEN_2D_Params = NETGEN_1D_2D.Parameters()
# Get a preset
preset = None
if params["preset"] == "VeryFine": preset = smeshBuilder.VeryFine
elif params["preset"] == "Fine": preset = smeshBuilder.Fine
elif params["preset"] == "Moderate": preset = smeshBuilder.Moderate
elif params["preset"] == "Coarse": preset = smeshBuilder.Coarse
elif params["preset"] == "VeryCoarse": preset = smeshBuilder.VeryCoarse
else: raise Exception("Netgen Preset " + params["preset"] + " not recognized")
NETGEN_2D_Params.SetFineness(preset)  
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
free_border_edges = mesh.MakeGroup("FreeBordersAll",
                                    SMESH.EDGE,
                                    SMESH.FT_FreeBorders)
# Subtract every single boundary condition we have
all_bc_edges = mesh.UnionListOfGroups(BCs, "AllBC_Edges") # We need a union group here
all_bc_edges = mesh.ConvertToStandalone(all_bc_edges)
free_borders_clean = mesh.CutGroups(
    free_border_edges, 
    all_bc_edges,
    "FreeBorders_NoBC"
)   
free_borders_clean = mesh.ConvertToStandalone(free_borders_clean)
# Pop out the x = 0 zero lines 
tol = 1e-9
def on_axis(node_id):
    x, y, z = mesh.GetNodeXYZ(node_id)
    return abs(x) < tol
edge_ids = free_borders_clean.GetIDs()
keep_edges = []
r0_edges = []
for eid in edge_ids:
    nodes = mesh.GetElemNodes(eid)
    if not (on_axis(nodes[0]) and on_axis(nodes[1])):
        keep_edges.append(eid)
    else:
        r0_edges.append(eid)

cleaned = mesh.CreateEmptyGroup(SMESH.EDGE, "BC_Cryostat")
cleaned.Add(keep_edges)
cleaned = mesh.ConvertToStandalone(cleaned)

cleaned2 = mesh.CreateEmptyGroup(SMESH.EDGE, "BC_r0")
cleaned2.Add(keep_edges)
cleaned2 = mesh.ConvertToStandalone(cleaned2)
# Delete helpers 
mesh.RemoveGroup(free_borders_clean)
mesh.RemoveGroup(all_bc_edges)
mesh.RemoveGroup(free_border_edges)


# --------------------------------- And save the mesh --------------------------------- 
if not os.path.isdir('mesh/'):
    os.mkdir("mesh")
mesh.ExportMED(
    base_path + "mesh/mesh.med",
    auto_groups=False,  
    minor=42,            # MED 4.2
)
#    return 

#if __name__ == '__main__':
#    main()
