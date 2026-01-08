import os
import sys
import yaml

# ============== Paths for path unaware TUI =====================
# FIXME: Find a better way of doing this
base_path = '/home/felix/MFEMElectrostatics/geometry/'
sys.path.insert(0, base_path)
from helper_functions import *

# ==============  Load Config File ===========================
with open(base_path + "config.yaml") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

if config['mesh']['geometry'] == "SR3":
    import SR3_constants as geometry
else:
    raise Exception("Geometry for " + config['mesh']['geometry'] + " not Implemented")
# ==============  Debug Options ===========================
makeface = True
stop_after_sketch = False

def main():
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

    if stop_after_sketch: return

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

    ## Check everything is named 
    any_unnamed = False
    for i in range(partition.result().numberOfSubs()):
        if "partition" in partition.result().subResult(i).name().lower():
            print("Partition result " + str(i+1) + " is not named (idx:"+str(i)+")" )
            any_unnamed = True
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
    # Make meshing groups
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
    # Boundary Condition Groups
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
    import time 
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
    used_names = (i.GetName() for i in (partition_surfaces, GXeMesh, LXeMesh, PTFEMesh, SuperMesh))
    BC_Groups = [g for g in all_groups if g.GetName() not in used_names]
    model.end()

    # Mesh : TODO Make yaml entry for this
    smesh = smeshBuilder.New()
    mesh = smesh.Mesh(SuperMesh)
    NETGEN_1D_2D = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
    # ---------- Configure Netgen Meshing rules -----------
    # https://docs.salome-platform.org/latest/gui/NETGENPLUGIN/netgen_2d_3d_hypo_page.html
    NETGEN_2D_Params = NETGEN_1D_2D.Parameters()
    # Get a preset
    NETGEN_2D_Params.SetFineness(smeshBuilder.VeryFine)          # VeryCoarse, Coarse, Moderate, Fine, VeryFine, Custom
    # --- General sizing ---
    # Min and max edge lengths
    NETGEN_2D_Params.SetMaxSize(ptfe_sketches["PTFEWall"]["Height"] / 100) # 1/100 of our largest part
    NETGEN_2D_Params.SetMinSize(0)        # Dont set a minimum
    NETGEN_2D_Params.SetGrowthRate(0.1)  # Mesh element growth rate relative to one another 
    # --- Curvature-driven sizing ---
    # Lots of curves we want to optimize against
    NETGEN_2D_Params.SetUseSurfaceCurvature(1)
    # Minimum segments per topological edge
    NETGEN_2D_Params.SetNbSegPerEdge(2)
    # For a circle the local target size is radius / N
    # Number of segments on a circle then 2piR / (R/N) = 2piN
    # So with N NbSegPerRadius we get ~2piN segments around the circle
    # So for 12 edges on a circle we want ~2=N
    NETGEN_2D_Params.SetNbSegPerRadius(1) # 8  
    # Alternatively we can use this with a maximum allowed deviation - not enabled
    NETGEN_2D_Params.SetChordalErrorEnabled(0)
    NETGEN_2D_Params.SetChordalError(electrode_sketches["AnodeElectrode"]["sub_sketches"]["AnodeWires"]["Radius"] * 0.1)
    # --- Local Size Control ---
    # This can be used to set local sizes on individual shapes
    # A source file can also be used
    #NETGEN_2D_Params.SetLocalSizeOnShape
    # --- Advanced Quality Algorithms ---
    NETGEN_2D_Params.SetElemSizeWeight(0.5) # Prioritize mesh quality over accuracy
    NETGEN_2D_Params.SetUseDelauney(1) # advancing front is a bit more fragile but sometimes nicer mesh
    NETGEN_2D_Params.SetCheckOverlapping(1) # Surface elements cant overlap
    NETGEN_2D_Params.SetCheckChartBoundary(2) # Strict checking
    NETGEN_2D_Params.SetFuseEdges(1)        # Avoid Duplicate edges 
    # --- Element Types and Optimization ---
    NETGEN_2D_Params.SetQuadAllowed(False)  # Dont allow quad meshes
    NETGEN_2D_Params.SetSecondOrder(0)      # Overkill for our purposes - also not sure MFEM supports this by default
    NETGEN_2D_Params.SetOptimize(True)      # Post meshing optimization of the mesh 

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
    for eid in edge_ids:
        nodes = mesh.GetElemNodes(eid)
        if not (on_axis(nodes[0]) and on_axis(nodes[1])):
            keep_edges.append(eid)
    cleaned = mesh.CreateEmptyGroup(SMESH.EDGE, "BC_Cryostat")
    cleaned.Add(keep_edges)
    cleaned = mesh.ConvertToStandalone(cleaned)
    # Delete helpers 
    mesh.RemoveGroup(free_borders_clean)
    mesh.RemoveGroup(all_bc_edges)
    mesh.RemoveGroup(free_border_edges)
    # And save the mesh 
    if not os.path.isdir('mesh/'):
        os.mkdir("mesh")
    mesh.ExportMED(
        base_path + "mesh/mesh.med",
        auto_groups=False,  
        minor=42,            # MED 4.2
    )
    return 

if __name__ == '__main__':
    main()
