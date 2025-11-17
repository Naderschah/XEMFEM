# ------------------------------------------------------------------
# PRECONDITIONS:
#   - LXeRegion already exists (cut LXe domain)
#   - imported_shapes["Cathode"] exists (electrode)
# ------------------------------------------------------------------

region_label = "LXeRegion"
region = LXeRegion
ename  = "Cathode"

E = imported_shapes[ename]

print("Region:", region_label, "->", region)
print("Electrode:", ename, "->", E)

# ------------------------------------------------------------------
# 1) Geometric intersection: Section(region, electrode)
# ------------------------------------------------------------------
sec = geompy.MakeSection(region, E)
geompy.addToStudy(sec, f"Section_{region_label}_{ename}")

print("WhatIs Section:")
print(geompy.WhatIs(sec))

sec_edges = geompy.ExtractShapes(sec, geompy.ShapeType["EDGE"], True)
print(f"Number of section edges: {len(sec_edges)}")

for i, e in enumerate(sec_edges, start=1):
    geompy.addToStudy(e, f"SecEdge_{region_label}_{ename}_{i}")

# ------------------------------------------------------------------
# 2) Map each section edge back onto the region using GetSame
# ------------------------------------------------------------------
collected = []

for i, e in enumerate(sec_edges, start=1):
    same_compound = geompy.GetSame(region, e)   # GEOM_Object or None
    if same_compound is None:
        print(f"Section edge {i}: no matching region edges")
        continue

    # Extract EDGE sub-shapes from the returned compound
    same_edges = geompy.ExtractShapes(same_compound, geompy.ShapeType["EDGE"], True)
    print(f"Section edge {i}: found {len(same_edges)} matching region edges")

    for j, ee in enumerate(same_edges, start=1):
        geompy.addToStudy(ee, f"BCedge_{region_label}_{ename}_{i}_{j}")
    collected.extend(same_edges)

print(f"Total matching region edges collected: {len(collected)}")

if collected:
    if len(collected) == 1:
        bc_shape = collected[0]
    else:
        bc_shape = geompy.MakeCompound(collected)

    bc_name = f"BC_{ename}_{region_label}"
    geompy.addToStudy(bc_shape, bc_name)
    print(f"Created BC geometry: {bc_name}")
else:
    print(f"WARNING: no matching region edges for electrode '{ename}' on {region_label}")
