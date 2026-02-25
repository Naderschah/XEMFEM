import FreeCAD as App
import importDXF
import os

path = r"/home/felix/MFEMElectrostatics/geometry/createGeometryFromCAD/DXF_slices_parts/slice_022.50deg/"

doc = App.newDocument("DXF_Project")

for f in os.listdir(path):
    if f.lower().endswith(".dxf"):
        if "XENT_TPC_A_Warm.Part__Feature088.Part__Feature088" in f:
            importDXF.insert(os.path.join(path, f), doc.Name)

doc.recompute()