import FreeCAD as App
import importDXF
import os

path = r"/absolute/path/to/dxf_folder"

doc = App.newDocument("DXF_Project")

for f in os.listdir(path):
    if f.lower().endswith(".dxf"):
        importDXF.insert(os.path.join(path, f), doc.Name)

doc.recompute()