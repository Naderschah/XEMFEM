
## Autogenerating the XENONnT simulation geometry from the CAD file


### 2D Slices

The file I received is STEP AP214, this contains many linked objects producing a whole host of issues when trying to resolve global placement of topological placements in freeCAD directly. 

To check which you have inspect the file header, it contains a line: `FILE_DESCRIPTION (( 'STEP AP214' ),`. 

So instead, each slice is saved as a DXF and then this is postprocessed to produce the simulation geometry.

The order of operations is 

```
# Produce all slices
python3 ExportSliceToDXF.py
# Distribute DXF in file structure rather than in DXF, also produces a hash lookup to allow for copy operations (faster to construct in salome)
python3 distribute_components_and_check_same.py  # Uses the naming_convention.csv file and checks if they are marked for meshing 
# We switch to the intermediary format saving the jsons in the same location, this will err on unclosed contours 
dxf_to_intermediary
```


### 3D Geometry 

We require a AP242 STEP file the AP214 will not open in reasonable time in Salome, the format can be created by opening and exporting in FreeCAD. Then open it in salome (this took on my device around 5-6 hours) and then immediately save it, format is hdf. 