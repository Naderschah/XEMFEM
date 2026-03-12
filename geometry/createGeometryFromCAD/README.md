
## Autogenerating the XENONnT simulation geometry from the CAD file


### 2D Slices

The file I received is STEP AP214, this contains many linked objects producing a whole host of issues when trying to resolve global placement of topological placements in freeCAD directly. 

To check which you have inspect the file header, it contains a line: `FILE_DESCRIPTION (( 'STEP AP214' ),`. 

So instead, each slice is saved as a DXF and then this is postprocessed to produce the simulation geometry.

The order of operations is 

```
pip install pyyaml ezdxf
# Produce all DXF component level slices, with audit file in case of problems
python3 ExportSliceToDXF.py \
  --workers 4 \
  --component-audit-csv DXF_slices_parts/full_audit.csv
# DXF to intermediary JSON format, SALOME capable but takes very long 
python3 createSerialization.py
# Finally group and cleanup with the cleanup config
python3 cleanup_jsons.py cleanup_config.yaml
```

The cylinder height is in y, the 0 degree slice corrsponds to the XY plane in this code. Unfortunately the angles are not in the xenon coordinate system, so I don't actually know which slice applies to data. So I picked an arbitrary PMT that forms a line with other PMTs and the central PMT, the picked PMT is at a distance of 607.75 mm in x and 80.01 mm in y from the central PMT, producing a rotation of 7.4998 degrees to align a PMT row with the XY plane.  

Note that DXF BSplines still have a weight field in my export, this is usually interpreted as a NURBS curve and then fails as there are no weights. For the DXF to Salome translator this is not an issue.

### 3D Geometry 

We require a AP242 STEP file the AP214 will not open in reasonable time in Salome, the format can be created by opening and exporting in FreeCAD. Then open it in salome (this took on my device around 5-6 hours) and then immediately save it, format is hdf. 