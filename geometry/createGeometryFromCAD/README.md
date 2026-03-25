
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

`cleanup_jsons.py` now supports explicit pre-group instructions in `cleanup_config.yaml`:

- `include`: merge extra YAML files into the current config
- `ignore_components`: regex rules for components that should be dropped before grouping
- `replacements`: replace one component by another before grouping
  - `target`: one JSON stem to replace
  - `targets`: several JSON stems that together form one faulty component
  - `output_name`: optional stem for the repaired component; defaults to the first target
  - `source`: JSON stem in the same slice, or `{path: ...}` for an external JSON
  - `periodicity`: optional `{axis, pitch, number}`
  - `slices`: optional regex list; omitted means all slices
- `materials`: per-slice material classification written to `cleanup_autogen.yaml`
  - currently used for `GXe`, `LXe`, and default `PTFE`
- `electrodes`: per-slice electrode hull definitions written to `cleanup_autogen.yaml`
  - each electrode becomes one hull with listed subcomponents
  - either define `components` directly for all slices, or use nested slice keys such as
    `slice_022.50:` / `slice_037.50:` for slice-specific component lists
  - nested `slice_*` keys are treated as exact slice shorthands; omitting the trailing `deg` is allowed
- `shape_match_tol`: strict tolerance for deciding whether two components are the same shape
- `periodicity_tol`: looser tolerance for detecting repeated spacing/pitch
- `grouping_mode`: `periodic` or `none`
  - `none` disables all grouping and override-periodicity handling; cleanup only applies ignore/replacement/manual-addition steps and writes the resulting JSONs plus slice Python
- `paths.in_root` / `paths.out_root` / `paths.slice_glob`: batch mode to process many `slice_*` directories in one run
  - `in_dir` / `out_dir` still work for a single slice

Execution order is:

1. Load and merge config includes
2. Apply slice-filtered `ignore_components`
3. Apply slice-filtered `replacements` into the working component set
4. Load slice-filtered `manual_additions` into the working set when they target the cleaned output directory
5. Group repeating patterns on the repaired component set
6. Write `cleanup_autogen.yaml` and `cleanup_autogen_constants.py` for debug/autogen loading

The cylinder height is in y, the 0 degree slice corrsponds to the XY plane in this code. Unfortunately the angles are not in the xenon coordinate system, so I don't actually know which slice applies to data. So I picked an arbitrary PMT that forms a line with other PMTs and the central PMT, the picked PMT is at a distance of 607.75 mm in x and 80.01 mm in y from the central PMT, producing a rotation of 7.4998 degrees to align a PMT row with the XY plane.  

Note that DXF BSplines still have a weight field in my export, this is usually interpreted as a NURBS curve and then fails as there are no weights. For the DXF to Salome translator this is not an issue.

### 3D Geometry 

We require a AP242 STEP file the AP214 will not open in reasonable time in Salome, the format can be created by opening and exporting in FreeCAD. Then open it in salome (this took on my device around 5-6 hours) and then immediately save it, format is hdf. 
