import os, json, glob

def build_sketch_dicts(_unused=None, angle_dir = "/work/geometry/createGeometryFromCAD/DXF_slices_parts/slice_015.00deg"):
    """
    Debug-only helper.

    Returns:
      ptfe_sketches       : dict
      electrode_sketches  : dict
      xenon_sketches      : dict
      manual_mapping      : dict (always empty)

    Each sketches dict has the form:
      { comp_key : { "pts": [...] }, ... }

    Components are arbitrarily distributed to guarantee
    each dict has at least one entry (if possible).
    """

    merged = {}

    # Load all JSON files in this angle directory
    for json_path in sorted(glob.glob(os.path.join(angle_dir, "*.json"))):
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            if not isinstance(data, dict):
                continue
            for k, v in data.items():
                merged[k] = v

    if not merged:
        raise RuntimeError(f"No component JSON files found in {angle_dir}")

    # Split arbitrarily but deterministically
    keys = list(merged.keys())

    ptfe_sketches = {}
    electrode_sketches = {}
    xenon_sketches = {}

    for i, k in enumerate(keys):
        if i % 3 == 0:
            ptfe_sketches[k] = merged[k]
        elif i % 3 == 1:
            electrode_sketches[k] = merged[k]
        else:
            xenon_sketches[k] = merged[k]

    # Ensure non-empty buckets if possible
    buckets = [ptfe_sketches, electrode_sketches, xenon_sketches]
    non_empty = [b for b in buckets if b]

    if len(non_empty) == 1:
        # all data landed in one bucket → copy one element into others
        k = next(iter(non_empty[0]))
        for b in buckets:
            if not b:
                b[k] = merged[k]

    manual_mapping = {}

    return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
  