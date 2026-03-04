def build_sketch_dicts(shrinkage_factor):
  wedge_base = 1.0
  wedge_height = 1.0

  y_offset = 0.05

  gap = 0.1
  groundedPlateH = 0.1

  # We define as if a tpc to reuse the existing infrastructure
  electrode_sketches = {}
  xenon_sketches = {}
  ptfe_sketches = {}

  electrode_sketches["Wedge"] = {
    'pts': [
      ['line', 0.0,0],
      ['line', wedge_base/2, wedge_height],
      ['line', -wedge_base/2,wedge_height],
    ]
  }
  
  # Dielectric - Infinite Space
  xenon_sketches["LXe"] = {
    'RadialPosition': 0,
    'VerticalPosition': 0,
    'Radius': 10*wedge_height,
  }

  # Simulation Volume
  xenon_sketches["GXe"] = {
    'RadialPosition': 0,
    'VerticalPosition': 0,
    'Radius': 4*wedge_height,
  }
  
  manual_mapping = {
    1-1: "Wedge",
    2-1: "GXe",
    3-1: "LXe"
  }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
