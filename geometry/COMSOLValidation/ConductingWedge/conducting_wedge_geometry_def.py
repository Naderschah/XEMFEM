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
      ['line', 0.0,-wedge_height+y_offset],
      ['line', wedge_base/2, y_offset],
      ['line', -wedge_base/2, y_offset],
    ]
  }

  electrode_sketches["GroundedPlate"] = {
    'RadialPosition': -wedge_base,
    'VerticalPosition': -wedge_height-gap-groundedPlateH,
    'Width': 2*wedge_base,
    'Height': groundedPlateH,
  }
  
  # Dielectric
  xenon_sketches["LXe"] = {
    'RadialPosition': -wedge_base/2,
    'VerticalPosition': -wedge_height-gap,
    'Width': wedge_base,
    'Height': gap,
  }

  # Simulation Volume
  xenon_sketches["GXe"] = {
    'RadialPosition': -wedge_base,
    'VerticalPosition': -wedge_height-gap-groundedPlateH,
    'Width': 2*wedge_base,
    'Height': +wedge_height+gap+groundedPlateH + 0.1,
  }
  
  manual_mapping = {
    1-1: "GroundedPlate",
    2-1: "Wedge",
    3-1: "GXe",
    4-1: "LXe",
  }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
