def build_sketch_dicts(shrinkage_factor):
  inner_cylinder_R = 0.2
  inner_offset_x = 0.1
  outer_cylinder_R = 1.0

  cylinder_thicknes = 0.2

  # Outer Simulation Volume
  x0 = - 1.5
  y0 = - 1.5
  W = 3
  H = 3

  # We define as if a tpc to reuse the existing infrastructure
  electrode_sketches = {}
  xenon_sketches = {}
  ptfe_sketches = {}

  electrode_sketches["InnerCylinder"] = {
    'RadialPosition': inner_offset_x,
    'VerticalPosition': 0.0,
    'Radius': inner_cylinder_R,
  }
  electrode_sketches["OuterCylinder"] = {
    'RadialPosition': 0.0,
    'VerticalPosition': 0.0,
    'Radius': outer_cylinder_R + cylinder_thicknes,
  }
  # Dielectric
  xenon_sketches["GXe"] = {
    'RadialPosition': 0.0,
    'VerticalPosition': 0.0,
    'Radius': outer_cylinder_R,
  }

  # Simulation Volume
  xenon_sketches["LXe"] = {
    'RadialPosition': x0,
    'VerticalPosition': y0,
    'Width': W,
    'Height': H,
  }
  
  manual_mapping = {
    1-1: "OuterCylinder",
    2-1: 'LXe_comp',
    3-1: "InnerCylinder",
    4-1: "GXe_comp"
  }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
