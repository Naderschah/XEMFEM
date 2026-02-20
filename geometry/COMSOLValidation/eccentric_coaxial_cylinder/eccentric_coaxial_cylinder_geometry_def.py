def build_sketch_dicts(shrinkage_factor):
  inner_cylinder_R = 0.5
  inner_offset_x = 0.2
  outer_cylinder_R = 1.0

  cylinder_thicknes = 0.1

  # Outer Simulation Volume
  x0 = - 1.5
  y0 = - 1.5
  W = 3
  H = 3


  # We define as if a tpc to reuse the existing infrastructure
  electrode_sketches = {}
  xenon_sketches = {}
  ptfe_sketches = {}

  # Should produce a cylinder not sure if the line segments will cause a problem
  electrode_sketches["InnerCylinder"] = {
    'pts': [
      ['arc', inner_cylinder_R,0, 0, 0, True],
      ['arc', -inner_cylinder_R,0, 0, 0, True],
      ['line', inner_cylinder_R, 0,],
      ['arc', inner_cylinder_R+cylinder_thicknes, 0, 0, 0, True],
      ['arc', -inner_cylinder_R+cylinder_thicknes, 0, 0, 0, True],
      ['line', inner_cylinder_R+cylinder_thicknes, 0],
    ]
  }
  electrode_sketches["Outer Cylinder"] = {
    'pts': [
      ['arc', outer_cylinder_R,0, 0, 0, True],
      ['arc', -outer_cylinder_R,0, 0, 0, True],
      ['line', outer_cylinder_R, 0,],
      ['arc', outer_cylinder_R+cylinder_thicknes, 0, 0, 0, True],
      ['arc', -outer_cylinder_R+cylinder_thicknes, 0, 0, 0, True],
      ['line', outer_cylinder_R+cylinder_thicknes, 0],
    ]
  }


  # Simulation Volume
  xenon_sketches["LXe"] = {
    'RadialPosition': x0,
    'VerticalPosition': y0,
    'Width': W,
    'Height': H,
  }
  # Dielectric
  xenon_sketches["GXe"] = {
    'RadialPosition': 0.0
    'VerticalPosition': 0.0
    'Radius': outer_cylinder_R
  }
  manual_mapping = {

    }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
