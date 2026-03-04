def build_sketch_dicts(shrinkage_factor):
  wire_center_x = 0.2
  wire_center_y = 0.0

  r1 = 0.1
  r2 = 0.15

  W = 2*wire_center_x +r1 + r2 + 2*max(r1, r2)
  x0 = - W/2
  H = max(r1, r2)*2 * 3/2
  y0 = -H/2

  # We define as if a tpc to reuse the existing infrastructure
  electrode_sketches = {}
  xenon_sketches = {}
  ptfe_sketches = {}

  electrode_sketches["Wire1"] = {
    'RadialPosition': wire_center_x,
    'VerticalPosition': wire_center_y,
    'Radius': r2
  }
  electrode_sketches["Wire2"] = {
    'RadialPosition': -wire_center_x,
    'VerticalPosition': wire_center_y,
    "Radius": r1
  }
  # Simulation Volume
  xenon_sketches["LXe"] = {
    'RadialPosition': 0,
    'VerticalPosition': 0,
    'Radius': 3*wire_center_x,
  }
  # Dielectric
  xenon_sketches["GXe"] = {
    'RadialPosition':-wire_center_x,
    'VerticalPosition':-min(r1, r2),
    'Width':2*wire_center_x,
    'Height':2*min(r1, r2),
  }
  ptfe_sketches["InfiniteSpace"] = {
    'RadialPosition': 0.0,
    'VerticalPosition': 0,
    'Radius': 10 * wire_center_x 
  }

  manual_mapping = {
    1-1: "InfiniteSpace",
    2-1: "LXe",
    3-1: "Wire2",
    4-1: "Wire2",
    5-1: "GXe",
    6-1: "Wire1",
    7-1: "Wire1"
    }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
