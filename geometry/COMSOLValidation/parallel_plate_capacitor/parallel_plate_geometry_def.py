def build_sketch_dicts(shrinkage_factor):
  dielectric_height = 0.1
  plate_x = -0.25
  plate_width = 0.5
  plate_height = 0.05

  # We define as if a tpc to reuse the existing infrastructure
  electrode_sketches = {}
  xenon_sketches = {}
  ptfe_sketches = {}

  electrode_sketches["TopPlate"] = {
    'RadialPosition':plate_x,
    'VerticalPosition': dielectric_height / 2,
    'Width':plate_width,
    'Height':plate_height,
  }
  electrode_sketches["BottomPlate"] = {
    'RadialPosition':plate_x,
    'VerticalPosition': -dielectric_height / 2 - plate_height,
    'Width':plate_width,
    'Height':plate_height,
  }
  # Simulation Volume
  xenon_sketches["LXe"] = {
    'RadialPosition': 0.0,
    'VerticalPosition': 0,
    'Radius': 3 * plate_width / 2
  }
  # Dielectric
  xenon_sketches["GXe"] = {
    'RadialPosition':plate_x,
    'VerticalPosition': -dielectric_height/2,
    'Width':plate_width,
    'Height':dielectric_height,
  }
  ptfe_sketches["InfiniteSpace"] = {
    'RadialPosition': 0.0,
    'VerticalPosition': 0,
    'Radius': 10 * plate_width / 2
  }

  manual_mapping = {
    1-1: "InfiniteSpace",
    2-1: "LXe",
    3-1: "BottomPlate",
    4-1: "GXe",
    5-1: "TopPlate"
    }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
