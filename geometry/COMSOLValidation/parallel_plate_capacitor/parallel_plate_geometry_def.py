def build_sketch_dicts(shrinkage_factor):
  dielectric_height = 0.1
  plate_x = 0.5
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
    'VerticalPosition': -(dielectric_height + plate_height * 2 + dielectric_height) * 1.5 /2,
    'Width': plate_x * 2 + plate_width,
    'Height': (dielectric_height + plate_height * 2 + dielectric_height) * 1.5 ,
  }
  # Dielectric
  xenon_sketches["GXe"] = {
    'RadialPosition':plate_x,
    'VerticalPosition': -dielectric_height/2,
    'Width':plate_width,
    'Height':dielectric_height,
  }
  manual_mapping = {
        1-1: "BottomPlate_part",
        2-1: "TopPlate_part",
        3-1: "GXe_part",
        4-1: "LXe_part",
    }
  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
