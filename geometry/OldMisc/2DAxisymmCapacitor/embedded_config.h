#pragma once
namespace embedded_config {
// Note the custom raw-string delimiter to avoid accidental terminators.
inline constexpr const char* kConfigYaml = R"EMBED(
schema_version: 1
geometry_id: 2DCapacitor

# Where the mesh is/should be produced
# FIXME Unused at the moment 
mesh:
  path: geometry.msh 


# Processing device TODO not implemented and will get extended
device:
  use_gpu: false
  threads: 20

debug:
  debug: false
  quick_mesh: false

# solver configuration
solver:
  order: 3
  atol: 1e-12
  rtol: 0.0
  maxiter: 100000
  printlevel: 1
  mesh_save_path: "simulation_mesh.msh"
  V_solution_path: "solution_V.gf"
  Emag_solution_path: "solution_Emag.gf"


# Geometry specifics
materials:
  dielectric:
    attr_id: 2002
    epsilon_r: 1.95
  air:
    attr_id: 2003
    epsilon_r: 1.0

boundaries:
  TopPlate:
    bdr_id: 1000
    type: dirichlet
    value: 1000.0
  BottomPlate:
    bdr_id: 1001
    type: dirichlet
    value: 0.0
  OuterBoundary:
    bdr_id: 1002
    type: dirichlet
    value: 0


)EMBED";
} // namespace embedded_config
