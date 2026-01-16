# ELectrostatics with MFEM

Generate Geometry:

run Salome docker image (see ./geometry), generate mesh, run the postprocessor, configure the yaml. 

Compile solver: 
```bash
mkdir build
cd build
cmake ../src && make 
```
Run solver 
```bash
./SWEEP --config /path/to/config.yaml
```

## Simulation configuration

The configuration is written as a yaml file. 

Not all options are used by all portions of the code, some options are meshing specific. Most are MFEM solver specific. 

An overview of the currently supported options is given below
```yaml
# Has no use here yet 
schema_version: 1
# Name used in metadata for simulations only 
geometry_id: SomeName 

mesh:
  # No need to modify this the mesh postprocessor writes the correct path
  path: "/home/felix/work/geometry/mesh/mesh22.msh"
  # Which geometry to generate in meshing, currently only SR3 is supported
  geometry: "SR3"
  # Wheter or not to autoappend the generated boundaries and attributes to the currently used config file
  autoappend: "True"

# Where Simulation results should be saved
save_path: "/home/felix/work/sim_results"
# Wheter to overwrite existing files (trashes the directory if true) - meta file is always overwritten
delete_files_present: true 

# Debug Settings
debug:
  # TODO: Make exhaustive list of what this does
  debug: false
  # TODO: Make exhaustive list of what this does
  dry_run: false
  # Print boundary conditions as loaded in MFEM simulation
  printBoundaryConditions: false
  # Wheter or not to print hypre warnings (there is always 1+ as Hypre misinterprets some solver internals)
  printHypreWarnings: false
  # TODO: Exhaustive list
  dumpdata: false
  # TODO: Exhaustive list 
  debug_single_seed: false

# Compute Settings
compute:
  # Multi processing Interface
  # Used for multimachine computation 
  # TODO Currently untested -> Config will need extending for this - might also be needed for multi cpu devices
  mpi: 
    enabled: false 
    ranks: auto 
    repartition_after_refine: true
  # Multithreading
  threads:
    enabled: true
    # Automatically assign number of threads based on what OpenMP sees as thread number for CPU 
    num: -1 
    # scatter or compact, how memory is to be managed between threads 
    affinity: scatter
  
  # Currently not respected
  # Will be once GPU and MPI is figured out 
  device: 
    type: none
    id: auto 
    per_rank: 1

# Sweep configuration 
# Sweeps performed if sweeps is non empty
# Represented as list of items to sweep over 
sweeps:
    # Some Name for the sweep 
  - name: "AnodeSweep"
    # Path inside the YAML from root (anything can be sweeped)
    path: "boundaries.BC_Anode.value" 
    # Specify a range sweep
    kind: "range"
    # Starting at 
    start: 1900.0
    # Ending At 
    end:   4900.0
    # in number steps 
    steps: 4
  # Or a discrete Sweep 
  - name: "order_sweep"
    path: "solver.order"
    kind: "discrete"
    values: [1, 2, 3]

# Optimization Configuration (minimize some metric relative to some parameters)
# Uses Nelder-Mead
optimize:
  enabled: false
  # TODO What does this do again
  metrics_only: false
  # Wheter to print per optimization results
  print_restults: true
  # Objective function
  # Available are: FieldSpread, CIV, self_weighting
  objective: "self_weighting"
  # Number of evaluations per iteration
  max_evals: 1
  # Number of iterations 
  max_iters: 1
  # TODO: What is this again 
  rel_tol_f: 1.0e-3
  # TODO: What is this again 
  rel_tol_x: 1.0e-3
  # TODO: What is this again 
  adaptive: true 

  # TPC dimension definition (parts of the detector we care about)
  r_min: 0.0         
  r_max: 0.664       # PanelRadialPosition
  z_min: -1.5025     # CathodeVerticalPosition + CathodeWireDiameter
  z_max: 0.004

  # Electron Tracing Configuration (Charge Insensitive Volume)
  # TODO Update with good parameters 
  # TODO Add proper docs
  trace_params:
    # Which propagation method to use
    # Euler-Cauchy, RK23, RK45
    method: "RK45"
    # The step size = c_step * hk, where hk is the shortest vertex to vertex distance inside the current mesh element 
    c_step: 0.25
    # these * hk produce the smallest and maximal step sizes 
    ds_min_factor: 0.05
    ds_max_factor: 2 
    # Step allowed if error_metric / hK  < tol_rel
    tol_rel: 0.1
    # Adaptive step resizing, factor by which estimate drops/grows
    adapt_shrink: 0.7
    adapt_grow: 1.2
    # Boundary check tolerance 
    geom_tol: 1e-5
    # Maximum integration steps not counting retries due to too high error
    max_steps: 1000000

  civ_params:
    # CIV method: RandomSample, ColumnSweep, InformedSweep
    method: "RandomSample" 
    # Samples to use for method RandomSample
    num_seed_elements: 512
    # Slices to take for ColumnSweep
    n_slices: 24
    # Wheter the defined boundary line is to be written to file or not
    dump_civ_boundary: false    
    # Informed Sweep: For initial column search define an extra point at the 
    # bottom min_col_pos away from the bottom, if all tested elements are charge 
    # sensitive the procedure terminates concluding there is no CIV (ie = 0)
    min_col_pos: 0.001


  # TODO Implement - whats the status on this
  fieldSpread_params:
    TopBottomDist: 0.14 
    BottomPercentile: 0.5
    UpperPercentile: 0.95 

  # Variables to optimize against
  # represented as a list,
  # TODO Is discrete supported?
  variables:
    - name: "BottomScreenVoltage"
      path: "boundaries.BC_BottomScreen.value"
      kind: "continuous"
      lower: -3000.0
      upper: 3000.0
      initial: 0.0

# Electrostatics Solver Configuration 
solver:
  # Run axisymmetric Simulation
  axisymmetric: true
  # Solution order
  order: 3
  # TODO  https://mfem.org/howto/assembly_levels/
  assembly_mode: partial              # full | partial | matrix_free
  # Solution absolute tolerance
  atol: 1e-50 #1e-25
  # Solution relative tolerance
  rtol: 0.0
  # Maximum number of iterations 
  maxiter: 100000
  # Should the solver print the status on each step (a lot of text)
  printlevel: 0
  # Names to save mesh and solutions under
  # TODO Is this still respected?
  mesh_save_path: "simulation_mesh.msh"
  V_solution_path:  "solution_V.gf"
  Emag_solution_path: "solution_Emag.gf"

# Plotting Options
# TODO Whats the status on this
plotting:
  show_mesh: false
  mesh_edge_width: 0.1
  img_width: 2160
  img_height: 3840


# Rest is autogenerated by mesh postprocessor
# Volumes
materials:
  # Some Name 
  GXe:
    # Internal attribute ID (inside the GMSH file)
    attr_id: 1
    # Relative permativity
    epsilon_r: 1.0
  # And so on
# Equipotential Surfaces
boundaries:
  # Some Name 
  BC_AllPMTs:
    # Internal boundary ID (inside the GMSH file)
    bdr_id: 4
    # Type of boundary condition (we only need dirichlet - TODO Is anything else currently supported?)
    type: dirichlet
    # At Voltage
    value: -1300
  # And So on 
  
# Voltages can also be autocomputed as a funcion of others
# Used only for field cage (as the name implies)
fieldcage_network:
  enabled: true
  # Equipotentials on the circuit 
  nodes:
    # List of 
    - name: "BC_FieldShapingRings0"
      # What boundary it corresponds to (same name as above)
      boundary: "BC_FieldShapingRings0"
      # If it is a fixed point 
      fixed: true
    - name: "BC_FieldShapingRings1"
      boundary: "BC_FieldShapingRings1"
      fixed: false
    # And so on 
  # Resistor Values in the network 
  resistors:
    R1: 1250000000.0
    R2: 2500000000.0
    R3: 5000000000.0
    R_C: 7000000000.0
  # How nodes are connected, list again 
  edges:
      # What is in between two nodes 
      # N1 - [R] - N2 
    - n1: "BC_FieldShapingRings0"
      n2: "BC_FieldShapingRings1"
      R: "R1"
    - n1: "BC_FieldShapingRings1"
      n2: "BC_FieldShapingRings2"
      R: "R1"

```
## Results Inspection

Plots can be made with the plot subcommand, point to the config file and everything will be plotted. 

Alternatively, the pvdu file found under Simulation in the results path can be opened in paraview (has a huge amount of viewing and postprocessing options) for inspection, open it from the terminal. Or in glvis using `glvis -m simulation_mesh.msh -g V.gf -k "AmaagcmRj"` which views the internally used file directly. 

## Software Stack

#### Electrostatics 
MFEM: 
- Version 4.8.0 
- Highest currently available

Hypre:
- version 2.33.0 
- version geq 3.0.0 does not work with MFEM 
- Previous versions 
  - tried with 2.32.0 : HypreBoomer AMG reports 1 OpenMP thread regardless of set number. Either checks in a funny way, or bugged  

TODO: GPU 


#### Meshing 

SALOME:
- Version 9.9.0
- Anything above this version should work, one version down does not due to a function used in tagging

TODO: Script should be available headless, and the two containers merged to make it possible to mesh and simulate in one shot 


# TODO

## Geometry 

Implement SR2 Geometry

Allow Headless running, and run postprocess directly

Merge docker containers

Netgen settngs inside config.yaml

## MFEM

use pragma stop with ifndef 

Also implement misc boundary conditions

GPU implementation & MPI testing (fairly certain wont work at the moment)

Need to fix indentation 

Cut down on the structs 

Generalized Path Handler 

Merge Plotting and Simulation (or did i do this already)

Add interpolated export

### Config File
Stubs to clean up:
- optimize.metrics_only : No longer used (once sweep executable gone)

### Optimization

Clean up dispatch logic 

CIV - Find good hyperparameters 

CIV Compare to electric field direction as proxy 

Implement difference of parameters as fixed




Field Optimization Strategy:

Two sets of metrics:
- One for Drift 
- One for Extraction 

Optimize Both seperately for some fixed parameters relative to themselves, and then check the superposition, then vary a little and see if this is indeed the best config for the combined


First steps figure out CIV - add injection point where simulation results exist and we simply compute 


#### Tracing

There are a number of functions that may be able to be merged/seperated or may be the same thing

FindElementForPointLocal
FindElementForPointGlobal
FindElementForPointRobust
TraceSingleCIProbe

Also all old methods should be updated and the code between them should be more shared


The informed sweep method should work but the initial column sweep fails due to some 
charge insensitive elements close to the wall -> For now addding wall surface charge
to eliminate this behavior to determine operating conditions

Once this is done I should investigate this and find a more suitable method.