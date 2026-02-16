# Electrostatics with MFEM

## Basic Usage
Generate Geometry:

- run Salome docker image (`./run_docker.sh` in geometry), 

Start salome on the command line, type Ctrl+t, run `build_geometry.py`, wait until the GUI becomes response, and then salome can be closed. 

Switch to the main docker file, and in geometry run: `python3 postprocess_mesh.py`, note this will rewrite the BC's in the config file and delete all comments present. Votages are sourced from the respective `<Prefix>_constants.py` file where prefix is configured in the yaml. The same goes for the geometry. 

Currenlty implemented geometries are:
- SR3 

Currently implemented Voltages:
- design
- SR0
- SR1
- SR2

Note that all currently use the same numbers for wall charge up, so for earlier runs this feature should be disabled (set the Neumann BC to null value or remove the list item entirely)

To solve the mesh.

Compile XEMFEM: 
```bash
mkdir build
cd build
cmake ../src && make 
```
And Run 
```bash
./XEMFEM sim 
```

config path can be specified but if you didnt change it all defaults should work.

Alongside this the following options are available:
```bash
 # Runs a basic simulation, a sweep if sweep parameters are configured, 
 # or a Nelder-Mead Optimization if parameters are configred (Optimization has priority over sweep)
./XEMFEM
./XEMFEM sim 
# Computes metrics used for optimization (currently CIV and Field homogeneity)
# for the last saved mesh and results files
./XEMFEM metrics
# Produces V and E Plots, Configuration of these via config.yaml not yet implemented
./XEMFEM plots

 # The config.yaml path defaults to geometry/config.yaml but can be provided with 
 # --config path/to/config.yaml

```

## Features

Geometry is created via points lists that specify the 2D geometry, these are constructed directly from slices of the master CAD file. Some cleanup has to be applied after the fact, if you have the same CAD file (I used XENT-TPC_20250428.STEP) you should get the exact same output, notably any other format is not guaranteed to work due to the way CAD file formats label their item trees and repeated items. See geometry/createGeometryFromCAD/README.md for more information. 

2D, 2D axisymmetric and 3D simulations are implemented, higher D is possible within MFEM, however, I am not sure the code is general enough to handle this. 3D also currently has limited support, it should work fine with the main simulation suite (haven't tried in a while however), but the metrics computation and plotting do not currently support higher D. 

The field Cage is written as a chain of edges (wires) and nodes (resistors), these are configured by the postprocess_mesh function. Prior to each simulation these are recomputed, agreement between the python script and MFEM implementation was verified. This is only relevant for sweeps.

Multithreading is implemented and supported, GPU support will (likely) be shipped in a seperate branch due to the extra heavy dependencies. OMP (threading) is active by default, so threads = 1 corresponds to the closest thing to single threaded available. MPI is active by default as I originally intended to fully support it, however, the config still must be heavily expanded to support this, alongside this not all code paths are MPI safe, this will eventually be added on the main release branch.

Field Homogeneity is a trivial computation and functions well, it is in its final form, barring code cleanup and possible minor refactoring for MPI and reducing computational cost. Electron tracing prooved nontrivial, with poorly conditioned meshes (TODO What exactly is currently the reason) producing circiling motion in the least variable region of the TPC (close to r = 0 central height range), hence a hard error is thrown in case of Timeout. VTK is the most stable integrator avaialble, however, it is intended for plotting so recording paths is mandatory, the alternate MFEM native implementation is less stable and will not integrate in the inner portion of the TPC, which will be problematic for some CIV strategies. 


## Simulation configuration

The configuration is written as a yaml file. 

Not all options are used by all portions of the code, some options are meshing specific and do not impact MFEM. 
Care should be taken not to modify any meshing specific fields during optimization to avoid result confusion.

An overview of the currently supported options is given below
```yaml
# Has no use here yet 
schema_version: 1
# Name used in metadata for simulations only 
geometry_id: SomeName 

mesh:
  # No need to modify this the mesh postprocessor writes the correct path
  path: "/work/geometry/mesh/mesh22.msh"
  # Which geometry to generate in meshing (uses <geometry>*.py for finding the file)
  geometry: SR3
  # Which set of voltages to use, also defined in the same file (mesh_postprocessor uses this)
  voltages: SR3
  # Wheter or not to autoappend the generated boundaries and attributes to the currently used config file -> Otherwise BC's and materials need to be copied from the config_autogen file 
  autoappend: true
  # Wheter to append fieldshaping (ring and guard) voltages to the generated config file
  # they are appended for the final written under save_path, but are not mandatory for 
  # the source config since the field cage network is solved later either way
  append_fieldShaping_Voltages_toConfig: false

  # Meshing Parameters (used in SALOME only)
  # If debugCoarse is given all else ignored else all are required
  NetgenParams:
    # For quick meshing (significantly reduced number of nodes)
    debugCoarse: false
    # Preset to go off, most if not all settings will be overwritten by below
    preset: VeryFine
    # Maximum Element size allowed (don't make this too large)
    maxSize: 0.015008
    # Minimum allowed element size (0 = it will pick 1e-3m), smallest feature is 2.16e-5 (Wires)
    minSize: 1.0e-6
    # How quickly elements may change size in adjacency space
    growthRate: 0.1
    # Use Surface Curvature as element size driver (ie CAD driven)
    UseSurfaceCurvature: 1
    # Number of segments per topological edge 
    NBSegmentsPerEdge: 4
    # Number of segments for a given radius ie acceptable mesh length is h = R/N where R is the radius and N is the below parameter. So for 12 edges along a circle we want N = 1.91, N=3 is approximately 19 edges
    NBSegmentsPerRadius: 3
    # Instead of Surface Curvature use chordal error (maximal deviation from CAD)
    UseChordalError: 0
    # Deviation allowed (meters)
    ChordalErrorValue: 1
    # Advanced Quality Algorithms 
    # Weight of mesh "quality" as compared to accuracy (0.5 both are equally important)
    AQ_ElemSizeWeight: 0.5
    # delauney vs advancing front algorithms 
    AQ_UseDelauney: 1
    # Use Strict overlap checks (overlapping elements = bad)
    AQ_CheckOverlapping: 1
    # Netgen splits surfaces into parametric shapes, if 1 it checks that the boudaries align properly 
    AQ_CheckChartBoundary: 1
    # Fuse duplicate edges (duplicat edges = bad)
    AQ_FuseEdges: 1
    # Allow quadrilateral elements (rather than triangles) dont know how well XEMFEM in its current form handles quadrilaterals (I think tracing will break), but its also not needed for our purposes
    Opt_QuadAllowed: false
    # if true a midpoint node is added to each edge and each edge line becomes curved (quadratic interpolation) - also overkill, also will likely break particle tracing (allthough it might work not 100% sure)
    Opt_SecondOrder: 0
    # Post meshing optimization
    Opt_Optimize: true

  # Adaptive Mesh Refinement (Indirectly Physics Driven) - XEMFEM
  # basic idea is solve the mesh, estimate a per element error, subdivide high error elements and coarsen low error elements. 
  AMR:
    # Enable it -> This will significantly increase computational time
    enable: true
    # maximal number of refinement/derefinement iterations
    max_iter: 10
    # Error estimate used to compute element-wise error
    # Kelly : Edge/flux jump based (robust for Poisson equation)
    #     L2 norm of the normal flux jump across faces 
    #     Each element has a different representation for its boundary, 
    #     this metric essentially bounds by how much the normal flux 
    #     (delta V dot n) may disagree between two adjacent elements
    #     at the boundary 
    #     In 2D it has units of Volts (in 3D Volts sqrt(m))
    # ZienkiewiczZhu : recovered graident estimator
    # L2ZienkiewiczZhu : Smoothed L2 projection variant
    estimator: Kelly
    # Absolute per element error tolerance, measured in the norm of the estimator used
    local_error_goal: 1e-3
    # Absolute gloabl error tolerance (p-norm of element errors), 0.0 for disabled
    total_error_goal: 0.0
    # Controls refinement aggression when using global error marking, 0.0 is purely local threshhold driven
    total_error_fraction: 0.0
    # p-value to compute global total error from element errors (1.0 -> L1 norm, 2.0 -> L2 norm) the larger the value the stronger peak element errors are emphasized
    total_error_norm_p: 2.0
    # Wheter or not to enable derefine (essentially do we care about elements having a smaller error)
    enable_derefine: true
    # Control parameter to determine if an element is eligible for coarsening, must be leq 1.0, so if an element has error e/allowed_err < derefine_hysterisis it will be considered for coarsening
    derefine_hysterisis: 0.5
    # How child errors are combined when deciding coarsening, 0->MinChildError, 1->SumOfChildErrors, 2->MaximumChildError
    derefine_op: 2
    # Cap on total number of elements, 0 is disabled
    max_elements: 0 
    # Maximal degrees of freedom (1 Triangle 2D = 3 DOF), 0 is disabled
    max_dofs: 0
    # Maximum allowed nonconforming refinement level difference between neighbors, 0 unlimited
    nc_limit: 0
    # Wheter to prefer conforming vs nonconforming refinement (for quad/hex this should be false)
    prefer_conforming_refinement: true
    # minimal number of elements that must be marked for refinement before loop early terminates
    min_marked_elements: 1
    # Minimum relative reduction in total error between iterations (fraction not percentage), 0 is disabled 
    min_rel_error_reduction: 0.0
    # Rebalance parallel mesh afterwards (redistribute over MPI nodes - if MPI it should be true, if not MPI it will not do anything so it should be true)
    enable_rebalance: true
    # Print AMR diagnostics
    verbose: false

# Where Simulation results should be saved
save_path: "/work/sim_results"
# Wheter to overwrite existing files (trashes the directory if true) - meta file is always overwritten
delete_files_present: true 

# Debug Settings
debug:
  # General debug prints
  debug: false
  # Executes Sweeps (and Optimization) without running simulation such that settings can be verified
  dry_run: false
  # Print boundary conditions as loaded in MFEM simulation (mesh gen sanity check)
  printBoundaryConditions: false
  # Wheter or not to print hypre warnings (there is always 1+ as Hypre misinterprets some solver internals)
  printHypreWarnings: false
  # Saves pathlines in particle tracing (if multiple sets are traced overwrites the old, best used with RandomSample CIV strategy)
  dumpdata: false

# Compute Settings
compute:
  # Multithreading (OMP)
  threads:
    # Number of Threads: -1 Automatically assign number of threads based on what OpenMP sees as thread number for CPU 
    num: -1 
    # scatter or compact, how memory is to be managed between threads 
    affinity: scatter

# Sweep configuration 
# Sweeps performed if sweeps (list) is non empty
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
  # Or a discrete Sweep (mesh sweeps are currently not possible TODO Should be easy modification)
  - name: "order_sweep"
    path: "solver.order"
    kind: "discrete"
    values: [1, 2, 3]
  # Or a sweep over fixed configurations (ie solve for all of these)
  - name: "Preset"
    kind: "fixed"
    configs:
      - label: "GateIsGround"
        set:
          boundaries.BC_GateElectrode.value: 0.0
          boundaries.BC_AnodeElectrode.value: 1000
      - label: "AnodeIsGround"
        set:
          boundaries.BC_GateElectrode.value: 1000
          boundaries.BC_AnodeElectrode.value: 0.0

# Optimization Configuration (minimize metric relative to some parameters)
optimize:
  enabled: false
  # Legacy executable setting
  metrics_only: false
  # Wheter to print per optimization results
  print_restults: true
  # Objective function
  # Available are: FieldSpread, CIV, self_weighting 
  # TODO Extra config block for drift v extraction computations
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

  # TPC ROI (This is Drift Region in SR3)
  r_min: 0.0         
  r_max: 0.664
  z_min: -1.5025
  z_max: 0.004

  # Electron Tracing Configuration (used for Charge Insensitive Volume)
  trace_params:
    # Which tracing provider is to be used 
    # VTK is from a plotting library (resistant to poor meshes)
    # XEMFEM is custom tracing implementation, sensitive to poor meshes but much quicker
    provider: VTK
    # Which stepping method to use, VTK support RK4 and Euler-Cauchy 
    # XEMFEM supports Euler-Cauchy, RK4, RK45 (Dormand-Prince 5), and RK54 (Cash-Karp), all pair methods do the standard error based step inc/dec logic
    method: "RK4"
    # step_size = c_step * smallest mesh triangle length
    c_step: 1e-4
    # Geometric tolerance 
    geom_tol: 1e-6
    # Particle may traverse drift region max_traversals times (max steps is computed from this and c_step)
    max_traversals: 2

  civ_params:
    # CIV method: Debugging methods include RandomSample, Grid, AdaptiveGrid, RowSweep. 
    method: "RandomSample" 
    # Samples to use for method RandomSample
    num_seed_elements: 512
    # grid divisions for grid based methods
    nr: 42
    nz: 126
    # Levels to reduce resolution in adaptive grid search
    max_levels: 2


  # In case other quantiles are desired for field spread
  fieldSpread_params:
    BottomPercentile: 0.05
    UpperPercentile: 0.95 

  # Variables to optimize against
  # represented as a list, only continuous implemented
  variables:
    - name: "BottomScreenVoltage"
      path: "boundaries.BC_BottomScreen.value"
      kind: "continuous"
      lower: -3000.0
      upper: 3000.0
      initial: 0.0

# Electrostatics Solver Configuration 
solver:
  # Run axisymmetric Simulation (No checks of the mesh are performed)
  axisymmetric: true
  # Solution order
  order: 3
  # full | partial | matrix_free :   https://mfem.org/howto/assembly_levels/
  assembly_mode: partial
  # Solution absolute tolerance
  atol: 1e-50
  # Solution relative tolerance
  rtol: 0.0
  # Maximum number of iterations 
  maxiter: 100000
  # Should the solver print the status on each step (a lot of text)
  printlevel: 0

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
    # dirichlet boundary condition enforces constant potential on surface
    type: dirichlet
    # At Voltage
    value: -1300
  # Some name for a surface charge
  BC_SurfaceCharge:
    bdr_id: 3
    # Neumann boundaries impose flux discontinuity over surface
    type: neumann 
    # If false applied constant over surface, otherwise linearly z dependent 
    depth_dependent: true
    # bottom z coordinate
    z_bot: -1.50016
    # Top z coordinate
    z_top: -0.0008
    # C/m at the bottom 
    value_bot: -1
    # C/m at the top
    value_top: -5
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

Plots can be made with the plot subcommand.

Alternatively, the pvdu file found under Simulation in the results path can be opened in paraview for inspection (has a huge amount of viewing and postprocessing options), open it from the terminal. Or in glvis using `glvis -m simulation_mesh.msh -g V.gf -k "AmaagcmRj"` which views the internally used file directly (MFEM native). 

## Software Stack

#### Meshing 

SALOME:
- Version 9.9.0
- Anything above this version should work, below does not due to the post partition tagging logic.

Note that sometimes the determinism of Salome in partitioning gets lost causing renaming indeces to mess up, just rerun the script and it will work.
This is always true in the graphical session on the second partition as the internal state doesnt properly reset. 

#### Electrostatics 
MFEM: 
- Version 4.8.0 

Hypre:
- version 2.33.0 
- version geq 3.0.0 does not work with MFEM 4.8
- Previous versions 
  - tried with 2.32.0 : HypreBoomer AMG reports 1 OpenMP thread regardless of set number. Either checks in a funny way, or bugged  

VTK:
- version 9.5.0
- Hardcoded opengl2 backend with OSMesa (v6) and no X (assumes no GPU)

GSLib:
- version 1.0.9


# Debugging notes

To debug MPI sections
```bash
mpirun -n 2 --allow-run-as-root xterm -hold -e gdb -ex run --args ./XEMFEM metrics
```

TO run in mpi: mpirun -np 4 --allow-run-as-root ./XEMFEM