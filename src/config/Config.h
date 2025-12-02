#pragma once
#include <string>
#include <unordered_map>
#include <vector>

// FIXME REMOVE THE DEFAULTS
// TODO Remove defaults - fix documentation  
// -------------------- Mesh Specific ----------------------------
struct Boundary {
    int bdr_id = -1;        // bdr_id
    std::string type;       // "dirichlet" | "neumann" | "robin"
    double value = 0.0;

};

struct Material {
    int id = -1;            // attr_id
    double epsilon_r = 1.0;
};

struct MeshSettings {
    std::string path = "geometry.msh"; // TODO Does mesh need more settings?
};

// -------------------- Compute / Runtime Settings ----------------------------
struct MPISettings {
    bool enabled = false;
    bool ranks_auto = true;     // true if "auto" was given
    int  ranks = 1;             // ignored if ranks_auto=true
    bool repartition_after_refine = true;
};

struct ThreadsSettings {
    bool enabled = true;
    bool num_auto = true;       // true if "auto" was given TODO is this used
    int  num = 0;               // 0 means "auto" when num_auto=true
    std::string affinity = "compact"; // "compact" | "scatter" | "none"
};

struct DeviceRuntime {
    std::string type = "none";  // "none" | "cuda" | "hip" | "occa" | "cpu"
    bool id_auto = true;        // true if "auto"
    int  id = 0;                // ignored if id_auto=true
    int  per_rank = 1;
};

struct ComputeSettings {
    MPISettings     mpi;
    ThreadsSettings threads;
    DeviceRuntime   device;
};

// If debug enabled but nothing else we default to maximal debug
struct DebugSettings
{
    bool debug                   = false;
    bool dry_run                 = false; // TODO Remove? 
    bool printBoundaryConditions = true;
    bool printHypreWarnings      = true;

    // Dump data to file where applicable
    bool dumpdata              = false;

    // If true, TraceElectronFieldLines (Optimization: CIV Computation) should:
    //   - only trace seed with index debug_single_seed_index,
    //   - optionally override c_step for that seed to a very small value.
    bool   debug_single_seed         = false;
    int    debug_single_seed_index   = -1; // Internal dont touch
    double debug_c_step_override = 0.; // Overwrite c_step in debug mode
};

// -------------------- Circuit Computation --------------------------
struct FieldCageNode {
    std::string name;      // Name to reference by 
    std::string boundary;  // Associated Boundary condition
    bool fixed = false;    // true if potential is fixed by the boundary condition (ie cathode & top field shaping)
};

struct FieldCageEdge { 
    std::string n1;      // node name from
    std::string n2;      // node name to
    std::string R_name;  // resistance from n1 to n2 
};

struct FieldCageNetworkSettings {
    bool enabled = false;
    std::vector<FieldCageNode> nodes; // BC 
    std::vector<FieldCageEdge> edges; // Resistors between two nodes
    std::unordered_map<std::string, double> R_values; // Named resistors 
};


// -------------------- Sweep description ----------------------------
struct SweepEntry { 
    enum class Kind { Discrete, Range };

    std::string name;   // optional label
    std::string path;   // dot-path into YAML, e.g. "boundaries.name.value"

    Kind kind = Kind::Discrete;

    std::vector<std::string> values; // Discrete interface

    // Range interface
    double start = 0.0;
    double end   = 0.0;
    int    steps = 0;   // inclusive
};

// ----------------- Optimization description ------------------------
// Currently no support for discrete variables implemented
struct OptimizeVar {
    std::string name;   // label
    std::string path;   // dot path into Config, like "solver.order" or "boundaries.BC_Anode.value"

    // For continuous (bounded) variables:
    double lower = 0.0;
    double upper = 0.0;
    double initial = 0.0;       // starting guess; if 0, you can derive from lower/upper

    // For discrete variables:
    std::vector<std::string> values;  // a small, finite set of allowed values
};
struct OptimizationSettings {
    bool enabled = false;
    bool metrics_only = false; // Triggers no mesh simulation but computes metrics for all available simulations

    bool print_results = false;

    // Optimization target
    std::string objective;

    // Evaluation domain
    double r_min          = 0.0;
    double r_max          = 0.664;
    double z_min          = -1.5028;
    double z_max          = 0.004;

    // Variables to tune:
    std::vector<OptimizeVar> variables;

    // Nelder-Mead Settings
    int    max_iters      = 200;
    int    max_fun_evals  = 200;
    double tol_fun        = 1e-3;
    double tol_x          = 1e-3;
    bool   adaptive       = true;
};
// -------------------- Compute / Runtime Settings ----------------------------
struct SolverSettings {
    // MFEM / solve controls
    bool axisymmetric = false;
    int order = 3;
    std::string assembly_mode = "partial";

    double atol = 1.0;
    double rtol = 0.0;
    int    maxiter = 100000;
    int    printlevel = 1;
};
// Field Line tracing parameters 
struct ElectronTraceParams
{
    // Select CIV tracing method 
    std::string method = "ColumnSweep"; // Or RandomSample or currently not working InformedSweep
    // Maximum number of accepted steps per electron
    int    max_steps      = 200;
    // Tolerance for r,z boundary checks (axis clamp, z-range, etc.)
    double geom_tol       = 1e-12;
    // Order for integration rules used when seeding CIV / volume weights
    int    ir_order       = 1;
    // Field-line step control (first-order drift integrator)
    // Base fraction of local element size used as nominal step:
    //   ds_nominal = c_step * hK
    double c_step         = 0.25;
    // Adaptive bounds (for error control):
    //   ds_min = ds_min_factor * hK; ds_max = ds_max_factor * hK
    double ds_min_factor  = 0.001;  // small steps near strong curvature
    double ds_max_factor  = 0.30;  // do not exceed frac of hK 
    // Error-control parameters for embedded RK
    double adapt_grow     = 1.5;   // max step growth per acceptance
    double adapt_shrink   = 0.25;   // min shrink factor on rejection
    double tol_rel        = 1e-3;  // relative position error tolerance
    // tracing maximum
    double tracing_z_max = 0.004;
    // If true, treat r <= geom_tol as an exit condition 
    bool   terminate_on_axis = true;
    int    num_seed_elements = 0;  // <= 0 means "use all elements"
    int    rng_seed = 67;       // We always want the same one - Currently not exposed

    // Parameters specific to informed sweep
    bool dump_civ_boundary = true;

    // TODO This should not be here 
    // It is only kept as i need to clean up CIV Tracing
    double r_min          = 0.0;
    double r_max          = 0.664;
    double z_min          = -1.5028;
    double z_max          = 0.004;
};

// ------------------------- Actual Config Struct -------------------------------
struct Config {
    int schema_version = 1;
    std::string geometry_id;

    std::string save_path; 
    bool delete_files_present = false;

    MeshSettings mesh;
    ComputeSettings compute;
    DebugSettings  debug;
    SolverSettings solver;

    std::unordered_map<std::string, Boundary> boundaries;
    std::unordered_map<std::string, Material> materials;

    // Sweep Config options
    std::vector<SweepEntry> sweeps;

    // Optimization Config Values
    OptimizationSettings optimize;

    // Configuration to solve circuits
    FieldCageNetworkSettings fieldcage_network;

    // TODO Not yet configuable (edit .cpp file)
    ElectronTraceParams tracing_params;

    // Load from path
    static Config Load(const std::string& path);
    // embed config in geometry binary
    static Config LoadFromString(const std::string& yaml_str);
};



void apply_fieldcage_network(Config &cfg);