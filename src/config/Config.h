#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <mpi.h>

// FIXME REMOVE THE DEFAULTS
// TODO Remove defaults - fix documentation  
// -------------------- Mesh Specific ----------------------------
struct Boundary {
    int bdr_id = -1;        // bdr_id
    std::string type;       // "dirichlet" | "neumann" | "neumann_internal" | "robin"
    double value = 0.0;
    // Used for depth dependent neumann only so far
    bool depth_dependent = false;
    double z_bot = 0.0;
    double z_top = 0.0;
    double value_bot = 0.0;
    double value_top = 0.0;
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
    bool enabled = true;
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

struct CIVSettings {
    std::string method = "ColumnSweep";
    // Random Sample config
    int num_seed_elements = 512;
    // For grid searches
    int nr = 100;
    int nz = 12;
    // For adaptive grid
    int max_levels = 5;
    // For row scanning 
    int block_size = 64;
    // ---------- Not Exposed
    // Integration order
    int    ir_order       = 1;
    // Random Sampling seed - Expose?
    int    rng_seed       = 67; 
};

struct FieldSpreadParams {
    double lower = 0.5;
    double upper = 0.95;
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
    std::string provider  = "VTK";
    std::string method    = "RK34";
    // Relative to the smallest mesh element 
    double c_step         = 2;
    // Tolerance for r,z boundary checks (axis clamp, z-range, etc.)
    double geom_tol       = 1e-12;

    // tracing region
    double r_min = 0.0;
    double r_max = 0.664;
    double z_min = -1.5028;
    double z_max = 0.004;
    // Overwrite for z_max to speed up tracing
    double tracing_z_max;

    // Steps to let electrons traverse n times, if 0 infinite loop
    double max_traversals = 10;
};
// ------------------------- Interpolation --------------
struct Interpolate
{
  int Nx = 1;
  int Ny = 1;
  int Nz = 1;
  bool H1_project;
  // Not Exposed
  bool accept_surface_projection = false;
  // TODO Integrate the below, requires knowledge of mesh dimensions
  double res_x = 0;
  double res_y = 0;
  double res_z = 0;

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

    // Tracing Parameters
    ElectronTraceParams tracing_params;

    // CIV Computation parameters
    CIVSettings civ_params;

    // FieldSpread Parameters
    FieldSpreadParams  field_spread;

    // Interpolation settings
    Interpolate interp;

    // Internal for runtime config checking
    std::string run_mode = "sim";

    // Load from path
    static Config Load(const std::string& path);
    static Config Load(const std::string &path, std::string run_mode);
    // embed config in geometry binary
    static Config LoadFromString(const std::string& yaml_str);
};



void apply_fieldcage_network(Config &cfg);