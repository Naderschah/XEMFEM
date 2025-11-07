#pragma once
#include <string>
#include <unordered_map>

// FIXME REMOVE THE DEFAULTS
// TODO Remove defaults - fix documentation  
// -------------------- Core Physics ----------------------------
struct Boundary {
    int bdr_id = -1;        // bdr_id
    std::string type;       // "dirichlet" | "neumann" | "robin"
    double value = 0.0;
};

struct Material {
    int id = -1;            // attr_id
    double epsilon_r = 1.0; // Epsilon
};

// -------------------- Compute / Runtime Settings ----------------------------
struct MPISettings {
    bool enabled = false;
    bool ranks_auto = true;     // true if "auto" was given
    int  ranks = 1;             // ignored if ranks_auto=true
    std::string hostfile;       // optional
    bool repartition_after_refine = true;
};


struct ThreadsSettings {
    bool enabled = true;
    bool num_auto = true;       // true if "auto" was given
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
    std::string preset;         // "", "serial", "threads", "mpi", "gpu", "mpi+gpu", "mpi+threads"
    MPISettings     mpi;
    ThreadsSettings threads;
    DeviceRuntime   device;
};


struct DebugSettings {
    bool debug = false;
    bool quick_mesh = false;
};

struct MeshSettings {
    std::string path = "geometry.msh";
};


// -------------------- Compute / Runtime Settings ----------------------------
struct SolverSettings {
    // MFEM / solve controls
    bool axisymmetric = false;
    int axisymmetric_r0_bd_attribute = 9999;
    int    order = 3;
    std::string assembly_mode = "partial";
    std::string solver = "pcg";
    std::string precond = "bommerang"; 

    double atol = 1.0;
    double rtol = 0.0;
    int    maxiter = 100000;
    int    printlevel = 1;

    // Outputs
    std::string mesh_save_path   = "simulation_mesh.msh";
    std::string V_solution_path  = "solution_V.gf";
    std::string Emag_solution_path = "solution_Emag.gf";
};

struct Config {
    int schema_version = 1;
    std::string geometry_id;

    MeshSettings mesh;
    ComputeSettings compute;
    DebugSettings  debug;
    SolverSettings solver;

    std::unordered_map<std::string, Boundary> boundaries;
    std::unordered_map<std::string, Material> materials;

    // Load from path
    static Config Load(const std::string& path);
    // embed config in geometry binary
    static Config LoadFromString(const std::string& yaml_str);
};

