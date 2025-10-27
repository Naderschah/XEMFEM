#ifndef CONSTANTS_H
#define CONSTANTS_H

// Paths
constexpr const char* msh_file = "geometry.msh";

// Device Constants
const bool USE_GPU     = false; // TODO Implement
constexpr int threads  = 20; // n threads to use for meshing and TODO simulating if USE_GPU is off  


// Electrical constants
constexpr double epsilonLXe  = 1.95; 
constexpr double epsilonGXe  = 1.0; 
constexpr double epsilonPTFE = 2.1; 
constexpr double epsilonDefault = 1; // Fallback
constexpr double VAnode   = 1000;             
constexpr double VGate    = 0;       
constexpr double VCathode = -1000.0; 

// Components to use
const bool USE_ANODE   = true; 
const bool USE_GATE    = true;
const bool USE_CATHODE = true;
const bool USE_PTFE    = true;
const bool USE_LXe     = false; // TODO Loosing reference to objects if i cut this out
const bool USE_GXe     = false;

// Geometry / mesh parameters 
constexpr double DriftRegionHeight      = 0.5;
constexpr int    n_wires_anode          = 10;
constexpr int    n_wires_gate           = 11;
constexpr int    n_wires_cathode        = 10;
constexpr double wire_diameter          = 0.01;
constexpr double inner_radius           = 0.2;
constexpr double outer_radius           = 0.25;
constexpr double ring_thickness         = 0.03;
constexpr double anode_wire_height      = 0.45;
constexpr double gate_wire_height       = 0.35;
constexpr double cathode_wire_height    = 0.05;

// Boundary Condition indeces
constexpr int anode_BC_index            = 1000;
constexpr int gate_BC_index             = 1001;
constexpr int cathode_BC_index          = 1002;
// Volume indices
constexpr int TPC_Volume_index          = 2001; // To be removed
constexpr int LXe_Volume_index          = 2002;
constexpr int GXe_Volume_index          = 2003;
constexpr int PTFE_Volume_index         = 2004;


// FEM parameters TODO Where is this used, is this used?
constexpr int fe_order = 3;             // linear finite elements

// Solver and tolerance parameters
constexpr double tol = 1;
constexpr int maxiter = 100000;
constexpr int printlevel = 1;

// DEBUG TODO Not properly usd also add geo vs solver specific
constexpr bool debug     = true; 
constexpr bool QuickMesh = true; // TO quickly generate meshes of poor quality

#endif
