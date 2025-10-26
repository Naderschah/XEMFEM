#ifndef CONSTANTS_H
#define CONSTANTS_H

// Paths
constexpr const char* msh_file = "geometry.msh";

// Physical constants
constexpr double epsilonLXe = 1; // vacuum permittivity (F/m)
constexpr double epsilonGXe = 1; // vacuum permittivity (F/m)
constexpr double VAnode = 0;             // Voltage on top plate (V)
constexpr double VCathode = -1000.0;             // Voltage on top plate (V)

// Geometry / mesh parameters (optional defaults)
constexpr double DriftRegionHeight      = 1;
constexpr int n_wires                   = 10;
constexpr double wire_diameter          = 0.025;
constexpr double inner_radius           = 0.2;
constexpr double outer_radius           = 0.25;
constexpr double anode_wire_height      = 0.95;
constexpr double cathode_wire_height    = 0.05;
constexpr double ring_thickness         = 0.05;

// Boundary Condition indeces
constexpr int anode_BC_index            = 1000;
constexpr int cathode_BC_index          = 1001;
constexpr int TPC_Volume_index          = 2001;



// TODO Is this used?
constexpr int nx = 200;                  // mesh divisions in x
constexpr int ny = 200;                 // mesh divisions in y

// FEM parameters
constexpr int fe_order = 3;             // linear finite elements

// Solver and tolerance parameters
constexpr double tol = 1e-10;
constexpr int maxiter = 100000;
constexpr int printlevel = 1;


// DEBUG
constexpr bool debug = true; 


#endif
