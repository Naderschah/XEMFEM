#ifndef CONSTANTS_H
#define CONSTANTS_H

// Paths
constexpr const char* msh_file = "geometry.msh";

// Physical constants
constexpr double epsilon0 = 1; // vacuum permittivity (F/m)
constexpr double V0 = 1000.0;             // Voltage on top plate (V)

// Geometry / mesh parameters (optional defaults)
constexpr double plate_length = 0.5;   // meters
constexpr double plate_gap = 0.1;      
constexpr double plate_thickness = 0.1;
constexpr double thickness_z = 0.05;

constexpr double padding_length = 1;    
constexpr double padding_gap = 1;    
constexpr double W = plate_length + 2*padding_length;    
constexpr double H = plate_gap + 2*padding_gap;    
constexpr int nx = 200;                  // mesh divisions in x
constexpr int ny = 200;                 // mesh divisions in y

// FEM parameters
constexpr int fe_order = 3;             // linear finite elements


// Solver and tolerance parameters
constexpr double tol = 1e-5;
constexpr int maxiter = 100000;
constexpr int printlevel = 1;


// DEBUG
constexpr bool debug = true; 


#endif
