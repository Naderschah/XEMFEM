/*
We have 2 options
- Interpolate a grid
- Integrate on the mesh <- Doing this 





Little GPT summary on the error propagation here

Field-Line Tracing Error Model
------------------------------

We trace trajectories of the discrete electric field E_h defined on the FEM mesh.
The exact field is E, the discrete field is E_h, and the numerical trajectory uses
a time step Δt.

1. Field Approximation Error
   Let x(t) solve ẋ = −E(x), and let x_h(t) solve ẋ = −E_h(x).
   If  ||E − E_h||_{L∞} ≤ C h^m  and E_h is Lipschitz with constant L, then

      ||x_h(t) − x(t)|| ≤ C_T ||E − E_h||_{L∞} ≤ C_T' h^m ,

   so field lines of E_h converge to those of E as the mesh is refined.

2. Time-Discretization Error
   Let x̃_h^n be the numerical solution of order p with step Δt. Then

      ||x̃_h(t_n) − x_h(t_n)|| ≤ C_T Δt^p ,

   giving total trajectory error

      ||x̃_h(t_n) − x(t_n)|| ≤ C_1 h^m + C_2 Δt^p .

3. Step-Size Condition
   We enforce

      Δt ≤ c h_K / ||v_h(x)|| ,

   where h_K is the size of the current element. Then

      ||x_{n+1} − x_n|| ≤ c h_K ,

   which guarantees x_{n+1} lies in the current element K or one of its
   immediate neighbors (shape-regularity). This ensures neighbor-only
   element lookup is exact and introduces no additional modeling error.

4. Comparison with Grid Interpolation
   If the field is additionally interpolated onto a structured grid
   (spacing H, interpolation order q), producing E_g, then

      ||x_g(t) − x(t)|| ≤ C_T ( h^m + H^q ),

   so grid-based tracing adds an extra spatial error term H^q that the
   direct FEM-based tracer does not incur.




*/
#pragma once

#include <mfem.hpp>
#include <vector>
// TODO Add SolverCore to cmake 
#include "solver_api.h"
#include "Config.h"

enum class ElectronExitCode
{
    None = 0,
    HitCathode,
    HitLiquidGas,
    HitAxis,
    HitWall,
    LeftVolume,
    WeakField,
    MaxSteps,
    InvalidSeed   // currently not produced by the tracer; kept for completeness
};

struct ElectronTraceResult
{
    std::vector<mfem::Vector> points; // (r,z) trajectory points
    ElectronExitCode          exit_code    = ElectronExitCode::None;
    int                       exit_element = -1;
};

// Shared seed container used by both CIV and tracing.
struct CivSeeds
{
    std::vector<mfem::Vector>          positions; // (r,z) seeds (physical coords)
    std::vector<double>                volumes;   // volume weight per seed
    std::vector<int>                   elements;  // element index for each seed
    std::vector<mfem::IntegrationPoint> ips;      // reference IP for each seed
};

// ============================================================================
// Public API
// ============================================================================

// Build volume-weighted seeds for CIV / optimization / tracing
// based on the mesh and cfg.tracing_params.
CivSeeds ExtractCivSeeds(const Config &cfg,
                         const SimulationResult &result);

// Trace all supplied seeds and fill out_results with one ElectronTraceResult
// per seed. Tracing configuration and geometry come from cfg.tracing_params.
void TraceElectronFieldLines(const SimulationResult           &sim,
                             const Config                     &cfg,
                             const CivSeeds                   &seeds,
                             std::vector<ElectronTraceResult> &out_results);