/*
Here we integrate field line trajectories of the electrons

We use RK4 with the physical equation expected rather than topological
to allow later extension if we choose to include other forces

dr/dt = v(r) = - mu E(r)

With mu the effective mobility 

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
    HitWall,
    LeftVolume,
    MaxSteps,
    DegenerateTimeStep,
    HitAxis
};

// Mesh connectivity / size info used by the tracer
struct ElementAdjacency
{
    std::vector<std::vector<int>> neighbors; // Neighbor idx
    std::vector<double>           h;         // Characteristic Height
};


struct DriftTimeEstimate
{
    double H            = 0.0;  // drift height
    double E_rms        = 0.0;  // RMS field in tracking region
    double a_typ        = 0.0;  // |q/m| * E_rms
    double t_drift      = 0.0;  // time to cross H under a_typ
    double dt_traverse  = 0.0;  // base dt giving T_trav drift crossings in max_steps
    double crossings_est= 0.0;  // estimated crossings with dt_traverse (= T_trav)
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
struct TpcGeometry
{
    double r_min;
    double r_max;
    double z_min;
    double z_max;

    explicit TpcGeometry(const Config &cfg)
    {
        const auto &tp = cfg.tracing_params;
        r_min = tp.r_min;
        r_max = tp.r_max;
        z_min = tp.z_min;
        z_max = tp.z_max;
    }

    bool Inside(double r, double z) const
    {
        return (r >= r_min && r <= r_max &&
                z >= z_min && z <= z_max);
    }

    ElectronExitCode ClassifyBoundary(double r, double z, double geom_tol) const
    // After lots of annoyance this became a hard cut 
    {
        // Assume did not hit the wall when leaving z bounds so could 
        // still hit the wall at the top/bottom most portion but we save two checks
        if (z >  z_max + geom_tol) { return ElectronExitCode::HitLiquidGas; }
        if (z <  z_min - geom_tol) { return ElectronExitCode::HitCathode;   }

        // And radial boundaries 
        if (r >  r_max) { return ElectronExitCode::HitWall;      }

        return ElectronExitCode::None; // still inside cylinder
    }
};

// Holds data for the binary CIV search
struct ElementGeometryData
{
    std::vector<double> r_centroid;
    std::vector<double> z_centroid;
    std::vector<double> volume;
};
struct ColumnRefineResult
{
    bool   found;
    double z_ci;              // refined highest CI height
    bool   is_cathode_ci;
    bool   is_wall_ci;
};
struct InterfacePoint
{
    double r;
    double z_ci;        // highest CI height at this radius
    int    elem;
    bool   cathode_ci;  // seen cathode-insensitive exits in this column
    bool   wall_ci;     // seen wall-insensitive exits in this column
};
struct CIVolumeResult
{
    double V_total; // total volume in the considered region
    double V_CI;    // charge-insensitive volume (cathode + wall)
};
struct InitialColumnCIResult
{
    bool   found;
    double r;
    double z_ci;          // highest CI height at r
    int    elem;          // element containing (r,z_ci)
    bool   cathode_ci;    // seen cathode CI
    bool   wall_ci;       // seen wall CI
};
struct ColumnCISearchResult
{
    bool   found;             // whether any CI was detected
    double z_ci;              // highest z with CI along this column
    bool   is_cathode_ci;     // true if CI dominated by cathode at this column
    bool   is_wall_ci;        // true if CI dominated by wall at this column
};

struct ProbeRequest
{
    double r;
    double z;
    int    elem_guess;   // optional (-1 means unknown)
};
struct ProbeResult
{
    bool   ok         = false; // mapping + tracing succeeded
    bool   is_ci      = false; // IsChargeInsensitive(...)
    bool   is_cathode = false; // IsCathodeInsensitive(...)
    bool   is_wall    = false; // IsWallInsensitive(...)

    double r          = 0.0;   // probed (r,z)
    double z          = 0.0;
    int    elem       = -1;    // element id returned/used by tracer
    int    exit_code  = -1;    // int(ElectronExitCode) for debugging
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

static double ComputeCIV_ColumnSweep(const SimulationResult &sim,
                                     const Config           &cfg);

double compute_civ(const Config &cfg, const SimulationResult &result);