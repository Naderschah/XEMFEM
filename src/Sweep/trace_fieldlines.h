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


struct ElectronTraceResult
{
    std::vector<mfem::Vector> points; // (r,z) trajectory points
    ElectronExitCode          exit_code    = ElectronExitCode::None;
    int                       exit_element = -1;
};


struct TpcGeometry
{
    double r_min;
    double r_max;
    double z_min;
    double z_max;

    explicit TpcGeometry(const ElectronTraceParams &tp)
    {
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


// Electric field coefficient wrapper around a vector GridFunction
struct ElectricFieldCoeff : public mfem::VectorCoefficient
{
    mfem::GridFunction &phi;  // H1 potential
    double sign;              // +1.0 for grad(phi), -1.0 for -grad(phi)

    ElectricFieldCoeff(mfem::GridFunction &phi_, double sign_ = -1.0)
        : mfem::VectorCoefficient(phi_.FESpace()->GetMesh()->SpaceDimension())
        , phi(phi_)
        , sign(sign_)
    { }

    void Eval(mfem::Vector &V,
              mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override
    {
        T.SetIntPoint(&ip);
        phi.GetGradient(T, V);
        V *= sign; // gives physical E = -grad(phi)
    }
};

struct InterfacePoint
{
    double r;
    double z_ci;        // highest CI height at this radius
    int    elem;
    bool   cathode_ci;  // seen cathode-insensitive exits in this column
    bool   wall_ci;     // seen wall-insensitive exits in this column
};
// Shared seed container used by both CIV and tracing.
struct CivSeeds
{
    std::vector<mfem::Vector>          positions; // (r,z) seeds (physical coords)
    std::vector<double>                volumes;   // volume weight per seed
    std::vector<int>                   elements;  // element index for each seed
    std::vector<mfem::IntegrationPoint> ips;      // reference IP for each seed
};


// Mesh connectivity / size info used by the tracer
struct ElementAdjacency
{
    std::vector<std::vector<int>> neighbors; // Neighbor idx
    std::vector<double>           h;         // Characteristic Height
};
// Build Adjacency
ElementAdjacency BuildAdjacency(mfem::ParMesh &mesh);
// Find mesh element given coordinate using neighbor search
bool FindElementForPointLocal(
    mfem::ParMesh            &mesh,
    const ElementAdjacency   &adj,
    int                       current_elem,
    const mfem::Vector       &x_new,
    int                      &out_elem,
    mfem::IntegrationPoint   &out_ip);

// Will hold the stepping method
using StepFunction = void (*)(mfem::ParMesh&,
                              ElectricFieldCoeff&,
                              const ElementAdjacency&,
                              const TpcGeometry&,
                              const ElectronTraceParams&,
                              bool,
                              int,
                              const mfem::IntegrationPoint&,
                              double,
                              mfem::Vector&,
                              int&,
                              mfem::IntegrationPoint&,
                              double&,
                              bool&,
                              bool&,
                              ElectronExitCode&);

// Trace  A single electron through the TPC
ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh                &mesh,
    ElectricFieldCoeff           &E_coeff,
    const ElementAdjacency       &adj,
    const TpcGeometry            &geom,
    const ElectronTraceParams    &params,
    int                           start_elem,
    const mfem::IntegrationPoint &start_ip,
    bool                          axisymmetric,
    bool                          save_pathlines);

// Trace all supplied seeds and fill out_results with one ElectronTraceResult
// per seed. Tracing configuration and geometry come from cfg.tracing_params.
// The Inner function does the work, the outer does object creation
void TraceElectronFieldLinesInner(mfem::ParMesh                   &global_mesh,
                                  mfem::GridFunction              &global_phi,
                                  const ElementAdjacency          &adj,
                                  const mfem::FiniteElementCollection *fec_phi,
                                  int                              ordering_phi,
                                  const ElectronTraceParams        &params,
                                  const Config                     &cfg,
                                  const CivSeeds                   &seeds,
                                  std::vector<ElectronTraceResult> &out_results,
                                  bool                             axisymmetric,
                                  bool                             save_paths,
                                  const double                     *z_max_overrides);
void TraceElectronFieldLines(const SimulationResult           &sim,
                             const Config                     &cfg,
                             const CivSeeds                   &seeds,
                             std::vector<ElectronTraceResult> &out_results);