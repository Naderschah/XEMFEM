#pragma once
#include <mfem.hpp>
#include <filesystem>
#include <vector>
// TODO Add SolverCore to cmake 
#include "solver_api.h"
#include "Config.h"

// Seed container used by both CIV and tracing.
struct Seeds
{
    std::vector<mfem::Vector>          positions; // (r,z) seeds (physical coords)
    std::vector<double>                volumes;   // volume weight per seed
};

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
        if (z >  z_max) { return ElectronExitCode::HitLiquidGas; }
        if (z <  z_min) { return ElectronExitCode::HitCathode;   }

        // And radial boundaries 
        if (r >  r_max) { return ElectronExitCode::HitWall;      }

        return ElectronExitCode::None; // still inside cylinder
    }
};


inline long long ComputeMaxStepsFromLimit(const ElectronTraceParams &p,
                                                 double ds,
                                                 int max_traversals)
{
    if (max_traversals <= 0) return 0; // 0 means unlimited
    const double height = std::fabs(p.z_max - p.z_min);
    if (!(height > 0.0) || !(ds > 0.0)) return 0;

    const double target = static_cast<double>(max_traversals) * height;
    const double n = target / ds;
    const long long max_steps = static_cast<long long>(std::ceil(n));
    return (max_steps > 0) ? max_steps : 1;
}

inline double ComputeGlobalMinEdgeLength(mfem::ParMesh &mesh, bool debug)
{
    using namespace mfem;
    double local_min = std::numeric_limits<double>::infinity();

    for (int e = 0; e < mesh.GetNE(); ++e)
    {
        const Element *el = mesh.GetElement(e);
        const int nv = el->GetNVertices();
        const int *v = el->GetVertices();

        for (int i = 0; i < nv; ++i)
        {
            const double *pi = mesh.GetVertex(v[i]);
            for (int j = i + 1; j < nv; ++j)
            {
                const double *pj = mesh.GetVertex(v[j]);
                const double dx = pi[0] - pj[0];
                const double dy = pi[1] - pj[1];
                const double dz = (mesh.SpaceDimension() > 2) ? (pi[2] - pj[2]) : 0.0;
                const double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d > 0.0 && d < local_min) local_min = d;
            }
        }
    }

    double global_min = local_min;
#ifdef MFEM_USE_MPI
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, mesh.GetComm());
#endif
    if (!std::isfinite(global_min) || global_min <= 0.0) global_min = 1.0;
    return global_min;
}

void DumpElectronPathsCSV(const Config                          &cfg,
                          const std::vector<ElectronTraceResult> &out_results);

                          