#pragma once
#include <mfem.hpp>
#include <filesystem>
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

void DumpElectronPathsCSV(const Config                          &cfg,
                          const std::vector<ElectronTraceResult> &out_results);
