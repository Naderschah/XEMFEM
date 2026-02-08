#pragma once

#include <mfem.hpp>
#include <mfem/fem/particleset.hpp>
#include <mpi.h>
#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <cstring>

#include "./common_tracing/tracing_objects.h"

namespace mpitracing
{
struct TraceSummaryPOD
{
    int32_t seed_id;    // global seed index (0..n_seeds-1)
    int32_t exit_code;  // ElectronExitCode as int32
    double  x[3];       // final position; x[2]=0 for 2D
};
static_assert(std::is_trivially_copyable<TraceSummaryPOD>::value,
              "TraceSummaryPOD must be POD for MPI gather.");

struct TrackPointPOD
{
    int32_t seed_id;
    int32_t step;
    double  x[3]; // use first sdim entries, pad with 0
};
static_assert(std::is_trivially_copyable<TrackPointPOD>::value, "POD required");


inline void SeedRangeForRank(std::size_t n, int rank, int size,
                            std::size_t &begin, std::size_t &end)
{
    const std::size_t base = n / (std::size_t)size;
    const std::size_t rem  = n % (std::size_t)size;

    begin = (std::size_t)rank * base + std::min<std::size_t>((std::size_t)rank, rem);
    end   = begin + base + ((std::size_t)rank < rem ? 1 : 0);
}

// Main driver: distributed Euler stepping + periodic redistribution.
// redistribution_every:
//   - 1 => redistribute every step
//   - k => redistribute every k steps
//   - 0 => never redistribute (debug only; generally incorrect for multi-rank)
void TraceDistributedEuler(
    mfem::ParMesh& mesh,
    mfem::FindPointsGSLIB& finder,           // per-rank instance, Setup(mesh) already done
    const mfem::ParGridFunction& E_gf,       // vector field, VectorDim == mesh.SpaceDimension()
    const ElectronTraceParams& params,       // uses c_step, geom_tol, bounds, max_traversals
    const Seeds& seeds,                      // global seeds (all ranks can see; each rank will take a slice)
    std::vector<ElectronTraceResult>& out_results, // rank 0 only filled
    double h_ref,
    bool axisymmetric,
    int redistribution_every,
    bool debug,
    int debug_every = 5,
    bool save_path = false);
} // namespace mpitracing

struct MPITraceContext
{
    std::unique_ptr<mfem::FindPointsGSLIB> finder;
    double h_ref = 1.0;

    void Build(mfem::ParMesh &mesh, bool debug);
    bool Ready() const { return (bool)finder; }
};