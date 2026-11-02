#pragma once

#include <mfem.hpp>
#include <mfem/fem/particleset.hpp>
#include <mpi.h>
#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <cstring>

#include "./common_tracing/tracing_objects.h"

namespace mpitracing
{
struct TraceProgressState;

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



/// Evaluates electric field E = -∇φ at arbitrary points on a distributed mesh.
/// Uses GSLIB for point location and MPI communication.
class FieldEvaluator {
public:
    /// Constructor.
    /// @param mesh   Distributed mesh (must be the same as used by finder).
    /// @param phi    Electrostatic potential grid function.
    /// @param finder GSLIB point locator already attached to the mesh.
    /// @param sdim   Spatial dimension (mesh.SpaceDimension()).
    FieldEvaluator(mfem::ParMesh& mesh,
                   const mfem::ParGridFunction& phi,
                   mfem::FindPointsGSLIB& finder,
                   int sdim,
                   TraceProgressState *progress = nullptr);

    /// Evaluate the electric field at `n_points` points.
    /// @param[in]  points   Coordinates, size = n_points * sdim, ordering = byVDIM.
    /// @param[in]  n_points Number of points.
    /// @param[out] E        Electric field, size = n_points * sdim, ordering = byVDIM.
    ///                      For points not found, E entries are set to NaN.
    /// Collective on the MPI communicator of the mesh.
    void EvaluateE(const mfem::Vector& points,
                   int n_points,
                   mfem::Vector& E,
                   mfem::Vector *element_sizes = nullptr);

private:
    mfem::ParMesh& mesh_;
    const mfem::ParGridFunction& field_;
    mfem::FindPointsGSLIB& finder_;
    const int sdim_;
    TraceProgressState *progress_;

    // Temporary buffers reused across calls to avoid reallocation.
    mfem::Array<unsigned int> recv_elem_;
    mfem::Vector recv_ref_;      // reference coordinates, size = n_recv * sdim
    mfem::Array<unsigned int> recv_code_;
    mfem::Vector recv_h_;        // current-element characteristic sizes
};


/// Abstract base class for ODE integrators (Euler, Runge‑Kutta, etc.)
class Integrator {
public:
    virtual ~Integrator() = default;

    /// Optional: called once after ParticleSet creation.
    /// Allows the integrator to add required per‑particle fields
    /// (e.g., "k1", "k2", ... for RK4) to the ParticleSet.
    virtual void Initialize(mfem::ParticleSet& particles) { }

    /// Perform one time step for all active particles.
    ///
    /// @param[in,out] particles   Current particle set.
    ///   On entry:  particles.Coords() contains current positions (byVDIM order).
    ///   On exit:   particles.Coords() must contain the positions after the step.
    /// @param[in]  step_parameter Global step size in legacy mode, or c_step
    ///                            multiplier in local-element mode.
    /// @param[in]  field         Field evaluator – used to compute E at arbitrary points.
    /// @param[out] to_remove     Indices of particles that must be removed due to
    ///                           degenerate field / NaN during this step.
    /// @param[out] exit_codes    Corresponding exit codes for each index in to_remove.
    ///                           The two arrays have the same length.
    ///
    /// The integrator may add entries to to_remove/exit_codes in any order,
    /// but each index must be valid (0 ≤ index < particles.GetNParticles()).
    /// The tracer will remove those particles and record the exit code.
    ///
    /// This method is called collectively by all MPI ranks that have particles.
    virtual void Step(mfem::ParticleSet& particles,
                      double step_parameter,
                      bool local_element_step,
                      FieldEvaluator& field,
                      mfem::Array<int>& to_remove,
                      std::vector<ElectronExitCode>& exit_codes,
                      mfem::Vector& step_lengths) = 0;
};


// redistribution_every:
//   - 1 => redistribute every step
//   - k => redistribute every k steps
//   - 0 => never redistribute (incorrect for a partitioned multi-rank mesh,
//          but correct for a complete mesh on MPI_COMM_SELF)
void TraceDistributed(
    mfem::ParMesh& mesh,
    mfem::FindPointsGSLIB& finder,           // per-rank instance, Setup(mesh) already done
    const mfem::ParGridFunction& E_gf,       // projected electric field
    const ElectronTraceParams& params,       // uses c_step, geom_tol, bounds, max_traversals
    const Seeds& seeds,                      // global seeds (all ranks can see; each rank will take a slice)
    std::vector<ElectronTraceResult>& out_results, // rank 0 only filled
    double h_ref,
    bool axisymmetric,
    int redistribution_every,
    bool debug,
    int debug_every = 5,
    bool save_path = false,
    double progress_interval_seconds = 0.0);

/// Trace statically strided seed jobs using a complete private mesh on every
/// world rank. The integration itself runs only on MPI_COMM_SELF; MPI world is
/// used once afterwards to gather results in global seed order.
void TraceReplicated(
    MPI_Comm world_comm,
    mfem::ParMesh& replicated_mesh,
    mfem::FindPointsGSLIB& replicated_finder,
    const mfem::ParGridFunction& replicated_phi, // projected electric field
    const ElectronTraceParams& params,
    const Seeds& seeds,
    std::vector<ElectronTraceResult>& out_results,
    double h_ref,
    bool axisymmetric,
    int redistribution_every,
    bool print_progress,
    bool save_path = false);
} // namespace mpitracing

struct MPITraceContext
{
    std::unique_ptr<mfem::FindPointsGSLIB> finder;
    double h_ref = 1.0;

    void Build(mfem::ParMesh &mesh, bool debug);
    bool Ready() const { return (bool)finder; }
};

/// A complete tracing mesh and projected field owned privately by one world
/// rank. Declaration order is intentional: dependent objects are destroyed
/// before their spaces and mesh.
struct ReplicatedMPITraceContext
{
    std::unique_ptr<mfem::ParMesh> replicated_mesh;
    std::unique_ptr<mfem::ParGridFunction> replicated_field;
    std::unique_ptr<mfem::FindPointsGSLIB> finder;
    double h_ref = 1.0;
    std::size_t serialized_mesh_bytes = 0;
    std::size_t serialized_field_bytes = 0;

    void Build(mfem::ParMesh &source_mesh,
               const mfem::ParGridFunction &source_field,
               bool debug);
    bool Ready() const
    {
        return replicated_mesh && replicated_field && finder;
    }
};
