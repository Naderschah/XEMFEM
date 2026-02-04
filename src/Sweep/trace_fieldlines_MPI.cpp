#include "trace_fieldlines_MPI.h"

namespace mpitracing
{
static inline void ClampAxisIfNeeded(mfem::Vector &x, double geom_tol, bool axisymmetric)
{
    if (!axisymmetric) { return; }
    if (x.Size() >= 1 && x[0] <= geom_tol) { x[0] = geom_tol; }
}

static inline ElectronExitCode ClassifyExitOrLeftVolume(const TpcGeometry &geom,
                                                        const mfem::Vector &x,
                                                        double geom_tol)
{
    ElectronExitCode c = geom.ClassifyBoundary(x[0], x[1], geom_tol);
    if (c == ElectronExitCode::None) { c = ElectronExitCode::LeftVolume; }
    return c;
}

static inline void AddFinished(std::vector<TraceSummaryPOD> &finished,
                               int seed_id,
                               ElectronExitCode code,
                               const mfem::Vector &x,
                               int sdim)
{
    TraceSummaryPOD s{};
    s.seed_id   = (int32_t)seed_id;
    s.exit_code = (int32_t)code;
    for (int d = 0; d < 3; ++d) { s.x[d] = 0.0; }
    for (int d = 0; d < sdim; ++d) { s.x[d] = x[d]; }
    finished.push_back(s);
}

static void GatherFinishedToRoot(MPI_Comm comm,
                                 int rank,
                                 int size,
                                 const std::vector<TraceSummaryPOD> &finished_local,
                                 std::size_t n_seeds,
                                 int sdim,
                                 std::vector<ElectronTraceResult> &out_results)
{
    const int sendcount = (int)finished_local.size();

    std::vector<int> recvcounts;
    if (rank == 0) { recvcounts.resize((std::size_t)size); }

    MPI_Gather(&sendcount, 1, MPI_INT,
               rank == 0 ? recvcounts.data() : nullptr, 1, MPI_INT,
               0, comm);

    std::vector<int> recvcountsB, displsB;
    std::vector<char> all_bytes;

    if (rank == 0)
    {
        recvcountsB.resize((std::size_t)size);
        displsB.resize((std::size_t)size);

        int offB = 0;
        for (int r = 0; r < size; ++r)
        {
            recvcountsB[(std::size_t)r] = recvcounts[(std::size_t)r] * (int)sizeof(TraceSummaryPOD);
            displsB[(std::size_t)r] = offB;
            offB += recvcountsB[(std::size_t)r];
        }
        all_bytes.resize((std::size_t)offB);
    }

    MPI_Gatherv((void*)finished_local.data(),
                sendcount * (int)sizeof(TraceSummaryPOD), MPI_BYTE,
                rank == 0 ? all_bytes.data() : nullptr,
                rank == 0 ? recvcountsB.data() : nullptr,
                rank == 0 ? displsB.data() : nullptr,
                MPI_BYTE, 0, comm);

    if (rank != 0)
    {
        out_results.clear();
        return;
    }

    // Reconstruct results on root
    out_results.assign(n_seeds, ElectronTraceResult{});

    const std::size_t n_all = all_bytes.size() / sizeof(TraceSummaryPOD);
    const auto *S = reinterpret_cast<const TraceSummaryPOD*>(all_bytes.data());

    // Fill from summaries
    std::vector<uint8_t> seen(n_seeds, 0);

    for (std::size_t i = 0; i < n_all; ++i)
    {
        const int sid = (int)S[i].seed_id;
        if (sid < 0 || (std::size_t)sid >= n_seeds) { continue; }

        out_results[(std::size_t)sid].exit_code = (ElectronExitCode)S[i].exit_code;
        out_results[(std::size_t)sid].points.clear();

        mfem::Vector v(sdim);
        for (int d = 0; d < sdim; ++d) { v[d] = S[i].x[d]; }
        out_results[(std::size_t)sid].points.push_back(std::move(v));

        seen[(std::size_t)sid] = 1;
    }

    // Safety: ensure all seeds produced a result
    for (std::size_t sid = 0; sid < n_seeds; ++sid)
    {
        if (!seen[sid])
        {
            // leave exit_code as None; caller may treat as error
            // but make it explicit here:
            out_results[sid].exit_code = ElectronExitCode::None;
            out_results[sid].points.clear();
        }
    }
}

void TraceDistributedEuler(
    mfem::ParMesh& mesh,
    mfem::FindPointsGSLIB& finder,
    const mfem::ParGridFunction& E_gf,
    const ElectronTraceParams& params,
    const Seeds& seeds,
    std::vector<ElectronTraceResult>& out_results,
    double h_ref,
    bool axisymmetric,
    int redistribution_every,
    bool debug,
    int debug_every)
{
    MPI_Comm comm = mesh.GetComm();
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int sdim = mesh.SpaceDimension();
    MFEM_VERIFY(E_gf.VectorDim() == sdim, "E_gf.VectorDim() must equal mesh.SpaceDimension().");

    const std::size_t n_seeds = seeds.positions.size();

    // Partition seeds so all ranks have initial particles
    std::size_t ib = 0, ie = 0;
    SeedRangeForRank(n_seeds, rank, size, ib, ie);
    const std::size_t n_local_init = ie - ib;

    const double ds = params.c_step * h_ref;
    MFEM_VERIFY(ds > 0.0 && std::isfinite(ds), "TraceDistributedEuler: ds must be positive finite.");

    const long long max_steps = ComputeMaxStepsFromLimit(params, ds, params.max_traversals);
    const long long N = (max_steps > 0) ? max_steps : std::numeric_limits<long long>::max();

    const TpcGeometry geom(params);

    // Active particle container (miniapp style)
    mfem::ParticleSet active(comm, (int)n_local_init, sdim, mfem::Ordering::byNODES);

    // Tags
    const int tag_seed = active.AddTag("seed_id");

    // Fields: E and direction, same ordering as coords (usually byNODES)
    const mfem::Ordering::Type ord = active.Coords().GetOrdering();
    const int fld_E   = active.AddField(sdim, ord, "E");
    const int fld_dir = active.AddField(sdim, ord, "dir");

    // Initialize coordinates and seed ids
    mfem::ParticleVector &X = active.Coords();
    if (debug && rank == 0)
    {
        std::cout << "[MPITracer:init] Coords.vdim=" << X.GetVDim()
                  << " Coords.Size=" << X.Size()
                  << " NParticles=" << active.GetNParticles()
                  << std::endl;
    }
    MFEM_VERIFY(X.Size() == sdim * (int)n_local_init,
            "Unexpected ParticleSet coordinate size (expected sdim*n_local_init).");
    MFEM_VERIFY(active.GetNParticles() == (int)n_local_init,
                "Unexpected ParticleSet particle count.");
    MFEM_VERIFY(X.GetVDim() == sdim, "Coords vdim mismatch.");
    MFEM_VERIFY(X.GetOrdering() == mfem::Ordering::byNODES,
                "Coords ordering must be byNODES for this tracer.");
    if (debug && rank == 0)
    {
        std::cout << "[MPITracer:init] sdim=" << sdim
                  << " mesh.SpaceDimension()=" << mesh.SpaceDimension()
                  << " E_gf.VectorDim()=" << E_gf.VectorDim()
                  << " n_local_init=" << n_local_init
                  << std::endl;
    }

    for (std::size_t li = 0; li < n_local_init; ++li)
    {
        const std::size_t gi = ib + li;
        mfem::Vector xi(sdim);
        for (int d = 0; d < sdim; ++d) { xi[d] = seeds.positions[gi][d]; }
        ClampAxisIfNeeded(xi, params.geom_tol, axisymmetric);
        for (int d = 0; d < sdim; ++d) { 
            X((int)li, d) = xi[d]; 
        }
        active.Tag(tag_seed)[(int)li] = (int)gi;
    }

    std::vector<TraceSummaryPOD> finished_local;
    finished_local.reserve(n_local_init);

    auto kill_indices =
        [&](const mfem::Array<int> &rm, const std::vector<ElectronExitCode> *codes_opt)
    {
        if (rm.Size() == 0) { return; }

        mfem::ParticleVector &Xk = active.Coords();

        for (int j = 0; j < rm.Size(); ++j)
        {
            const int idx = rm[j];
            const int sid = active.Tag(tag_seed)[idx];

            mfem::Vector xi(sdim);
            for (int d = 0; d < sdim; ++d) { xi[d] = Xk(idx, d); }

            ElectronExitCode code = ElectronExitCode::LeftVolume;
            if (codes_opt) { code = (*codes_opt)[(std::size_t)j]; }

            AddFinished(finished_local, sid, code, xi, sdim);
        }

        active.RemoveParticles(rm);
    };

    // Main stepping loop
    for (long long step = 0; step < N; ++step)
    {
        const int n_active = active.GetNParticles();

        int global_active = 0;
        MPI_Allreduce(&n_active, &global_active, 1, MPI_INT, MPI_SUM, comm);
        if (global_active == 0) { break; }

        if (debug && debug_every > 0 && (step % debug_every == 0))
        {
            std::cout << "[rank " << rank << "] step=" << step
                      << " local_active=" << n_active
                      << " global_active=" << global_active
                      << std::endl;
        }

        
        // Locate current points (even if none are available for MPI Interaction)
        finder.FindPoints(active.Coords(), active.Coords().GetOrdering());

        const mfem::Array<int>       &codes0 = finder.GetCode();
        const mfem::Array<long long> &elem0  = finder.GetElem();

        // Interpolate E into particle field
        mfem::ParticleVector &Evals = active.Field(fld_E);
        finder.Interpolate(E_gf, Evals);

        mfem::ParticleVector &dir = active.Field(fld_dir);
        if (n_active > 0)
        {
          mfem::Array<int> rm_now; rm_now.Reserve(n_active);
          mfem::Array<int> rm_deg; rm_deg.Reserve(n_active);

          for (int i = 0; i < n_active; ++i)
          {
              const bool ok = (elem0[i] >= 0 && codes0[i] == 0);
              if (!ok) { rm_now.Append(i); continue; }

              double n2 = 0.0;
              for (int d = 0; d < sdim; ++d)
              {
                  const double e = Evals(i, d);        // <-- FIX
                  n2 += e * e;
              }
              const double nrm = std::sqrt(n2);
              if (!(nrm > 0.0) || !std::isfinite(nrm)) { rm_deg.Append(i); continue; }

              for (int d = 0; d < sdim; ++d)
              {
                  dir(i, d) = -Evals(i, d) / nrm;      // <-- FIX
              }
          }

          if (rm_now.Size() > 0) { kill_indices(rm_now, nullptr); }

          if (rm_deg.Size() > 0)
          {
              std::vector<ElectronExitCode> codes;
              codes.resize((std::size_t)rm_deg.Size(), ElectronExitCode::DegenerateTimeStep);
              kill_indices(rm_deg, &codes);
          }

          const int n_active2 = active.GetNParticles();
          if (n_active2 > 0)
          {
              finder.FindPoints(active.Coords(), active.Coords().GetOrdering());

              const mfem::Array<int>       &codes1 = finder.GetCode();
              const mfem::Array<long long> &elem1  = finder.GetElem();

              mfem::ParticleVector &E2 = active.Field(fld_E);
              finder.Interpolate(E_gf, E2);

              mfem::ParticleVector &dir2 = active.Field(fld_dir);

              mfem::Array<int> rm_bad; rm_bad.Reserve(n_active2);

              for (int i = 0; i < n_active2; ++i)
              {
                  const bool ok = (elem1[i] >= 0 && codes1[i] == 0);
                  if (!ok) { rm_bad.Append(i); continue; }

                  double n2 = 0.0;
                  for (int d = 0; d < sdim; ++d)
                  {
                      const double e = E2(i, d);        // <-- FIX
                      n2 += e * e;
                  }
                  const double nrm = std::sqrt(n2);
                  if (!(nrm > 0.0) || !std::isfinite(nrm)) { rm_bad.Append(i); continue; }

                  for (int d = 0; d < sdim; ++d)
                  {
                      dir2(i, d) = -E2(i, d) / nrm;     // <-- FIX
                  }
              }

              if (rm_bad.Size() > 0) { kill_indices(rm_bad, nullptr); }
          }

          const int n_active3 = active.GetNParticles();
          if (n_active3 > 0)
          {
              mfem::ParticleVector &X3   = active.Coords();
              mfem::ParticleVector &dir3 = active.Field(fld_dir);

              // Euler update
              for (int i = 0; i < n_active3; ++i)
              {
                  mfem::Vector xi(sdim);
                  for (int d = 0; d < sdim; ++d) { xi[d] = X3(i, d); } 

                  for (int d = 0; d < sdim; ++d) { xi[d] += ds * dir3(i, d); }

                  ClampAxisIfNeeded(xi, params.geom_tol, axisymmetric);

                  for (int d = 0; d < sdim; ++d) { X3(i, d) = xi[d]; }
              }

              mfem::Array<int> rm_geom; rm_geom.Reserve(n_active3);
              std::vector<ElectronExitCode> rm_geom_codes;
              rm_geom_codes.reserve((std::size_t)n_active3);

              for (int i = 0; i < n_active3; ++i)
              {
                  const double r = X3(i, 0);
                  const double z = X3(i, 1);

                  if (!geom.Inside(r, z))
                  {
                      rm_geom.Append(i);
                      mfem::Vector xi(sdim);
                      for (int d = 0; d < sdim; ++d) { xi[d] = X3(i, d); } 
                      rm_geom_codes.push_back(ClassifyExitOrLeftVolume(geom, xi, params.geom_tol));
                  }
              }

              if (rm_geom.Size() > 0) { kill_indices(rm_geom, &rm_geom_codes); }
          }
        }
        // Redistribute periodically (collective)
        if (redistribution_every > 0 && ((step % redistribution_every) == 0))
        {
            // Decide collectively whether there is anything to redistribute
            int n_local  = active.GetNParticles();
            int n_global = 0;
            MPI_Allreduce(&n_local, &n_global, 1, MPI_INT, MPI_SUM, comm);

            if (n_global > 0)
            {
                // ---- 1) Remove not-found particles (local operation) ----
                finder.FindPoints(active.Coords(), active.Coords().GetOrdering());
                const auto &codesR = finder.GetCode();
                const auto &elemR  = finder.GetElem();

                if (n_local > 0)
                {
                    mfem::Array<int> rm_nf;
                    rm_nf.Reserve(n_local);

                    for (int i = 0; i < n_local; ++i)
                    {
                        if (elemR[i] < 0 || codesR[i] != 0)
                        {
                            rm_nf.Append(i);
                        }
                    }

                    if (rm_nf.Size() > 0)
                    {
                        kill_indices(rm_nf, nullptr);
                    }
                }

                // ---- 2) Recompute owners after removals ----
                finder.FindPoints(active.Coords(), active.Coords().GetOrdering());
                const auto &procs = finder.GetProc();

                mfem::Array<unsigned int> procs2(active.GetNParticles());
                for (int i = 0; i < active.GetNParticles(); ++i)
                {
                    procs2[i] = procs[i];
                }

                // ---- 3) Collective redistribute (ALL ranks call this) ----
                active.Redistribute(procs2);
            }
        }
    }
    // Mark remaining active particles as MaxSteps
    if (active.GetNParticles() > 0)
    {
        mfem::ParticleVector &Xf = active.Coords();
        for (int i = 0; i < active.GetNParticles(); ++i)
        {
            const int sid = active.Tag(tag_seed)[i];
            mfem::Vector xi(sdim);
            for (int d = 0; d < sdim; ++d) { xi[d] = Xf(i, d); }
            AddFinished(finished_local, sid, ElectronExitCode::MaxSteps, xi, sdim);
        }

        mfem::Array<int> rm_all(active.GetNParticles());
        for (int i = 0; i < active.GetNParticles(); ++i) { rm_all[i] = i; }
        active.RemoveParticles(rm_all);
    }

    // Gather finished to root and construct out_results
    GatherFinishedToRoot(comm, rank, size, finished_local, n_seeds, sdim, out_results);
}

}

