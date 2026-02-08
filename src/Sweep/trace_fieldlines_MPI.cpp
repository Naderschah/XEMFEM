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
                                 const std::vector<TrackPointPOD>   *track_log_local, // nullable
                                 std::size_t n_seeds,
                                 int sdim,
                                 std::vector<ElectronTraceResult> &out_results)
{
    // Decide tracks collectively (IMPORTANT)
    int do_tracks_local = (track_log_local != nullptr) ? 1 : 0;
    int do_tracks = 0;
    MPI_Allreduce(&do_tracks_local, &do_tracks, 1, MPI_INT, MPI_MAX, comm);
    // If do_tracks==1, then every rank should participate in track gather.
    // Ranks with no local tracks just send 0.

    // -----------------------------
    // 1) Gather finished summaries
    // -----------------------------
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

    // MPI requires a valid pointer even for 0 counts on some builds; keep it safe
    const void *sendbuf_sum = finished_local.empty() ? (const void*)nullptr : (const void*)finished_local.data();

    MPI_Gatherv((void*)sendbuf_sum,
                sendcount * (int)sizeof(TraceSummaryPOD), MPI_BYTE,
                rank == 0 ? all_bytes.data() : nullptr,
                rank == 0 ? recvcountsB.data() : nullptr,
                rank == 0 ? displsB.data() : nullptr,
                MPI_BYTE, 0, comm);

    // Root: reconstruct results (last point always)
    if (rank == 0)
    {
        out_results.assign(n_seeds, ElectronTraceResult{});

        const std::size_t n_all = all_bytes.size() / sizeof(TraceSummaryPOD);
        const auto *S = reinterpret_cast<const TraceSummaryPOD*>(all_bytes.data());

        std::vector<uint8_t> seen(n_seeds, 0);

        for (std::size_t i = 0; i < n_all; ++i)
        {
            const int sid = (int)S[i].seed_id;
            if (sid < 0 || (std::size_t)sid >= n_seeds) { continue; }

            auto &res = out_results[(std::size_t)sid];
            res.exit_code = (ElectronExitCode)S[i].exit_code;
            res.points.clear();

            mfem::Vector v(sdim);
            for (int d = 0; d < sdim; ++d) { v[d] = S[i].x[d]; }
            res.points.push_back(std::move(v));

            seen[(std::size_t)sid] = 1;
        }

        for (std::size_t sid = 0; sid < n_seeds; ++sid)
        {
            if (!seen[sid])
            {
                out_results[sid].exit_code = ElectronExitCode::None;
                out_results[sid].points.clear();
            }
        }
    }
    else
    {
        // keep it empty for now; we'll clear at end after track collectives too
        out_results.clear();
    }

    // -----------------------------------------
    // 2) Optional: Gather and merge track points
    // -----------------------------------------
    if (do_tracks)
    {
        const int send_n_track = (track_log_local ? (int)track_log_local->size() : 0);

        std::vector<int> recv_n_track;
        if (rank == 0) { recv_n_track.resize((std::size_t)size); }

        MPI_Gather((void*)&send_n_track, 1, MPI_INT,
                   rank == 0 ? recv_n_track.data() : nullptr, 1, MPI_INT,
                   0, comm);

        std::vector<int> recvcountsTB, displsTB;
        std::vector<char> all_track_bytes;

        if (rank == 0)
        {
            recvcountsTB.resize((std::size_t)size);
            displsTB.resize((std::size_t)size);

            int offB = 0;
            for (int r = 0; r < size; ++r)
            {
                recvcountsTB[(std::size_t)r] = recv_n_track[(std::size_t)r] * (int)sizeof(TrackPointPOD);
                displsTB[(std::size_t)r] = offB;
                offB += recvcountsTB[(std::size_t)r];
            }
            all_track_bytes.resize((std::size_t)offB);
        }

        const void *sendbuf_tr =
            (track_log_local && !track_log_local->empty()) ? (const void*)track_log_local->data() : (const void*)nullptr;

        MPI_Gatherv((void*)sendbuf_tr,
                    send_n_track * (int)sizeof(TrackPointPOD), MPI_BYTE,
                    rank == 0 ? all_track_bytes.data() : nullptr,
                    rank == 0 ? recvcountsTB.data() : nullptr,
                    rank == 0 ? displsTB.data() : nullptr,
                    MPI_BYTE, 0, comm);

        if (rank == 0)
        {
            const std::size_t nT = all_track_bytes.size() / sizeof(TrackPointPOD);
            const auto *P = reinterpret_cast<const TrackPointPOD*>(all_track_bytes.data());

            std::vector<std::vector<TrackPointPOD>> per_seed(n_seeds);

            for (std::size_t k = 0; k < nT; ++k)
            {
                const int sid = (int)P[k].seed_id;
                if (sid < 0 || (std::size_t)sid >= n_seeds) { continue; }
                per_seed[(std::size_t)sid].push_back(P[k]);
            }

            for (std::size_t sid = 0; sid < n_seeds; ++sid)
            {
                auto &vec = per_seed[sid];
                if (vec.empty()) { continue; } // keep summary last-point only

                std::sort(vec.begin(), vec.end(),
                          [](const TrackPointPOD &a, const TrackPointPOD &b)
                          {
                              return a.step < b.step;
                          });

                auto &pts = out_results[sid].points;
                pts.clear();
                pts.reserve(vec.size());

                for (const auto &tp : vec)
                {
                    mfem::Vector x(sdim);
                    for (int d = 0; d < sdim; ++d) { x[d] = tp.x[d]; }
                    pts.push_back(std::move(x));
                }
            }
        }
    }

    // Ensure non-root ends empty (by contract)
    if (rank != 0)
    {
        out_results.clear();
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
    int debug_every,
    bool save_path)
{
    MPI_Comm comm = mesh.GetComm();
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int sdim = mesh.SpaceDimension();
    MFEM_VERIFY(E_gf.VectorDim() == sdim, "E_gf.VectorDim() must equal mesh.SpaceDimension().");

    const std::size_t n_seeds = seeds.positions.size();

    std::size_t ib = 0, ie = 0;
    SeedRangeForRank(n_seeds, rank, size, ib, ie);
    const std::size_t n_local_init = ie - ib;

    const double ds = params.c_step * h_ref;
    MFEM_VERIFY(ds > 0.0 && std::isfinite(ds), "TraceDistributedEuler: ds must be positive finite.");

    const long long max_steps = ComputeMaxStepsFromLimit(params, ds, params.max_traversals);
    const long long N = (max_steps > 0) ? max_steps : std::numeric_limits<long long>::max();

    const TpcGeometry geom(params);

    mfem::ParticleSet active(comm, (int)n_local_init, sdim, mfem::Ordering::byVDIM);

    const int tag_seed = active.AddTag("seed_id");

    const mfem::Ordering::Type ord = active.Coords().GetOrdering();
    const int fld_E   = active.AddField(sdim, ord, "E");
    const int fld_dir = active.AddField(sdim, ord, "dir");

    mfem::ParticleVector &X = active.Coords();

    for (std::size_t li = 0; li < n_local_init; ++li)
    {
        const std::size_t gi = ib + li;
        mfem::Vector xi(sdim);
        for (int d = 0; d < sdim; ++d) { xi[d] = seeds.positions[gi][d]; }
        ClampAxisIfNeeded(xi, params.geom_tol, axisymmetric);
        for (int d = 0; d < sdim; ++d) { X((int)li, d) = xi[d]; }
        active.Tag(tag_seed)[(int)li] = (int)gi;
    }

    std::vector<TraceSummaryPOD> finished_local;
    finished_local.reserve(n_local_init);

    // Pathline logging
    std::vector<TrackPointPOD> track_log_local;
    if (save_path) track_log_local.reserve(N);

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

            ElectronExitCode code = ElectronExitCode::None;
            if (codes_opt) { code = (*codes_opt)[(std::size_t)j]; }

            AddFinished(finished_local, sid, code, xi, sdim);
        }

        active.RemoveParticles(rm);
    };
    
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

        // Collective point location + remote evaluation (GSLIB path)
        finder.FindPoints(active.Coords(), active.Coords().GetOrdering());
        const mfem::Array<int>       &codes = finder.GetCode();
        const mfem::Array<long long> &elem  = finder.GetElem();

        mfem::ParticleVector &Evals = active.Field(fld_E);
        finder.Interpolate(E_gf, Evals);

        mfem::ParticleVector &dir = active.Field(fld_dir);
        mfem::ParticleVector &Xc  = active.Coords();

        if (n_active > 0)
        {
            // Single-batch removal (preserves indexing correctness)
            mfem::Array<int> rm_all; rm_all.Reserve(n_active);
            std::vector<ElectronExitCode> rm_codes;
            rm_codes.reserve((std::size_t)n_active);

            // Single per-particle loop for the step
            for (int i = 0; i < n_active; ++i)
            {
                // Not found / invalid location
                if (elem[i] < 0 || codes[i] == 2)
                {
                    rm_all.Append(i);
                    rm_codes.push_back(ElectronExitCode::LeftVolume); 
                    continue;
                }

                // Compute direction from E
                double n2 = 0.0;
                bool ok = true;
                for (int d = 0; d < sdim; ++d)
                {
                    const double e = Evals(i, d);
                    ok = ok && std::isfinite(e);
                    n2 += e * e;
                }
                ok = ok && std::isfinite(n2) && (n2 > 1e-30);
                if (!ok)
                {
                    rm_all.Append(i);
                    rm_codes.push_back(ElectronExitCode::DegenerateTimeStep);
                    continue;
                }

                const double inv_nrm = 1.0 / std::sqrt(n2);
                double d2 = 0.0;
                for (int d = 0; d < sdim; ++d)
                {
                    const double v = -Evals(i, d) * inv_nrm;
                    ok = ok && std::isfinite(v);
                    dir(i, d) = v;
                    d2 += v * v;
                }
                if (!ok || !(d2 > 1e-30))
                {
                    rm_all.Append(i);
                    rm_codes.push_back(ElectronExitCode::DegenerateTimeStep);
                    continue;
                }

                // Euler update (only for valid particles)
                mfem::Vector xi(sdim);
                for (int d = 0; d < sdim; ++d) { xi[d] = Xc(i, d); }
                for (int d = 0; d < sdim; ++d) { xi[d] += ds * dir(i, d); }

                ClampAxisIfNeeded(xi, params.geom_tol, axisymmetric);

                for (int d = 0; d < sdim; ++d) { Xc(i, d) = xi[d]; }

                // Geometry exit check (uses updated position)
                const double r = Xc(i, 0);
                const double z = Xc(i, 1);
                if (!geom.Inside(r, z))
                {
                    rm_all.Append(i);
                    rm_codes.push_back(ClassifyExitOrLeftVolume(geom, xi, params.geom_tol));
                }
            }

            // Remove once (rm_all indices are ascending; rm_codes aligned 1:1 with rm_all)
            if (rm_all.Size() > 0)
            {
                MFEM_VERIFY((std::size_t)rm_all.Size() == rm_codes.size(),
                            "rm_all / rm_codes size mismatch");
                kill_indices(rm_all, &rm_codes);
            }
        }

        // Periodic collective redistribution
        if (redistribution_every > 0 && ((step % redistribution_every) == 0))
        {
            int n_local  = active.GetNParticles();
            int n_global = 0;
            MPI_Allreduce(&n_local, &n_global, 1, MPI_INT, MPI_SUM, comm);

            if (n_global > 0)
            {
                // 1) Remove not-found particles (uses current coords)
                finder.FindPoints(active.Coords(), active.Coords().GetOrdering());
                const auto &codesR = finder.GetCode();
                const auto &elemR  = finder.GetElem();

                n_local = active.GetNParticles();
                if (n_local > 0)
                {
                    mfem::Array<int> rm_nf;
                    rm_nf.Reserve(n_local);
                    for (int i = 0; i < n_local; ++i)
                    {
                        if (elemR[i] < 0 || codesR[i] == 2) { rm_nf.Append(i); }
                    }
                    if (rm_nf.Size() > 0) { kill_indices(rm_nf, nullptr); }
                }

                // 2) Compute owners and redistribute (collective)
                finder.FindPoints(active.Coords(), active.Coords().GetOrdering());
                const auto &procs = finder.GetProc();

                mfem::Array<unsigned int> procs2(active.GetNParticles());
                for (int i = 0; i < active.GetNParticles(); ++i) { procs2[i] = procs[i]; }

                active.Redistribute(procs2);
            }
            // Pathline logging
            if (save_path)
            {
                mfem::ParticleVector &X = active.Coords();
                const int np = active.GetNParticles();
                for (int i = 0; i < np; ++i)
                {
                    TrackPointPOD p;
                    p.seed_id = active.Tag(tag_seed)[i];
                    p.step    = (int32_t)step;
                    p.x[0] = X(i, 0);
                    p.x[1] = X(i, 1);
                    p.x[2] = (sdim == 3) ? X(i, 2) : 0.0;
                    track_log_local.push_back(p);
                }
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

    GatherFinishedToRoot(comm, rank, size,
                     finished_local,
                     save_path ? &track_log_local : nullptr,
                     n_seeds, sdim, out_results);
}

}

