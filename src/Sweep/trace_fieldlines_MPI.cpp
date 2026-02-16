#include "trace_fieldlines_MPI.h"

#include <chrono>

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

// -------------------------------- Field Eval -----------------------------
FieldEvaluator::FieldEvaluator(mfem::ParMesh& mesh,
                               const mfem::ParGridFunction& phi,
                               mfem::FindPointsGSLIB& finder,
                               int sdim)
    : mesh_(mesh), phi_(phi), finder_(finder), sdim_(sdim)
{
    // The finder must already be attached to the mesh.
    // No further initialisation required.
}

void FieldEvaluator::EvaluateE(const mfem::Vector& points,
                               int n_points,
                               mfem::Vector& E)
{
    // 1. Locate all points on the distributed mesh.
    //    Assumes points are in byVDIM ordering.
    finder_.FindPoints(points, mfem::Ordering::byVDIM);

    // 2. Distribute point information to the ranks that own the elements.
    finder_.DistributePointInfoToOwningMPIRanks(recv_elem_, recv_ref_, recv_code_);
    const int n_recv = recv_elem_.Size();
    MFEM_VERIFY(recv_ref_.Size() == n_recv * sdim_, "Invalid reference coordinates size");

    // 3. Compute electric field on the owning ranks.
    recv_E_.SetSize(n_recv * sdim_);
    mfem::Vector grad(sdim_);

    for (int i = 0; i < n_recv; ++i) {
        if (recv_code_[i] == 2) { // point not found
            for (int d = 0; d < sdim_; ++d)
                recv_E_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // Element is local to this rank (guaranteed by the distribute call)
        int elem_id = (int)recv_elem_[i];
        mfem::ElementTransformation* T = mesh_.GetElementTransformation(elem_id);
        MFEM_VERIFY(T != nullptr, "Null element transformation");

        // Reference coordinates (GSLIB returns in [0,1] for MFEM mesh)
        const double* ref = recv_ref_.GetData() + i * sdim_;
        mfem::IntegrationPoint ip;
        if (sdim_ == 2) ip.Set2(ref[0], ref[1]);
        else            ip.Set3(ref[0], ref[1], ref[2]);

        T->SetIntPoint(&ip);
        phi_.GetGradient(*T, grad);
        grad *= -1.0;   // E = -∇φ

        for (int d = 0; d < sdim_; ++d)
            recv_E_[i * sdim_ + d] = grad[d];
    }

    // 4. Distribute computed E values back to the original ranks.
    //    The output 'E' will be ordered consistently with the input 'points'.
    finder_.DistributeInterpolatedValues(recv_E_, sdim_,
                                         mfem::Ordering::byVDIM, E);
}


// ---------------------------- Steppers/Integrators --------------------------
class EulerIntegrator : public Integrator {
public:
    /// No special initialisation required.
    void Initialize(mfem::ParticleSet&) override { }  

    /// Perform one Euler step for all active particles.
    void Step(mfem::ParticleSet& particles,
          double dt,
          FieldEvaluator& field,
          mfem::Array<int>& to_remove,
          std::vector<ElectronExitCode>& exit_codes) override
    {
        const int n_active = particles.GetNParticles();
        const int sdim = particles.GetDim();

        // Always call EvaluateE (collective), even with n_active == 0.
        mfem::Vector& X = particles.Coords();  // may be empty
        mfem::Vector E_flat;
        field.EvaluateE(X, n_active, E_flat);

        to_remove.SetSize(0);
        exit_codes.clear();

        if (n_active == 0) return;

        MFEM_VERIFY(E_flat.Size() == n_active * sdim, "Field evaluator returned wrong size");

        to_remove.Reserve(n_active);
        exit_codes.reserve(n_active);

        double* X_data = X.GetData();           // layout: [x0,y0, x1,y1, ...]
        const double* E_data = E_flat.GetData(); // layout: [Ex0,Ey0, Ex1,Ey1, ...]

        for (int i = 0; i < n_active; ++i) {
            double* xi = X_data + i * sdim;
            const double* Ei = E_data + i * sdim;

            // 3a. Check for degenerate field.
            bool ok = true;
            double n2 = 0.0;
            for (int d = 0; d < sdim; ++d) {
                const double e = Ei[d];
                ok = ok && std::isfinite(e);
                n2 += e * e;
            }
            ok = ok && std::isfinite(n2) && (n2 > 1e-30);
            if (!ok) {
                to_remove.Append(i);
                exit_codes.push_back(ElectronExitCode::DegenerateTimeStep);
                continue;
            }

            const double inv_nrm = 1.0 / std::sqrt(n2);
            double v[3] = {0.0, 0.0, 0.0}; 
            double d2 = 0.0;
            for (int d = 0; d < sdim; ++d) {
                const double vd = -Ei[d] * inv_nrm;
                v[d] = vd;
                d2 += vd * vd;
            }
            if (!(d2 > 1e-30)) {
                to_remove.Append(i);
                exit_codes.push_back(ElectronExitCode::DegenerateTimeStep);
                continue;
            }

            // 3c. Euler update: x = x + dt * v.
            for (int d = 0; d < sdim; ++d) {
                xi[d] += dt * v[d];
            }

        }
    }
};

/// 4th‑order Runge–Kutta integrator for drift velocity v = -E/|E|.
class RK4Integrator : public Integrator {
private:
    int k1_idx, k2_idx, k3_idx, k4_idx;

public:
    void Initialize(mfem::ParticleSet& particles) override {
        const int sdim = particles.GetDim();
        const auto ord = particles.Coords().GetOrdering();
        k1_idx = particles.AddField(sdim, ord, "k1");
        k2_idx = particles.AddField(sdim, ord, "k2");
        k3_idx = particles.AddField(sdim, ord, "k3");
        k4_idx = particles.AddField(sdim, ord, "k4");
    }

    void Step(mfem::ParticleSet& particles,
              double dt,
              FieldEvaluator& field,
              mfem::Array<int>& to_remove,
              std::vector<ElectronExitCode>& exit_codes) override
    {
        const int n_active = particles.GetNParticles();
        const int sdim = particles.GetDim();

        if (n_active == 0) {
            // Collective empty calls
            mfem::Vector empty;
            mfem::Vector dummy;
            for (int i = 0; i < 4; ++i) field.EvaluateE(empty, 0, dummy);
            to_remove.SetSize(0);
            exit_codes.clear();
            return;
        }

        // Access fields using stored indices
        mfem::ParticleVector& k1 = particles.Field(k1_idx);
        mfem::ParticleVector& k2 = particles.Field(k2_idx);
        mfem::ParticleVector& k3 = particles.Field(k3_idx);
        mfem::ParticleVector& k4 = particles.Field(k4_idx);

        mfem::Vector X_temp(n_active * sdim);
        mfem::Vector E_flat;

        std::vector<char> valid(n_active, 1);

        auto compute_k = [&](const mfem::Vector& E_flat, mfem::ParticleVector& k) {
            const double* E_data = E_flat.GetData();
            for (int i = 0; i < n_active; ++i) {
                double* ki = &k(i, 0);
                double norm2 = 0.0;
                bool ok = true;
                for (int d = 0; d < sdim; ++d) {
                    double e = E_data[i * sdim + d];
                    ok = ok && std::isfinite(e);
                    norm2 += e * e;
                }
                ok = ok && (norm2 > 1e-30);
                if (!ok) {
                    valid[i] = 0;
                    for (int d = 0; d < sdim; ++d) ki[d] = std::numeric_limits<double>::quiet_NaN();
                } else {
                    double inv_norm = 1.0 / std::sqrt(norm2);
                    for (int d = 0; d < sdim; ++d) ki[d] = -E_data[i * sdim + d] * inv_norm;
                }
            }
        };

        // Stage 1
        field.EvaluateE(particles.Coords(), n_active, E_flat);
        compute_k(E_flat, k1);

        // Stage 2
        {
            double* Xt = X_temp.GetData();
            const double* X0 = particles.Coords().GetData();
            for (int i = 0; i < n_active; ++i) {
                const double* k1i = &k1(i, 0);
                double* xt = Xt + i * sdim;
                const double* x0 = X0 + i * sdim;
                for (int d = 0; d < sdim; ++d)
                    xt[d] = x0[d] + 0.5 * dt * k1i[d];
            }
            field.EvaluateE(X_temp, n_active, E_flat);
            compute_k(E_flat, k2);
        }

        // Stage 3
        {
            double* Xt = X_temp.GetData();
            const double* X0 = particles.Coords().GetData();
            for (int i = 0; i < n_active; ++i) {
                const double* k2i = &k2(i, 0);
                double* xt = Xt + i * sdim;
                const double* x0 = X0 + i * sdim;
                for (int d = 0; d < sdim; ++d)
                    xt[d] = x0[d] + 0.5 * dt * k2i[d];
            }
            field.EvaluateE(X_temp, n_active, E_flat);
            compute_k(E_flat, k3);
        }

        // Stage 4
        {
            double* Xt = X_temp.GetData();
            const double* X0 = particles.Coords().GetData();
            for (int i = 0; i < n_active; ++i) {
                const double* k3i = &k3(i, 0);
                double* xt = Xt + i * sdim;
                const double* x0 = X0 + i * sdim;
                for (int d = 0; d < sdim; ++d)
                    xt[d] = x0[d] + dt * k3i[d];
            }
            field.EvaluateE(X_temp, n_active, E_flat);
            compute_k(E_flat, k4);
        }

        // Final combination into X_temp
        double* Xnew = X_temp.GetData();
        const double* X0 = particles.Coords().GetData();
        for (int i = 0; i < n_active; ++i) {
            if (valid[i]) {
                const double* k1i = &k1(i, 0);
                const double* k2i = &k2(i, 0);
                const double* k3i = &k3(i, 0);
                const double* k4i = &k4(i, 0);
                double* xnew = Xnew + i * sdim;
                const double* x0 = X0 + i * sdim;
                for (int d = 0; d < sdim; ++d)
                    xnew[d] = x0[d] + (dt / 6.0) *
                            (k1i[d] + 2.0 * k2i[d] + 2.0 * k3i[d] + k4i[d]);
            } else {
                double* xnew = Xnew + i * sdim;
                for (int d = 0; d < sdim; ++d)
                    xnew[d] = std::numeric_limits<double>::quiet_NaN();
            }
        }

        // Copy back
        double* X_dst = particles.Coords().GetData();
        for (int i = 0; i < n_active; ++i) {
            const double* src = Xnew + i * sdim;
            double* dst = X_dst + i * sdim;
            for (int d = 0; d < sdim; ++d)
                dst[d] = src[d];
        }

        // Post-update NaN check
        for (int i = 0; i < n_active; ++i) {
            if (valid[i]) {
                double* xi = X_dst + i * sdim;
                bool ok = true;
                for (int d = 0; d < sdim; ++d)
                    ok = ok && std::isfinite(xi[d]);
                if (!ok) valid[i] = 0;
            }
        }

        // Build removal list
        to_remove.SetSize(0);
        exit_codes.clear();
        for (int i = 0; i < n_active; ++i) {
            if (!valid[i]) {
                to_remove.Append(i);
                exit_codes.push_back(ElectronExitCode::DegenerateTimeStep);
            }
        }
    }
};

// ---------------------------- Particle Handling ----------------------------
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


// ------------------------------- Tracing ----------------------------------
void TraceDistributed(
    mfem::ParMesh& mesh,
    mfem::FindPointsGSLIB& finder,
    const mfem::ParGridFunction& phi,
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
    const std::size_t n_seeds = seeds.positions.size();

    // ----- Distribute seeds among ranks -----
    std::size_t ib = 0, ie = 0;
    SeedRangeForRank(n_seeds, rank, size, ib, ie);
    const std::size_t n_local_init = ie - ib;

    // ----- Step size and max steps -----
    const double ds = params.c_step * h_ref;
    MFEM_VERIFY(ds > 0.0 && std::isfinite(ds),
                "TraceDistributed: ds must be positive finite.");

    const long long max_steps = ComputeMaxStepsFromLimit(params, ds, params.max_traversals);
    const long long N = (max_steps > 0) ? max_steps : std::numeric_limits<long long>::max();

    // ----- Geometry -----
    const TpcGeometry geom(params);

    // ----- Particle set initialisation -----
    mfem::ParticleSet active(comm, (int)n_local_init, sdim, mfem::Ordering::byVDIM);
    const int tag_seed    = active.AddTag("seed_id");
    const int tag_owner   = active.AddTag("owner_rank");   // sticky owner
    const int tag_last_elem = active.AddTag("last_elem");  // (may be removed later)
    const int tag_has_last  = active.AddTag("has_last");

    // Initialise particles from local seeds
    mfem::ParticleVector& X = active.Coords();
    for (std::size_t li = 0; li < n_local_init; ++li)
    {
        const std::size_t gi = ib + li;
        mfem::Vector xi(sdim);
        for (int d = 0; d < sdim; ++d) { xi[d] = seeds.positions[gi][d]; }
        ClampAxisIfNeeded(xi, params.geom_tol, axisymmetric);
        for (int d = 0; d < sdim; ++d) { X((int)li, d) = xi[d]; }

        active.Tag(tag_seed)[(int)li] = (int)gi;
        active.Tag(tag_owner)[(int)li] = rank;
        active.Tag(tag_has_last)[(int)li]  = 0;
        active.Tag(tag_last_elem)[(int)li] = -1;
    }

    // ----- Field evaluator -----
    FieldEvaluator evaluator(mesh, phi, finder, sdim);

    // ----- Integrator selection -----
    std::unique_ptr<Integrator> integrator;
    // Assume params contains an enum IntegratorType; default to Euler.
    if (params.method == "Euler-Cauchy") {
      integrator = std::make_unique<EulerIntegrator>();
    }
    else if (params.method == "RK4") {
      integrator = std::make_unique<RK4Integrator>();
    }
    // else if (...) { integrator = std::make_unique<RK4Integrator>(); }

    integrator->Initialize(active);   // let integrator add any required fields

    // ----- Termination storage -----
    std::vector<TraceSummaryPOD> finished_local;
    finished_local.reserve(n_local_init);

    // ----- Path logging (optional) -----
    std::vector<TrackPointPOD> track_log_local;
    if (save_path) {
        track_log_local.reserve(N);   // may be large; consider not reserving
    }

    // ----- Helper lambda to remove particles and record summaries -----
    auto kill_indices =
        [&](const mfem::Array<int> &rm, const std::vector<ElectronExitCode> *codes_opt)
    {
        if (rm.Size() == 0) return;

        mfem::ParticleVector& Xk = active.Coords();
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

    // ----- Main time stepping loop -----
    for (long long step = 0; step < N; ++step)
    {
        int n_active = active.GetNParticles();
        int global_active = 0;
        MPI_Allreduce(&n_active, &global_active, 1, MPI_INT, MPI_SUM, comm);
        if (global_active == 0) break;

        if (debug && debug_every > 0 && (step % debug_every == 0))
        {
            std::cout << "[rank " << rank << "] step=" << step
                    << " local_active=" << n_active
                    << " global_active=" << global_active
                    << std::endl;
        }

        // ----- 1. Integrator step -----
        mfem::Array<int> rm_fail;
        std::vector<ElectronExitCode> fail_codes;
        integrator->Step(active, ds, evaluator, rm_fail, fail_codes);

        if (rm_fail.Size() > 0)
        {
            kill_indices(rm_fail, &fail_codes);
            n_active = active.GetNParticles();
        }

        // ----- 2. Axis clamping -----
        if (axisymmetric && n_active > 0)
        {
            mfem::ParticleVector& Xc = active.Coords();
            for (int i = 0; i < n_active; ++i)
            {
                double& r = Xc(i, 0);
                if (r <= params.geom_tol) { r = params.geom_tol; }
            }
        }
        // ----- 3. Geometry exit check -----
        if (n_active > 0)
        {
            mfem::Array<int> rm_geom;
            std::vector<ElectronExitCode> geom_codes;
            rm_geom.Reserve(n_active);
            geom_codes.reserve(n_active);

            mfem::ParticleVector& Xc = active.Coords();
            for (int i = 0; i < n_active; ++i)
            {
                const double r = Xc(i, 0);
                const double z = Xc(i, 1);
                if (!geom.Inside(r, z))
                {
                    rm_geom.Append(i);
                    mfem::Vector xi(&Xc(i,0), sdim);
                    geom_codes.push_back(ClassifyExitOrLeftVolume(geom, xi, params.geom_tol));
                }
            }

            if (rm_geom.Size() > 0)
            {
                kill_indices(rm_geom, &geom_codes);
                n_active = active.GetNParticles();
            }
        }

        // ----- 4. Path logging -----
        if (save_path && n_active > 0)
        {
            mfem::ParticleVector& Xlog = active.Coords();
            for (int i = 0; i < n_active; ++i)
            {
                TrackPointPOD p;
                p.seed_id = active.Tag(tag_seed)[i];
                p.step    = (int32_t)step;
                p.x[0]    = Xlog(i, 0);
                p.x[1]    = Xlog(i, 1);
                p.x[2]    = (sdim == 3) ? Xlog(i, 2) : 0.0;
                track_log_local.push_back(p);
            }
        }

        // ----- 5. Redistribution -----
        if (redistribution_every > 0 && (step % redistribution_every) == 0)
        {
            // Locate all remaining particles (collective)
            finder.FindPoints(active.Coords(), mfem::Ordering::byVDIM);
            const auto& codes = finder.GetCode();
            const auto& procs = finder.GetProc();

            // Remove particles that are now outside the mesh (local operation)
            mfem::Array<int> rm_outside;
            const int np = active.GetNParticles();
            for (int i = 0; i < np; ++i)
            {
                if (codes[i] == 2) { rm_outside.Append(i); }
            }
            if (rm_outside.Size() > 0)
            {
                kill_indices(rm_outside, nullptr);
            }

            // After removal, get fresh location data (collective – must be called by all ranks)
            finder.FindPoints(active.Coords(), mfem::Ordering::byVDIM);

            // Redistribute particles to their element owners (collective – all ranks call it)
            active.Redistribute(finder.GetProc());
        }
        n_active = active.GetNParticles();

        if (debug && debug_every > 0 && (step % debug_every == 0))
        {
            std::cout << "[rank " << rank << "] step " << step
                    << std::endl;
        }
    } // end of time loop

    // ----- After loop: remove any particles still active (MaxSteps) -----
    if (active.GetNParticles() > 0)
    {
        mfem::ParticleVector& Xf = active.Coords();
        for (int i = 0; i < active.GetNParticles(); ++i)
        {
            const int sid = active.Tag(tag_seed)[i];
            mfem::Vector xi(sdim);
            for (int d = 0; d < sdim; ++d) { xi[d] = Xf(i, d); }
            AddFinished(finished_local, sid, ElectronExitCode::MaxSteps, xi, sdim);
        }

        mfem::Array<int> all(active.GetNParticles());
        for (int i = 0; i < active.GetNParticles(); ++i) all[i] = i;
        active.RemoveParticles(all);
    }

    // ----- Gather all results to root -----
    GatherFinishedToRoot(comm, rank, size,
                         finished_local,
                         save_path ? &track_log_local : nullptr,
                         n_seeds, sdim, out_results);
}
}

