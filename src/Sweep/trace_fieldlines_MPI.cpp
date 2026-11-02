#include "trace_fieldlines_MPI.h"

#include <chrono>
#include <climits>
#include <atomic>
#include <condition_variable>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <thread>

namespace
{
void BroadcastString(MPI_Comm comm, int root, std::string &value)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    unsigned long long size =
        rank == root ? static_cast<unsigned long long>(value.size()) : 0ULL;
    MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG_LONG, root, comm);

    if (rank != root)
    {
        value.resize(static_cast<std::size_t>(size));
    }

    unsigned long long offset = 0;
    while (offset < size)
    {
        const int chunk = static_cast<int>(
            std::min<unsigned long long>(size - offset, INT_MAX));
        MPI_Bcast(value.data() + static_cast<std::size_t>(offset),
                  chunk, MPI_CHAR, root, comm);
        offset += static_cast<unsigned long long>(chunk);
    }
}
}

namespace mpitracing
{
namespace
{
enum class TracePhase : int
{
    Idle = 0,
    FindPoints,
    DistributePointInfo,
    EvaluateGradients,
    DistributeValues
};

const char *TracePhaseName(TracePhase phase)
{
    switch (phase)
    {
        case TracePhase::FindPoints:          return "FindPoints";
        case TracePhase::DistributePointInfo: return "DistributePointInfo";
        case TracePhase::EvaluateGradients:   return "EvaluateGradients";
        case TracePhase::DistributeValues:    return "DistributeValues";
        case TracePhase::Idle:                return "between_field_calls";
    }
    return "unknown";
}

long long SteadyNowNanoseconds()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}
}

struct TraceProgressState
{
    std::atomic<long long> step{0};
    std::atomic<int> active{0};
    std::atomic<int> points{0};
    std::atomic<unsigned long long> field_calls_started{0};
    std::atomic<unsigned long long> field_calls_completed{0};
    std::atomic<int> phase{static_cast<int>(TracePhase::Idle)};
    std::atomic<long long> phase_started_ns{SteadyNowNanoseconds()};
    std::atomic<unsigned long long> find_points_ns{0};
    std::atomic<unsigned long long> distribute_info_ns{0};
    std::atomic<unsigned long long> gradients_ns{0};
    std::atomic<unsigned long long> distribute_values_ns{0};

    void SetPhase(TracePhase next, int current_points)
    {
        points.store(current_points, std::memory_order_relaxed);
        phase_started_ns.store(SteadyNowNanoseconds(), std::memory_order_relaxed);
        phase.store(static_cast<int>(next), std::memory_order_release);
    }
};

namespace
{
class TraceProgressHeartbeat
{
public:
    TraceProgressHeartbeat(TraceProgressState *state,
                           double interval_seconds,
                           long long max_steps,
                           bool variable_steps,
                           std::size_t initial_active)
        : state_(state), interval_seconds_(interval_seconds),
          max_steps_(max_steps), variable_steps_(variable_steps),
          initial_active_(initial_active),
          started_(std::chrono::steady_clock::now())
    {
        if (!state_ || interval_seconds_ <= 0.0) { return; }
        worker_ = std::thread([this]() { Run(); });
    }

    ~TraceProgressHeartbeat() { Stop(); }

    void Stop()
    {
        if (!worker_.joinable()) { return; }
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        wake_.notify_one();
        worker_.join();
    }

private:
    void Run()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        const auto interval = std::chrono::duration<double>(interval_seconds_);
        while (!wake_.wait_for(lock, interval, [this]() { return stop_; }))
        {
            lock.unlock();
            Print();
            lock.lock();
        }
    }

    void Print() const
    {
        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - started_).count();
        const long long step = state_->step.load(std::memory_order_relaxed);
        const int phase_value = state_->phase.load(std::memory_order_acquire);
        const auto phase = static_cast<TracePhase>(phase_value);
        const long long phase_started =
            state_->phase_started_ns.load(std::memory_order_relaxed);
        const double phase_elapsed = static_cast<double>(
            SteadyNowNanoseconds() - phase_started) * 1.0e-9;
        const auto completed =
            state_->field_calls_completed.load(std::memory_order_relaxed);

        std::cout << "[MPITracer heartbeat] local_step=" << step;
        if (max_steps_ > 0)
        {
            std::cout << "/" << max_steps_;
            if (step > 0)
            {
                const double remaining = elapsed *
                    static_cast<double>(max_steps_ - step) /
                    static_cast<double>(step);
                std::cout << " step_limit_eta_s=" << remaining;
            }
        }
        else
        {
            std::cout << (variable_steps_
                ? "/variable(local-element)"
                : "/unlimited");
        }
        std::cout << " active="
                  << state_->active.load(std::memory_order_relaxed)
                  << "/" << initial_active_
                  << " field_call="
                  << state_->field_calls_started.load(std::memory_order_relaxed)
                  << " completed_field_calls=" << completed
                  << " phase=" << TracePhaseName(phase)
                  << " phase_elapsed_s=" << phase_elapsed
                  << " points=" << state_->points.load(std::memory_order_relaxed)
                  << " elapsed_s=" << elapsed
                  << " cumulative_findpoints_s="
                  << state_->find_points_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << " cumulative_point_info_s="
                  << state_->distribute_info_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << " cumulative_gradients_s="
                  << state_->gradients_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << " cumulative_value_distribution_s="
                  << state_->distribute_values_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << std::endl;
    }

    TraceProgressState *state_ = nullptr;
    double interval_seconds_ = 0.0;
    long long max_steps_ = 0;
    bool variable_steps_ = false;
    std::size_t initial_active_ = 0;
    std::chrono::steady_clock::time_point started_;
    mutable std::mutex mutex_;
    std::condition_variable wake_;
    bool stop_ = false;
    std::thread worker_;
};
}

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
                               const mfem::ParGridFunction& field,
                               mfem::FindPointsGSLIB& finder,
                               int sdim,
                               TraceProgressState *progress)
    : mesh_(mesh), field_(field), finder_(finder), sdim_(sdim),
      progress_(progress)
{
    MFEM_VERIFY(field_.VectorDim() == sdim_,
                "MPI tracer electric field dimension mismatch.");
}

void FieldEvaluator::EvaluateE(const mfem::Vector& points,
                               int n_points,
                               mfem::Vector& E,
                               mfem::Vector *element_sizes)
{
    if (progress_)
    {
        progress_->field_calls_started.fetch_add(1, std::memory_order_relaxed);
        progress_->SetPhase(TracePhase::FindPoints, n_points);
    }
    auto phase_start = std::chrono::steady_clock::now();

    // 1. Locate all points on the distributed mesh.
    //    Assumes points are in byVDIM ordering.
    finder_.FindPoints(points, mfem::Ordering::byVDIM);
    if (progress_)
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - phase_start).count();
        progress_->find_points_ns.fetch_add(
            static_cast<unsigned long long>(elapsed), std::memory_order_relaxed);
        progress_->SetPhase(TracePhase::DistributePointInfo, n_points);
    }
    phase_start = std::chrono::steady_clock::now();

    // 2. Interpolate the projected electric field. This is deliberately the
    // same GridFunction consumed by the BOOST/VTK paths. In particular, L2
    // averaging at element interfaces is handled by FindPointsGSLIB rather
    // than selecting one raw H1 potential gradient arbitrarily.
    finder_.Interpolate(field_, E);
    const auto field_ordering = field_.FESpace()->GetOrdering();
    if (field_ordering != mfem::Ordering::byVDIM)
    {
        mfem::Ordering::Reorder(E, sdim_, field_ordering,
                                mfem::Ordering::byVDIM);
    }
    if (progress_)
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - phase_start).count();
        progress_->distribute_info_ns.fetch_add(
            static_cast<unsigned long long>(elapsed), std::memory_order_relaxed);
        progress_->SetPhase(TracePhase::DistributePointInfo, n_points);
    }
    phase_start = std::chrono::steady_clock::now();

    // 3. The local element identity is only needed for element-relative step
    // sizes. Keep this custom distribution separate from field interpolation.
    if (element_sizes)
    {
        finder_.DistributePointInfoToOwningMPIRanks(
            recv_elem_, recv_ref_, recv_code_);
        const int n_recv = recv_elem_.Size();
        recv_h_.SetSize(n_recv);
        for (int i = 0; i < n_recv; ++i)
        {
            if (recv_code_[i] == 2)
            {
                recv_h_[i] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            const int elem_id = static_cast<int>(recv_elem_[i]);
            // Use the minimum singular-value element dimension. The default
            // determinant-based characteristic size can be much too large in
            // anisotropic elements and step across a field-direction reversal.
            const double h = mesh_.GetElementSize(elem_id, 1);
            recv_h_[i] = (h > 0.0 && std::isfinite(h))
                ? h
                : std::numeric_limits<double>::quiet_NaN();
        }
        finder_.DistributeInterpolatedValues(recv_h_, 1,
                                             mfem::Ordering::byVDIM,
                                             *element_sizes);
    }
    if (progress_)
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - phase_start).count();
        progress_->gradients_ns.fetch_add(
            static_cast<unsigned long long>(elapsed), std::memory_order_relaxed);
        progress_->SetPhase(TracePhase::DistributeValues, n_points);
    }
    phase_start = std::chrono::steady_clock::now();

    // Field and optional element sizes are now in original point order.
    if (progress_)
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - phase_start).count();
        progress_->distribute_values_ns.fetch_add(
            static_cast<unsigned long long>(elapsed), std::memory_order_relaxed);
        progress_->field_calls_completed.fetch_add(1, std::memory_order_relaxed);
        progress_->SetPhase(TracePhase::Idle, n_points);
    }
}


// ---------------------------- Steppers/Integrators --------------------------
class EulerIntegrator : public Integrator {
public:
    /// No special initialisation required.
    void Initialize(mfem::ParticleSet&) override { }  

    /// Perform one Euler step for all active particles.
    void Step(mfem::ParticleSet& particles,
          double step_parameter,
          bool local_element_step,
          FieldEvaluator& field,
          mfem::Array<int>& to_remove,
          std::vector<ElectronExitCode>& exit_codes,
          mfem::Vector& step_lengths) override
    {
        const int n_active = particles.GetNParticles();
        const int sdim = particles.GetDim();

        // Always call EvaluateE (collective), even with n_active == 0.
        mfem::Vector& X = particles.Coords();  // may be empty
        mfem::Vector E_flat;
        mfem::Vector element_sizes;
        field.EvaluateE(X, n_active, E_flat,
                        local_element_step ? &element_sizes : nullptr);

        to_remove.SetSize(0);
        exit_codes.clear();
        step_lengths.SetSize(n_active);

        if (n_active == 0) return;

        MFEM_VERIFY(E_flat.Size() == n_active * sdim, "Field evaluator returned wrong size");

        to_remove.Reserve(n_active);
        exit_codes.reserve(n_active);

        double* X_data = X.GetData();           // layout: [x0,y0, x1,y1, ...]
        const double* E_data = E_flat.GetData(); // layout: [Ex0,Ey0, Ex1,Ey1, ...]

        for (int i = 0; i < n_active; ++i) {
            double* xi = X_data + i * sdim;
            const double* Ei = E_data + i * sdim;
            const double dt = local_element_step
                ? step_parameter * element_sizes[i]
                : step_parameter;
            step_lengths[i] =
                (dt > 0.0 && std::isfinite(dt)) ? dt : 0.0;

            // 3a. Check for degenerate field.
            bool ok = step_lengths[i] > 0.0;
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
                xi[d] += step_lengths[i] * v[d];
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
              double step_parameter,
              bool local_element_step,
              FieldEvaluator& field,
              mfem::Array<int>& to_remove,
              std::vector<ElectronExitCode>& exit_codes,
              mfem::Vector& step_lengths) override
    {
        const int n_active = particles.GetNParticles();
        const int sdim = particles.GetDim();
        mfem::ParticleVector& k1 = particles.Field(k1_idx);
        mfem::ParticleVector& k2 = particles.Field(k2_idx);
        mfem::ParticleVector& k3 = particles.Field(k3_idx);
        mfem::ParticleVector& k4 = particles.Field(k4_idx);
        mfem::Vector X_temp(n_active * sdim);
        mfem::Vector E_flat;
        mfem::Vector element_sizes;
        step_lengths.SetSize(n_active);
        std::vector<char> valid(n_active, 1);
        auto compute_k = [&](const mfem::Vector& values,
                             mfem::ParticleVector& k) {
            const double* E_data = values.GetData();
            for (int i = 0; i < n_active; ++i) {
                double norm2 = 0.0;
                bool ok = true;
                for (int d = 0; d < sdim; ++d) {
                    const double e = E_data[i * sdim + d];
                    ok = ok && std::isfinite(e);
                    norm2 += e * e;
                }
                ok = ok && std::isfinite(norm2) && norm2 > 1e-30;
                if (!ok) {
                    valid[i] = 0;
                    for (int d = 0; d < sdim; ++d)
                        k(i, d) = std::numeric_limits<double>::quiet_NaN();
                } else {
                    const double inv_norm = 1.0 / std::sqrt(norm2);
                    for (int d = 0; d < sdim; ++d)
                        k(i, d) = -E_data[i * sdim + d] * inv_norm;
                }
            }
        };

        field.EvaluateE(particles.Coords(), n_active, E_flat,
                        local_element_step ? &element_sizes : nullptr);
        compute_k(E_flat, k1);
        for (int i = 0; i < n_active; ++i) {
            const double dt = local_element_step
                ? step_parameter * element_sizes[i] : step_parameter;
            step_lengths[i] = (dt > 0.0 && std::isfinite(dt)) ? dt : 0.0;
            if (step_lengths[i] == 0.0) { valid[i] = 0; }
        }

        auto stage = [&](const mfem::ParticleVector& k_in,
                         double scale,
                         mfem::ParticleVector& k_out) {
            const double* x0 = particles.Coords().GetData();
            double* xt = X_temp.GetData();
            for (int i = 0; i < n_active; ++i)
                for (int d = 0; d < sdim; ++d) {
                    const int j = i * sdim + d;
                    xt[j] = x0[j] + scale * step_lengths[i] * k_in(i, d);
                }
            field.EvaluateE(X_temp, n_active, E_flat, nullptr);
            compute_k(E_flat, k_out);
        };

        stage(k1, 0.5, k2);
        stage(k2, 0.5, k3);
        stage(k3, 1.0, k4);

        const double* x0 = particles.Coords().GetData();
        double* xnew = X_temp.GetData();
        for (int i = 0; i < n_active; ++i)
            for (int d = 0; d < sdim; ++d) {
                const int j = i * sdim + d;
                xnew[j] = valid[i]
                    ? x0[j] + (step_lengths[i] / 6.0) *
                        (k1(i, d) + 2.0*k2(i, d) +
                         2.0*k3(i, d) + k4(i, d))
                    : x0[j];
            }
        particles.Coords() = X_temp;
        to_remove.SetSize(0);
        exit_codes.clear();
        for (int i = 0; i < n_active; ++i)
            if (!valid[i]) {
                to_remove.Append(i);
                exit_codes.push_back(ElectronExitCode::DegenerateTimeStep);
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
    bool save_path,
    double progress_interval_seconds)
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

    // ----- Step policy and traversal limit -----
    const bool local_element_step = params.local_element_step;
    const double legacy_ds = params.c_step * h_ref;
    MFEM_VERIFY(local_element_step ||
                (legacy_ds > 0.0 && std::isfinite(legacy_ds)),
                "TraceDistributed: global step size must be positive finite.");

    const long long max_steps = local_element_step
        ? 0
        : ComputeMaxStepsFromLimit(params, legacy_ds, params.max_traversals);
    const long long N = (max_steps > 0)
        ? max_steps
        : std::numeric_limits<long long>::max();
    const double max_path_length = params.max_traversals > 0.0
        ? params.max_traversals * std::fabs(params.z_max - params.z_min)
        : 0.0;
    const double step_parameter = local_element_step
        ? params.c_step
        : legacy_ds;

    // ----- Geometry -----
    const TpcGeometry geom(params);

    // ----- Particle set initialisation -----
    mfem::ParticleSet active(comm, (int)n_local_init, sdim, mfem::Ordering::byVDIM);
    const int tag_seed    = active.AddTag("seed_id");
    const int tag_owner   = active.AddTag("owner_rank");   // sticky owner
    const int tag_last_elem = active.AddTag("last_elem");  // (may be removed later)
    const int tag_has_last  = active.AddTag("has_last");
    int path_length_idx = -1;
    if (local_element_step)
    {
        path_length_idx = active.AddField(
            1, mfem::Ordering::byVDIM, "path_length");
        active.Field(path_length_idx) = 0.0;
    }

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

    // ----- Field evaluator and optional independent progress heartbeat -----
    TraceProgressState progress_state;
    TraceProgressState *progress =
        progress_interval_seconds > 0.0 ? &progress_state : nullptr;
    if (progress)
    {
        progress->active.store(static_cast<int>(n_local_init),
                               std::memory_order_relaxed);
    }
    TraceProgressHeartbeat heartbeat(progress,
                                     progress_interval_seconds,
                                     max_steps,
                                     local_element_step,
                                     n_local_init);
    FieldEvaluator evaluator(mesh, phi, finder, sdim, progress);

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
        if (max_steps > 0) { track_log_local.reserve(max_steps); }
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

        if (progress)
        {
            progress->step.store(step, std::memory_order_relaxed);
            progress->active.store(global_active, std::memory_order_relaxed);
        }

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
        mfem::Vector step_lengths;
        std::vector<double> coordinates_before_step;
        const mfem::ParticleVector &before = active.Coords();
        coordinates_before_step.assign(before.GetData(),
            before.GetData() + static_cast<std::size_t>(n_active * sdim));
        integrator->Step(active,
                         step_parameter,
                         local_element_step,
                         evaluator,
                         rm_fail,
                         fail_codes,
                         step_lengths);

        mfem::Vector actual_step_lengths(n_active);
        const mfem::ParticleVector &after_step = active.Coords();
        for (int i = 0; i < n_active; ++i)
        {
            double displacement2 = 0.0;
            bool finite = true;
            for (int d = 0; d < sdim; ++d)
            {
                const double delta = after_step(i, d) -
                    coordinates_before_step[static_cast<std::size_t>(i * sdim + d)];
                finite = finite && std::isfinite(delta);
                displacement2 += delta * delta;
            }
            actual_step_lengths[i] = finite ? std::sqrt(displacement2) : 0.0;
        }

        if (local_element_step)
        {
            MFEM_VERIFY(step_lengths.Size() == active.GetNParticles(),
                        "Integrator returned the wrong local-step count.");
            mfem::ParticleVector &path_length = active.Field(path_length_idx);
            for (int i = 0; i < active.GetNParticles(); ++i)
            {
                path_length(i, 0) += actual_step_lengths[i];
            }
        }

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

        // With element-local step sizes, max_traversals is a physical path
        // length limit rather than a global iteration count.
        if (local_element_step && max_path_length > 0.0 && n_active > 0)
        {
            mfem::Array<int> rm_limit;
            std::vector<ElectronExitCode> limit_codes;
            const mfem::ParticleVector &path_length =
                active.Field(path_length_idx);
            for (int i = 0; i < n_active; ++i)
            {
                if (path_length(i, 0) >= max_path_length)
                {
                    rm_limit.Append(i);
                    limit_codes.push_back(ElectronExitCode::MaxSteps);
                }
            }
            if (rm_limit.Size() > 0)
            {
                kill_indices(rm_limit, &limit_codes);
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
            AddFinished(finished_local, sid, ElectronExitCode::MaxSteps,
                        xi, sdim);
        }

        mfem::Array<int> all(active.GetNParticles());
        for (int i = 0; i < active.GetNParticles(); ++i) all[i] = i;
        active.RemoveParticles(all);
    }

    heartbeat.Stop();
    if (progress)
    {
        std::cout << "[MPITracer phase totals] field_calls="
                  << progress->field_calls_completed.load(std::memory_order_relaxed)
                  << " findpoints_s="
                  << progress->find_points_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << " point_info_s="
                  << progress->distribute_info_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << " gradients_s="
                  << progress->gradients_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << " value_distribution_s="
                  << progress->distribute_values_ns.load(std::memory_order_relaxed) * 1.0e-9
                  << std::endl;
    }

    // ----- Gather all results to root -----
    GatherFinishedToRoot(comm, rank, size,
                         finished_local,
                         save_path ? &track_log_local : nullptr,
                         n_seeds, sdim, out_results);
}

void TraceReplicated(
    MPI_Comm world_comm,
    mfem::ParMesh& replicated_mesh,
    mfem::FindPointsGSLIB& replicated_finder,
    const mfem::ParGridFunction& replicated_phi,
    const ElectronTraceParams& params,
    const Seeds& seeds,
    std::vector<ElectronTraceResult>& out_results,
    double h_ref,
    bool axisymmetric,
    int redistribution_every,
    bool print_progress,
    bool save_path)
{
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(world_comm, &world_rank);
    MPI_Comm_size(world_comm, &world_size);

    int replica_size = 0;
    MPI_Comm_size(replicated_mesh.GetComm(), &replica_size);
    MFEM_VERIFY(replica_size == 1,
                "TraceReplicated requires a mesh on MPI_COMM_SELF.");

    const std::size_t n_seeds = seeds.positions.size();
    MFEM_VERIFY(n_seeds <= static_cast<std::size_t>(INT32_MAX),
                "TraceReplicated seed count exceeds its 32-bit result id.");

    // Assign complete jobs in a strided pattern. Fixed-grid seeds are ordered,
    // so striding avoids systematically grouping the longest traces together.
    Seeds local_seeds;
    std::vector<std::size_t> global_ids;
    const std::size_t local_capacity =
        (n_seeds + static_cast<std::size_t>(world_size) - 1U) /
        static_cast<std::size_t>(world_size);
    local_seeds.positions.reserve(local_capacity);
    global_ids.reserve(local_capacity);
    const bool have_volumes = seeds.volumes.size() == n_seeds;
    if (have_volumes) { local_seeds.volumes.reserve(local_capacity); }

    for (std::size_t global_id = static_cast<std::size_t>(world_rank);
         global_id < n_seeds;
         global_id += static_cast<std::size_t>(world_size))
    {
        local_seeds.positions.push_back(seeds.positions[global_id]);
        if (have_volumes) { local_seeds.volumes.push_back(seeds.volumes[global_id]); }
        global_ids.push_back(global_id);
    }

    constexpr double progress_interval_seconds = 30.0;
    if (print_progress && world_rank == 0)
    {
        std::cout << "[MPITracer] mode=replicated-per-rank"
                  << " global_electrons=" << n_seeds
                  << " world_ranks=" << world_size
                  << " rank0_electrons=" << local_seeds.positions.size()
                  << " field_source="
                  << (params.use_l2_electric_field
                      ? "canonical_L2_discrete_gradient"
                      : "projected_H1_electric_field")
                  << " step_mode="
                  << (params.local_element_step
                      ? "current_element_minimum_dimension"
                      : "global_minimum_edge")
                  << " c_step=" << params.c_step
                  << " global_min_edge_reference=" << h_ref;
        if (params.local_element_step)
        {
            const double max_path_length = params.max_traversals > 0.0
                ? params.max_traversals *
                    std::fabs(params.z_max - params.z_min)
                : 0.0;
            std::cout << " step_size=variable"
                      << " max_path_length=";
            if (max_path_length > 0.0) { std::cout << max_path_length; }
            else { std::cout << "unlimited"; }
        }
        else
        {
            const double ds = params.c_step * h_ref;
            const long long max_steps =
                ComputeMaxStepsFromLimit(params, ds, params.max_traversals);
            std::cout << " step_size=" << ds
                      << " max_steps=";
            if (max_steps > 0) { std::cout << max_steps; }
            else { std::cout << "unlimited"; }
        }
        std::cout << " progress_interval_s=" << progress_interval_seconds
                  << " configured_redistribution_every=" << redistribution_every
                  << " effective_redistribution_every=0"
                  << " eta_basis="
                  << (params.local_element_step
                      ? "physical_path_limit_variable_steps"
                      : "rank0_step_rate_to_step_limit")
                  << std::endl;
    }

    const auto trace_start = std::chrono::steady_clock::now();
    std::vector<ElectronTraceResult> local_results;

    // This invokes the existing, known-correct evaluator/integrator unchanged.
    // Every collective inside it is confined to the rank's MPI_COMM_SELF.
    TraceDistributed(replicated_mesh,
                     replicated_finder,
                     replicated_phi,
                     params,
                     local_seeds,
                     local_results,
                     h_ref,
                     axisymmetric,
                     /*redistribution_every=*/0,
                     /*debug=*/false,
                     /*debug_every=*/0,
                     save_path,
                     (print_progress && world_rank == 0)
                         ? progress_interval_seconds
                         : 0.0);

    const auto trace_stop = std::chrono::steady_clock::now();
    if (print_progress && world_rank == 0)
    {
        const double seconds =
            std::chrono::duration<double>(trace_stop - trace_start).count();
        std::cout << "[MPITracer replicated rank " << world_rank << "] seeds="
                  << local_seeds.positions.size()
                  << " local_trace_s=" << seconds << std::endl;
    }

    MFEM_VERIFY(local_results.size() == local_seeds.positions.size(),
                "Replicated local trace returned the wrong result count.");

    std::vector<TraceSummaryPOD> finished_local;
    std::vector<TrackPointPOD> track_log_local;
    finished_local.reserve(local_results.size());

    for (std::size_t local_id = 0; local_id < local_results.size(); ++local_id)
    {
        const auto &result = local_results[local_id];
        if (result.points.empty()) { continue; }

        const int global_id = static_cast<int>(global_ids[local_id]);
        AddFinished(finished_local, global_id, result.exit_code,
                    result.points.back(), replicated_mesh.SpaceDimension());

        if (save_path)
        {
            for (std::size_t step = 0; step < result.points.size(); ++step)
            {
                MFEM_VERIFY(step <= static_cast<std::size_t>(INT32_MAX),
                            "Replicated trace path exceeds its 32-bit step id.");
                TrackPointPOD point{};
                point.seed_id = static_cast<int32_t>(global_id);
                point.step = static_cast<int32_t>(step);
                for (int d = 0; d < replicated_mesh.SpaceDimension(); ++d)
                {
                    point.x[d] = result.points[step][d];
                }
                track_log_local.push_back(point);
            }
        }
    }

    // This is the only MPI_COMM_WORLD synchronization in the tracing call.
    GatherFinishedToRoot(world_comm, world_rank, world_size,
                         finished_local,
                         save_path ? &track_log_local : nullptr,
                         n_seeds,
                         replicated_mesh.SpaceDimension(),
                         out_results);
}
}

void MPITraceContext::Build(mfem::ParMesh &mesh, bool debug)
{
    h_ref = ComputeGlobalMinEdgeLength(mesh, debug);

    if (mesh.GetNodes() == nullptr)
    {
        const int order = 1;
        const bool discont = false;
        const int space_dim = mesh.SpaceDimension();
        const int ordering = mfem::Ordering::byVDIM;
        mesh.SetCurvature(order, discont, space_dim, ordering);
    }

    finder = std::make_unique<mfem::FindPointsGSLIB>(mesh.GetComm());
    finder->Setup(mesh);
    finder->SetDefaultInterpolationValue(
        std::numeric_limits<double>::quiet_NaN());
    // Average the two polynomial representations on partition interfaces.
    finder->SetL2AvgType(
        mfem::FindPointsGSLIB::AvgType::ARITHMETIC);
}

void ReplicatedMPITraceContext::Build(
    mfem::ParMesh &source_mesh,
    const mfem::ParGridFunction &source_field,
    bool debug)
{
    MPI_Comm world_comm = source_mesh.GetComm();
    int world_rank = 0;
    MPI_Comm_rank(world_comm, &world_rank);

    // Match the existing MPI context preparation before serializing. GSLIB's
    // point locator must see the same nodal representation in both modes.
    if (source_mesh.GetNodes() == nullptr)
    {
        const int order = 1;
        const bool discont = false;
        const int space_dim = source_mesh.SpaceDimension();
        const int ordering = mfem::Ordering::byVDIM;
        source_mesh.SetCurvature(order, discont, space_dim, ordering);
    }

    // Preserve the existing step-size policy exactly. This collective occurs
    // once during setup, before any independent tracing begins.
    h_ref = ComputeGlobalMinEdgeLength(source_mesh, debug);

    constexpr int serialization_root = 0;
    std::string mesh_data;
    std::string field_data;
    std::string serialization_error;
    int serialization_ok = 1;
    {
        // The canonical electric field has already been computed in parallel
        // on source_mesh. Gather it once, then distribute the completed field
        // instead of assembling a full-mesh gradient independently per rank.
        mfem::Mesh serial_mesh =
            source_mesh.GetSerialMesh(serialization_root);
        mfem::GridFunction serial_field =
            source_field.GetSerialGridFunction(
                serialization_root, serial_mesh);

        if (world_rank == serialization_root)
        {
            try
            {
                std::ostringstream mesh_stream;
                mesh_stream << std::setprecision(17);
                serial_mesh.Print(mesh_stream);
                mesh_data = mesh_stream.str();

                std::ostringstream field_stream;
                field_stream << std::setprecision(17);
                serial_field.Save(field_stream);
                field_data = field_stream.str();
            }
            catch (const std::exception &e)
            {
                serialization_ok = 0;
                serialization_error = e.what();
            }
        }
    }

    MPI_Bcast(&serialization_ok, 1, MPI_INT,
              serialization_root, world_comm);
    BroadcastString(world_comm, serialization_root, serialization_error);
    if (!serialization_ok)
    {
        throw std::runtime_error(
            "Unable to serialize replicated tracing state: " +
            serialization_error);
    }

    BroadcastString(world_comm, serialization_root, mesh_data);
    BroadcastString(world_comm, serialization_root, field_data);
    serialized_mesh_bytes = mesh_data.size();
    serialized_field_bytes = field_data.size();

    std::string local_error;
    int local_ok = 1;
    try
    {
        {
            std::istringstream mesh_stream(mesh_data);
            mfem::Mesh local_serial_mesh(mesh_stream,
                                         /*generate_edges=*/1,
                                         /*refine=*/1,
                                         /*fix_orientation=*/true);
            replicated_mesh = std::make_unique<mfem::ParMesh>(
                MPI_COMM_SELF, local_serial_mesh);
        }
        mesh_data.clear();
        mesh_data.shrink_to_fit();

        {
            std::istringstream field_stream(field_data);
            replicated_field = std::make_unique<mfem::ParGridFunction>(
                replicated_mesh.get(), field_stream);
        }
        field_data.clear();
        field_data.shrink_to_fit();
        replicated_field->ExchangeFaceNbrData();

        finder = std::make_unique<mfem::FindPointsGSLIB>(MPI_COMM_SELF);
        finder->Setup(*replicated_mesh);
        finder->SetDefaultInterpolationValue(
            std::numeric_limits<double>::quiet_NaN());
        finder->SetL2AvgType(
            mfem::FindPointsGSLIB::AvgType::ARITHMETIC);
    }
    catch (const std::exception &e)
    {
        local_ok = 0;
        local_error = e.what();
        finder.reset();
        replicated_field.reset();
        replicated_mesh.reset();
    }

    int all_ok = 0;
    MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, world_comm);
    if (!all_ok)
    {
        if (!local_ok)
        {
            std::cerr << "[MPITracer replicated rank " << world_rank
                      << "] reconstruction failed: " << local_error
                      << std::endl;
        }
        throw std::runtime_error(
            "Replicated tracing state reconstruction failed on at least one MPI rank.");
    }

    if (world_rank == 0)
    {
        const auto payload_mib =
            static_cast<double>(serialized_mesh_bytes + serialized_field_bytes) /
            (1024.0 * 1024.0);
        std::cout << "[MPITracer] mode=replicated-per-rank"
                  << " elements_per_rank=" << replicated_mesh->GetNE()
                  << " field_dofs_per_rank=" << replicated_field->Size()
                  << " serialized_payload_MiB=" << payload_mib
                  << std::endl;
    }
}

