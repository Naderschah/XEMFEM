#include <mfem.hpp>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <indicators/progress_bar.hpp>
#include <atomic>

#include <stdexcept>

#include "trace_fieldlines.h"

using namespace mfem;


// ------------------------------------- Commmon ----------------
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


static void ValidateElectronExitCodesAllSeeds(const std::vector<ElectronTraceResult> &out_results)
{
    bool had_max_steps = false;
    bool had_none_exit = false;

    // If you truly guarantee "all particles are traced", then any missing/invalid
    // trace should be represented by an exit_code you can detect (often None).
    for (const auto &res : out_results)
    {
        if (res.exit_code == ElectronExitCode::MaxSteps) { had_max_steps = true; }
        else if (res.exit_code == ElectronExitCode::None) { had_none_exit = true; }
    }

    if (had_max_steps)
    {
        std::cout << "[WARNING] Max Steps Exit Condition for electron tracing" << std::endl;
    }
    if (had_none_exit)
    {
        std::cout << "[WARNING] None Exit Condition for electron tracing" << std::endl;
    }
}


// ----------------------------------------- Own Tracer ----------

using state_t = std::vector<double>;

static inline void ClampAxis(state_t &x, double geom_tol) noexcept
{
    if (x[0] <= geom_tol) x[0] = geom_tol;
}

static inline void MaybePush(ElectronTraceResult &out, const state_t &x, bool save_pathlines)
{
    const int sdim = x.size();
    if (!save_pathlines) return;
    mfem::Vector p(sdim);
    for (int d = 0; d < sdim; ++d) { p[d] = x[d]; }
    out.points.push_back(std::move(p));
}

static inline ElectronExitCode ClassifyExit(const TpcGeometry &geom,
                                     const state_t &x,
                                     double geom_tol) noexcept
{
    ElectronExitCode c = geom.ClassifyBoundary(x[0], x[1], geom_tol);
    if (c == ElectronExitCode::None) c = ElectronExitCode::LeftVolume;
    return c;
}

// Evaluate -E(x) using FindPointsGSLIB in "batched" form (npts=1).
struct FieldRHS
{
    mfem::FindPointsGSLIB       &finder;
    const mfem::ParGridFunction &E_gf;

    mfem::Vector points;   // size = sdim
    mfem::Vector E;        // size = vdim

    const int sdim;
    const int vdim;
    const double geom_tol;

    FieldRHS(mfem::FindPointsGSLIB &f, const mfem::ParGridFunction &gf,
             int sdim_, double tol)
      : finder(f), E_gf(gf), points(sdim_), E(gf.VectorDim()),
        sdim(sdim_), vdim(gf.VectorDim()), geom_tol(tol) {}

    inline void operator()(const state_t &xin, state_t &dxds, double)
    {
        dxds.resize(sdim);

        for (int d = 0; d < sdim; ++d) points[d] = xin[d];

        finder.FindPoints(points, mfem::Ordering::byNODES);

        const auto &elem_ids = finder.GetElem();
        const auto &codes    = finder.GetCode();
        
        // FIXME Handle boundary (exit code 1)
        if (elem_ids.Size() != 1 || codes.Size() != 1 || elem_ids[0] < 0 || codes[0] == 2)
        {
            std::fill(dxds.begin(), dxds.end(), 0.0);
            // Get Debug information
            // Reference coordinates (MFEM + GSLIB)
            const mfem::Vector &ref  = finder.GetReferencePosition();
            const mfem::Vector &gref = finder.GetGSLIBReferencePosition();
            const mfem::Vector &dist = finder.GetDist();

            std::ostringstream oss;
            oss << std::setprecision(17);

            oss << "[WARNING] FindPoints failed | x=(";
            for (int d = 0; d < sdim; ++d)
                oss << points[d] << (d + 1 < sdim ? "," : "");

            oss << ") | elem_ids=[";
            for (int i = 0; i < elem_ids.Size(); ++i)
                oss << elem_ids[i] << (i + 1 < elem_ids.Size() ? "," : "");
            oss << "] | codes=[";

            for (int i = 0; i < codes.Size(); ++i)
                oss << codes[i] << (i + 1 < codes.Size() ? "," : "");
            oss << "]";

            if (ref.Size() >= sdim)
            {
                oss << " | ref=(";
                for (int d = 0; d < sdim; ++d)
                    oss << ref[d] << (d + 1 < sdim ? "," : "");
                oss << ")";
            }

            if (gref.Size() >= sdim)
            {
                oss << " | gref=(";
                for (int d = 0; d < sdim; ++d)
                    oss << gref[d] << (d + 1 < sdim ? "," : "");
                oss << ")";
            }

            if (dist.Size() > 0) { oss << " | dist=" << dist[0]; }
            std::cout << oss.str() << std::endl;
            return;
        }

        finder.Interpolate(E_gf, E);

        double n2 = 0.0;
        for (int d = 0; d < sdim; ++d) n2 += E[d]*E[d];
        const double n = std::sqrt(n2);
        if (!(n > 0.0) || !std::isfinite(n))
        {
            std::fill(dxds.begin(), dxds.end(), 0.0);
            std::cout << "[WARNING] Field Magnitude ^2 is negative or non finite" << std::endl;
            return;
        }

        for (int d = 0; d < sdim; ++d) dxds[d] = -E[d] / n;
    }

    // Evaluate E at a given physical point.
    // Returns true if the point was successfully located and interpolated.
    //
    // Optional outputs:
    //   ref_out        : MFEM reference coordinates (size sdim)
    //   gslib_ref_out  : GSLIB internal reference coords in [-1,1] (size sdim)
    //   dist_out       : distance-to-border (meaningful especially for code==1)
    bool EvaluateFieldAt(const state_t &x,
                         mfem::Vector &E_out,
                         double &norm_out,
                         int &elem_id_out,
                         int &code_out,
                         mfem::Vector *ref_out       = nullptr,
                         mfem::Vector *gslib_ref_out = nullptr,
                         double *dist_out            = nullptr)
    {
        MFEM_VERIFY((int)x.size() == sdim, "State dimension mismatch.");

        for (int d = 0; d < sdim; ++d) points[d] = x[d];

        finder.FindPoints(points, mfem::Ordering::byNODES);

        const auto &elem_ids = finder.GetElem();
        const auto &codes    = finder.GetCode();

        if (elem_ids.Size() != 1 || codes.Size() != 1)
        {
            elem_id_out = -1;
            code_out    = -1;
            norm_out    = 0.0;

            if (ref_out)       { ref_out->SetSize(sdim);       (*ref_out) = 0.0; }
            if (gslib_ref_out) { gslib_ref_out->SetSize(sdim); (*gslib_ref_out) = 0.0; }
            if (dist_out)      { *dist_out = 0.0; }
            return false;
        }

        elem_id_out = (int)elem_ids[0];
        code_out    = (int)codes[0];

        // Reference coordinates: MFEM reference space
        if (ref_out)
        {
            const mfem::Vector &R = finder.GetReferencePosition(); // length = npts*sdim
            ref_out->SetSize(sdim);
            if (R.Size() >= sdim)
            {
                for (int d = 0; d < sdim; ++d) { (*ref_out)[d] = R[d]; }
            }
            else
            {
                (*ref_out) = 0.0;
            }
        }

        // Reference coordinates: GSLIB internal [-1,1]
        if (gslib_ref_out)
        {
            const mfem::Vector &RG = finder.GetGSLIBReferencePosition(); // length = npts*sdim
            gslib_ref_out->SetSize(sdim);
            if (RG.Size() >= sdim)
            {
                for (int d = 0; d < sdim; ++d) { (*gslib_ref_out)[d] = RG[d]; }
            }
            else
            {
                (*gslib_ref_out) = 0.0;
            }
        }

        if (dist_out)
        {
            const mfem::Vector &D = finder.GetDist(); // length = npts
            *dist_out = (D.Size() > 0) ? D[0] : 0.0;
        }

        if (elem_id_out < 0 || code_out != 0)
        {
            norm_out = 0.0;
            return false;
        }

        finder.Interpolate(E_gf, E_out);

        double n2 = 0.0;
        for (int d = 0; d < sdim; ++d) { n2 += E_out[d] * E_out[d]; }

        norm_out = std::sqrt(n2);
        return (std::isfinite(norm_out) && norm_out >= 0.0);
    }
};


// Exception used to break out of ode integrate
struct TraceExitEvent {};

// Observer: records points and terminates on leaving domain.
struct TraceObserver
{
    ElectronTraceResult      &out;
    const TpcGeometry        &geom;
    const ElectronTraceParams &p;
    const bool               save;

    TraceObserver(ElectronTraceResult &o,
                  const TpcGeometry &g,
                  const ElectronTraceParams &pp,
                  bool sp)
        : out(o), geom(g), p(pp), save(sp) {}

    void operator()(const state_t &x, double /*s*/)
    {
        MaybePush(out, x, save);
        if (!geom.Inside(x[0], x[1]))
        {
            out.exit_code = ClassifyExit(geom, x, p.geom_tol);
            // Need the final position 
            if (not save) {MaybePush(out, x, true);}
            throw TraceExitEvent{};
        }
    }
};

// Fixed-step integration using integrate_n_steps (no adaptivity).
template <class Stepper>
static inline ElectronTraceResult RunFixedNSteps(Stepper stepper,
                                                 FieldRHS &rhs,
                                                 const TpcGeometry &geom,
                                                 const ElectronTraceParams &p,
                                                 const double ds,
                                                 const long long max_steps, // 0 => unlimited
                                                 const bool save_pathlines,
                                                 state_t x)
{
    ElectronTraceResult out;
    out.exit_code = ElectronExitCode::None;

    if (!geom.Inside(x[0], x[1]))
    {
        out.exit_code = ClassifyExit(geom, x, p.geom_tol);
        MaybePush(out, x, true);
        return out;
    }

    // Push initial point
    MaybePush(out, x, save_pathlines);

    const std::size_t N =
        (max_steps > 0) ? static_cast<std::size_t>(max_steps)
                        : static_cast<std::size_t>(std::numeric_limits<std::size_t>::max());

    TraceObserver obs(out, geom, p, save_pathlines);

    try
    {
        // x is updated in-place by odeint.
        boost::numeric::odeint::integrate_n_steps(stepper, rhs, x, 0.0, ds, N, obs);
        // If we get here, we did all N steps without leaving the domain.
        out.exit_code = ElectronExitCode::MaxSteps;
        MaybePush(out, x, true);
    }
    catch (const TraceExitEvent &)
    {
        // out.exit_code already set by observer
    }

    return out;
}
template <class ControlledStepper>
static inline ElectronTraceResult RunAdaptive(
    ControlledStepper controlled_stepper,
    FieldRHS &rhs,
    const TpcGeometry &geom,
    const ElectronTraceParams &p,
    double ds_init,
    long long max_steps,
    bool save_pathlines,
    state_t x)
{
    ElectronTraceResult out;
    out.exit_code = ElectronExitCode::None;

    if (!geom.Inside(x[0], x[1]))
    {
        out.exit_code = ClassifyExit(geom, x, p.geom_tol);
        MaybePush(out, x, save_pathlines);
        return out;
    }

    MaybePush(out, x, save_pathlines);

    TraceObserver obs(out, geom, p, save_pathlines);

    double s = 0.0;
    double ds = ds_init;
    long long steps = 0;
    try
    {
        while (max_steps <= 0 || steps < max_steps)
        {
            // On success x, s and ds are modified
            // On failure only ds
            const boost::numeric::odeint::controlled_step_result res = 
                                    controlled_stepper.try_step(rhs, x, s, ds);
            if (res == boost::numeric::odeint::success)
            {
                s += ds;
                ++steps;
                obs(x, s);
            }
            else if (!(ds > 0.0) || !std::isfinite(ds))
            {
                out.exit_code = ElectronExitCode::DegenerateTimeStep;
                break;
            }
        }

        out.exit_code = ElectronExitCode::MaxSteps;
    }
    catch (const TraceExitEvent &)
    {
        // exit_code already set
    }

    return out;
}


ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh                    &mesh,
    mfem::FindPointsGSLIB            &finder,
    const mfem::ParGridFunction      &E_gf,
    const TpcGeometry                &geom,
    const ElectronTraceParams        &params,
    const mfem::Vector               &x0_in,          // (r,z,*) or (x,y,z)
    bool                              axisymmetric,
    bool                              save_pathlines,
    int                               max_traversals,
    mfem::Vector                     &pos_scratch,
    mfem::Vector                     &E_scratch,
    const double                      h_ref)
{
    using namespace boost::numeric::odeint;

    const int sdim = mesh.SpaceDimension();
    state_t x(sdim, 0.0);
    for (int d = 0; d < sdim; ++d){ x[d] = x0_in[d]; }

    ClampAxis(x, params.geom_tol);

    const double ds = params.c_step * h_ref;
    if (!(ds > 0.0))
        throw std::runtime_error("TraceSingleElectronLine: ds <= 0 (check c_step and mesh length scale).");

    FieldRHS rhs(finder, E_gf, sdim, params.geom_tol);

    const long long max_steps = ComputeMaxStepsFromLimit(params, ds, max_traversals);

    const std::string &m = params.method;

    if (m == "Euler-Cauchy")
        return RunFixedNSteps(euler<state_t>{}, rhs, geom, params, ds, max_steps, save_pathlines, x);

    else if (m == "RK4")
        return RunFixedNSteps(runge_kutta4<state_t>{}, rhs, geom, params, ds, max_steps, save_pathlines, x);

    else if (m == "RK45")
        return RunAdaptive(make_controlled(1e-6, 1e-6, runge_kutta_dopri5<state_t>{}), rhs, geom, params, ds, max_steps, save_pathlines, x);

    else if (m == "RK54")
        return RunAdaptive(make_controlled(1e-6, 1e-6, runge_kutta_cash_karp54<state_t>{}), rhs, geom, params, ds, max_steps, save_pathlines, x);
    throw std::invalid_argument("Integration method '" + m + "' unknown.");
}
// -----------------------
// Debug Helpers
// -----------------------

static void PrintExitConditionSummary(const Config                          &cfg,
                                      const Seeds                           &seeds,
                                      const std::vector<ElectronTraceResult> &out_results)
{
    if (!cfg.debug.debug) { return; }

    // Only rank 0 prints. Other ranks return safely even if out_results is empty.
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, size = 1;
#ifdef MFEM_USE_MPI
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
#endif

    // Local counts/volumes (ranks with empty out_results contribute 0)
    unsigned long long lc_hit_lgi     = 0;
    unsigned long long lc_hit_cathode = 0;
    unsigned long long lc_hit_wall    = 0;
    unsigned long long lc_left_volume = 0;
    unsigned long long lc_max_steps   = 0;
    unsigned long long lc_deg_dt      = 0;
    unsigned long long lc_none        = 0;
    unsigned long long lc_hit_axis    = 0;

    double lV_total       = 0.0;
    double lV_hit_lgi     = 0.0;
    double lV_hit_cathode = 0.0;
    double lV_hit_wall    = 0.0;
    double lV_left_volume = 0.0;
    double lV_max_steps   = 0.0;
    double lV_deg_dt      = 0.0;
    double lV_none        = 0.0;
    double lV_hit_axis    = 0.0;

    // Use the safe overlap length (handles ranks where out_results is empty or partial)
    const std::size_t n_res   = out_results.size();
    const std::size_t n_seedV = seeds.volumes.size();
    const std::size_t n_loc   = std::min(n_res, n_seedV);

    for (std::size_t i = 0; i < n_loc; ++i)
    {
        const auto &res = out_results[i];

        // Counts
        switch (res.exit_code)
        {
            case ElectronExitCode::HitLiquidGas:       lc_hit_lgi++;     break;
            case ElectronExitCode::HitCathode:         lc_hit_cathode++; break;
            case ElectronExitCode::HitWall:            lc_hit_wall++;    break;
            case ElectronExitCode::LeftVolume:         lc_left_volume++; break;
            case ElectronExitCode::MaxSteps:           lc_max_steps++;   break;
            case ElectronExitCode::DegenerateTimeStep: lc_deg_dt++;      break;
            case ElectronExitCode::HitAxis:            lc_hit_axis++;    break;
            case ElectronExitCode::None:               lc_none++;        break;
        }

        // Volume fractions (only if we have volumes for this index)
        const double dV = seeds.volumes[i];
        lV_total += dV;

        switch (res.exit_code)
        {
            case ElectronExitCode::HitLiquidGas:       lV_hit_lgi     += dV; break;
            case ElectronExitCode::HitCathode:         lV_hit_cathode += dV; break;
            case ElectronExitCode::HitWall:            lV_hit_wall    += dV; break;
            case ElectronExitCode::LeftVolume:         lV_left_volume += dV; break;
            case ElectronExitCode::MaxSteps:           lV_max_steps   += dV; break;
            case ElectronExitCode::DegenerateTimeStep: lV_deg_dt      += dV; break;
            case ElectronExitCode::HitAxis:            lV_hit_axis    += dV; break;
            case ElectronExitCode::None:               lV_none        += dV; break;
        }
    }

    // Reduce to rank 0
    unsigned long long gc_hit_lgi     = lc_hit_lgi;
    unsigned long long gc_hit_cathode = lc_hit_cathode;
    unsigned long long gc_hit_wall    = lc_hit_wall;
    unsigned long long gc_left_volume = lc_left_volume;
    unsigned long long gc_max_steps   = lc_max_steps;
    unsigned long long gc_deg_dt      = lc_deg_dt;
    unsigned long long gc_none        = lc_none;
    unsigned long long gc_hit_axis    = lc_hit_axis;

    double gV_total       = lV_total;
    double gV_hit_lgi     = lV_hit_lgi;
    double gV_hit_cathode = lV_hit_cathode;
    double gV_hit_wall    = lV_hit_wall;
    double gV_left_volume = lV_left_volume;
    double gV_max_steps   = lV_max_steps;
    double gV_deg_dt      = lV_deg_dt;
    double gV_none        = lV_none;
    double gV_hit_axis    = lV_hit_axis;

#ifdef MFEM_USE_MPI
    if (size > 1)
    {
        MPI_Reduce(&lc_hit_lgi,     &gc_hit_lgi,     1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_hit_cathode, &gc_hit_cathode, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_hit_wall,    &gc_hit_wall,    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_left_volume, &gc_left_volume, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_max_steps,   &gc_max_steps,   1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_deg_dt,      &gc_deg_dt,      1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_none,        &gc_none,        1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&lc_hit_axis,    &gc_hit_axis,    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);

        MPI_Reduce(&lV_total,       &gV_total,       1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_hit_lgi,     &gV_hit_lgi,     1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_hit_cathode, &gV_hit_cathode, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_hit_wall,    &gV_hit_wall,    1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_left_volume, &gV_left_volume, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_max_steps,   &gV_max_steps,   1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_deg_dt,      &gV_deg_dt,      1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_none,        &gV_none,        1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&lV_hit_axis,    &gV_hit_axis,    1, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
#endif

    if (rank != 0) { return; }

    const double N = static_cast<double>(gc_hit_lgi + gc_hit_cathode + gc_hit_wall +
                                         gc_hit_axis + gc_left_volume + gc_max_steps +
                                         gc_deg_dt + gc_none);

    auto frac = [&](unsigned long long c) -> double
    {
        return (N > 0.0 ? static_cast<double>(c) / N : 0.0);
    };

    auto vfrac = [&](double V) -> double
    {
        return (gV_total > 0.0 ? V / gV_total : 0.0);
    };

    std::cout << "\n---------------- Electron Tracing Summary ----------------\n";
    std::cout << "Total seeds traced: " << static_cast<unsigned long long>(N) << "\n\n";

    std::cout << std::left
              << std::setw(22) << "Exit Type"
              << std::setw(18) << "Seed Fraction"
              << std::setw(18) << "Volume Fraction\n";

    std::cout << "-----------------------------------------------------------------\n";

    auto row = [&](const char *label, double f_seed, double f_vol)
    {
        std::cout << std::left
                  << std::setw(22) << label
                  << std::setw(18) << f_seed
                  << std::setw(18) << f_vol
                  << "\n";
    };

    row("Hit Liquid-Gas",    frac(gc_hit_lgi),     vfrac(gV_hit_lgi));
    row("Hit Cathode",       frac(gc_hit_cathode), vfrac(gV_hit_cathode));
    row("Hit Wall",          frac(gc_hit_wall),    vfrac(gV_hit_wall));
    row("Hit Axis",          frac(gc_hit_axis),    vfrac(gV_hit_axis));
    row("Left Volume",       frac(gc_left_volume), vfrac(gV_left_volume));
    row("Max Steps",         frac(gc_max_steps),   vfrac(gV_max_steps));
    row("Degenerate dt",     frac(gc_deg_dt),      vfrac(gV_deg_dt));
    row("None (unexpected)", frac(gc_none),        vfrac(gV_none));

    std::cout << "-----------------------------------------------------------------\n";
}

// We need a POD for MPI to share data (we dont own mfem::vector)
struct ElectronTraceSummary
{
    int32_t exit_code;
    double  x[3];   // store up to 3D; for 2D use x[2]=0
};
static_assert(std::is_trivially_copyable<ElectronTraceSummary>::value,
              "Summary must be POD for MPI");
// Generate Seeds for each rank
static inline void SeedRangeForRank(std::size_t n, int rank, int size,
                                    std::size_t &begin, std::size_t &end)
{
    const std::size_t base = n / (std::size_t)size;
    const std::size_t rem  = n % (std::size_t)size;

    begin = (std::size_t)rank * base + std::min<std::size_t>((std::size_t)rank, rem);
    end   = begin + base + ((std::size_t)rank < rem ? 1 : 0);
}

static void TraceElectronFieldLinesInnerPreparedBOOST(
    mfem::ParMesh                   &mesh,
    mfem::FindPointsGSLIB           &finder,      // cached
    const mfem::ParGridFunction     &E_gf,
    const ElectronTraceParams       &params,
    const Config                    &cfg,
    const Seeds                     &seeds,
    std::vector<ElectronTraceResult> &out_results,
    bool axisymmetric,
    bool save_paths,
    const double *z_max_overrides,              // size=n_seeds or nullptr
    mfem::Vector                    &pos,        // cached scratch
    mfem::Vector                    &Eout,       // cached scratch
    double                           h_ref)       // cached
{
    using namespace mfem;

    const std::size_t n_seeds = seeds.positions.size();
    out_results.resize(n_seeds);

    const int sdim = mesh.SpaceDimension();
    MFEM_VERIFY(E_gf.VectorDim() == sdim,
                "E_gf VectorDim must match mesh.SpaceDimension()");
    MFEM_VERIFY(sdim == E_gf.ParFESpace()->GetMesh()->SpaceDimension(),
                "ODE state dimension must equal mesh SpaceDimension()");

    // --- MPI info ---
    int rank = 0, size = 1;
    MPI_Comm comm = mesh.GetComm();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Ensure scratch sizes are correct (in case mesh changed).
    if (pos.Size()  != sdim) pos.SetSize(sdim);
    if (Eout.Size() != sdim) Eout.SetSize(sdim);

    // 1) Compute local seed range [ib, ie)
    std::size_t ib = 0, ie = 0;
    SeedRangeForRank(n_seeds, rank, size, ib, ie);
    const std::size_t n_local = ie - ib;

    // 2) Compute local summaries
    std::vector<ElectronTraceSummary> local_sum;
    local_sum.resize(n_local);

    for (std::size_t ui = ib; ui < ie; ++ui)
    {
        ElectronTraceParams local_params = params;
        if (z_max_overrides) { local_params.z_max = z_max_overrides[ui]; }
        TpcGeometry local_geom(local_params);

        Vector x0(sdim);
        for (int d = 0; d < sdim; ++d) { x0[d] = seeds.positions[ui][d]; }

        ElectronTraceResult r = TraceSingleElectronLine(mesh, finder, E_gf,
                                                        local_geom, local_params,
                                                        x0, axisymmetric,
                                                        /*save_pathlines=*/false,
                                                        cfg.tracing_params.max_traversals,
                                                        pos, Eout, h_ref);

        ElectronTraceSummary s;
        s.exit_code = (int32_t)r.exit_code;
        // TODO Remove when ready
        MFEM_VERIFY(!r.points.empty(), "Expected final point stored when save_paths=false.");

        const mfem::Vector &last = r.points.back();
        s.x[0] = last.Size() > 0 ? last[0] : 0.0;
        s.x[1] = last.Size() > 1 ? last[1] : 0.0;
        s.x[2] = (sdim == 3 && last.Size() > 2) ? last[2] : 0.0;

        local_sum[ui - ib] = s;
    }
    // Populate for single MPI Node
    if (size == 1)
    {
        MFEM_VERIFY(local_sum.size() == n_seeds,
                    "size==1 but local_sum does not cover all seeds");

        out_results.resize(n_seeds);

        for (std::size_t i = 0; i < n_seeds; ++i)
        {
            const ElectronTraceSummary &s = local_sum[i];

            out_results[i].exit_code = static_cast<ElectronExitCode>(s.exit_code);
            out_results[i].points.clear();

            mfem::Vector v(sdim);
            v[0] = s.x[0];
            v[1] = s.x[1];
            if (sdim == 3) { v[2] = s.x[2]; }

            // single-point path (final position)
            out_results[i].points.push_back(std::move(v));
        }

        return;
    }
    // For more than one MPI Node

    // Gather counts (in number of ElectronTraceSummary elements)
    int sendcount = static_cast<int>(local_sum.size());

    std::vector<int> recvcounts;
    std::vector<unsigned long long> recv_ib; // use unsigned long long for portability with MPI_UNSIGNED_LONG_LONG
    if (rank == 0)
    {
        recvcounts.resize(static_cast<std::size_t>(size));
        recv_ib.resize(static_cast<std::size_t>(size));
    }

    // Gather sendcounts
    MPI_Gather(&sendcount, 1, MPI_INT,
               rank == 0 ? recvcounts.data() : nullptr, 1, MPI_INT,
               0, comm);

    // Gather each rank's starting global index (ib)
    unsigned long long ib_ull = static_cast<unsigned long long>(ib);
    MPI_Gather(&ib_ull, 1, MPI_UNSIGNED_LONG_LONG,
               rank == 0 ? recv_ib.data() : nullptr, 1, MPI_UNSIGNED_LONG_LONG,
               0, comm);

    // Build byte counts/displacements for MPI_Gatherv (in rank-order receive buffer)
    std::vector<int> recvcountsB, displsB;
    std::vector<char> all_bytes;

    if (rank == 0)
    {
        recvcountsB.resize(static_cast<std::size_t>(size));
        displsB.resize(static_cast<std::size_t>(size));

        int total_elems = 0;
        for (int r = 0; r < size; ++r) { total_elems += recvcounts[static_cast<std::size_t>(r)]; }
        MFEM_VERIFY(static_cast<std::size_t>(total_elems) == n_seeds,
                    "Total gathered summaries != n_seeds");

        int offB = 0;
        for (int r = 0; r < size; ++r)
        {
            recvcountsB[static_cast<std::size_t>(r)] =
                recvcounts[static_cast<std::size_t>(r)] * static_cast<int>(sizeof(ElectronTraceSummary));
            displsB[static_cast<std::size_t>(r)] = offB;
            offB += recvcountsB[static_cast<std::size_t>(r)];
        }

        all_bytes.resize(n_seeds * sizeof(ElectronTraceSummary));
    }

    // Gather summaries as bytes (rank-ordered in all_bytes)
    MPI_Gatherv(local_sum.data(),
                sendcount * static_cast<int>(sizeof(ElectronTraceSummary)), MPI_BYTE,
                rank == 0 ? all_bytes.data() : nullptr,
                rank == 0 ? recvcountsB.data() : nullptr,
                rank == 0 ? displsB.data() : nullptr,
                MPI_BYTE, 0, comm);

    if (rank == 0)
    {
        out_results.resize(n_seeds);

        // Place each rank's received block into its global offset explicitly using recv_ib[]
        for (int r = 0; r < size; ++r)
        {
            const std::size_t r_count = static_cast<std::size_t>(recvcounts[static_cast<std::size_t>(r)]);
            const std::size_t r_ib    = static_cast<std::size_t>(recv_ib[static_cast<std::size_t>(r)]);

            MFEM_VERIFY(r_ib + r_count <= n_seeds, "Rank block out of bounds");

            const char *rank_base = all_bytes.data() + displsB[static_cast<std::size_t>(r)];
            const auto *S = reinterpret_cast<const ElectronTraceSummary*>(rank_base);

            for (std::size_t k = 0; k < r_count; ++k)
            {
                const std::size_t gi = r_ib + k;

                out_results[gi].exit_code = static_cast<ElectronExitCode>(S[k].exit_code);
                out_results[gi].points.clear();

                mfem::Vector v(sdim);
                for (int d = 0; d < sdim; ++d) { v[d] = S[k].x[d]; } // one-line loop over sdim

                out_results[gi].points.push_back(std::move(v));
            }
        }
    }
    else
    {
        out_results.clear();
    }
}

// ------------------------ VTK Tracing ----------------------------------------
// TODO Did I remove the callback?
static constexpr int VTK_REASON_LEFT_TPC_DOMAIN = 20001;
struct TerminateOutsideContext
{
    TpcGeometry geom;
    explicit TerminateOutsideContext(const ElectronTraceParams &tp) : geom(tp) {}
};

// Terminate when the supplied point is no longer inside the domain.
// IMPORTANT: This checks the point already in the streamline polyline.
static bool TerminateIfOutsideDomain(void *clientdata,
                                    vtkPoints *pts,
                                    vtkDataArray*,
                                    int idx)
{
    auto *ctx = static_cast<TerminateOutsideContext*>(clientdata);

    double p[3] = {0,0,0};
    pts->GetPoint(idx, p);

    const bool inside = ctx->geom.Inside(p[0], p[1]);

    if (idx <= 2) // don’t spam
    {
        std::cerr << "[CB] idx=" << idx
                  << " p=(" << p[0] << "," << p[1] << "," << p[2] << ")"
                  << " inside=" << inside << "\n";
    }

    return !inside;
}

static vtkIdType InsertVTKCellFromMFEMElement(
    vtkUnstructuredGrid* ug,
    const mfem::Mesh& mesh,
    int el_id)
{
    const mfem::Element* el = mesh.GetElement(el_id);
    const int nv = el->GetNVertices();
    const int* v = el->GetVertices();

    // VTK expects vtkIdType connectivity
    std::vector<vtkIdType> conn(nv);
    for (int i = 0; i < nv; ++i) conn[i] = static_cast<vtkIdType>(v[i]);

    switch (el->GetGeometryType())
    {
        case mfem::Geometry::SEGMENT:
            if (nv != 2) return -1;
            ug->InsertNextCell(VTK_LINE, 2, conn.data());
            break;

        case mfem::Geometry::TRIANGLE:
            if (nv != 3) return -1;
            ug->InsertNextCell(VTK_TRIANGLE, 3, conn.data());
            break;

        case mfem::Geometry::SQUARE:
            if (nv != 4) return -1;
            ug->InsertNextCell(VTK_QUAD, 4, conn.data());
            break;

        case mfem::Geometry::TETRAHEDRON:
            if (nv != 4) return -1;
            ug->InsertNextCell(VTK_TETRA, 4, conn.data());
            break;

        case mfem::Geometry::CUBE:
            if (nv != 8) return -1;
            ug->InsertNextCell(VTK_HEXAHEDRON, 8, conn.data());
            break;

        case mfem::Geometry::PRISM:
            // MFEM "PRISM" is a wedge (6 verts)
            if (nv != 6) return -1;
            ug->InsertNextCell(VTK_WEDGE, 6, conn.data());
            break;

        case mfem::Geometry::PYRAMID:
            if (nv != 5) return -1;
            ug->InsertNextCell(VTK_PYRAMID, 5, conn.data());
            break;

        default:
            // Add more mappings as needed (e.g. mfem::Geometry::...)
            return -1;
    }

    return ug->GetNumberOfCells() - 1;
}

static vtkSmartPointer<vtkUnstructuredGrid>
BuildVTKGridFromMFEM(mfem::Mesh &mesh,
                     const mfem::GridFunction &E_gf,
                     bool axisymmetric)
{
    const int dim = mesh.SpaceDimension();
    MFEM_VERIFY(dim == 2 || dim == 3, "Only 2D/3D supported.");

    // ---------------------------------------------------------------------
    // 0) Points = MFEM vertices (linear grid)
    // ---------------------------------------------------------------------
    auto pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(mesh.GetNV());

    for (int vi = 0; vi < mesh.GetNV(); ++vi)
    {
        const double *X = mesh.GetVertex(vi);
        double p[3] = {0.0, 0.0, 0.0};

        if (dim == 3)
        {
            p[0] = X[0]; p[1] = X[1]; p[2] = X[2];
        }
        else
        {
            p[0] = X[0];
            p[1] = X[1];
            p[2] = 0.0;
        }

        pts->SetPoint(vi, p);
    }

    auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->SetPoints(pts);
    ug->Allocate(mesh.GetNE());

    // ---------------------------------------------------------------------
    // 1) Cells (MFEM -> VTK ordering using mfem/vtk.hpp)
    // ---------------------------------------------------------------------
    for (int el_id = 0; el_id < mesh.GetNE(); ++el_id)
    {
        const mfem::Element *el = mesh.GetElement(el_id);
        const mfem::Geometry::Type geom = el->GetGeometryType();

        const int vtk_cell_type = mfem::VTKGeometry::Map[geom];
        MFEM_VERIFY(vtk_cell_type > 0, "Unsupported geometry for VTK linear map.");

        const int nv = el->GetNVertices();
        const int *v = el->GetVertices();

        std::vector<vtkIdType> conn(nv);

        const int *perm = mfem::VTKGeometry::VertexPermutation[geom];
        if (perm)
        {
            for (int i = 0; i < nv; ++i) { conn[i] = static_cast<vtkIdType>(v[perm[i]]); }
        }
        else if (geom == mfem::Geometry::PRISM)
        {
            for (int i = 0; i < nv; ++i) { conn[i] = static_cast<vtkIdType>(v[mfem::VTKGeometry::PrismMap[i]]); }
        }
        else
        {
            for (int i = 0; i < nv; ++i) { conn[i] = static_cast<vtkIdType>(v[i]); }
        }

        ug->InsertNextCell(vtk_cell_type, nv, conn.data());
    }

    // ---------------------------------------------------------------------
    // 2) Point-data vector array "E" (3 comps) sampled at vertices
    //    Fast adjacency lookup using Mesh::GetVertexToElementTable()
    // ---------------------------------------------------------------------
    auto E = vtkSmartPointer<vtkDoubleArray>::New();
    E->SetName("E");
    E->SetNumberOfComponents(3);
    E->SetNumberOfTuples(mesh.GetNV());

    std::unique_ptr<mfem::Table> v2e(mesh.GetVertexToElementTable());
    MFEM_VERIFY(v2e, "GetVertexToElementTable returned null.");

    mfem::Vector eval(dim);

    for (int vi = 0; vi < mesh.GetNV(); ++vi)
    {
        const int n_adj = v2e->RowSize(vi);
        MFEM_VERIFY(n_adj > 0, "Vertex has no adjacent elements.");

        const int adj_el = v2e->GetRow(vi)[0];

        mfem::ElementTransformation *T = mesh.GetElementTransformation(adj_el);

        const mfem::Geometry::Type base_geom = mesh.GetElementBaseGeometry(adj_el);
        const mfem::Element *el = mesh.GetElement(adj_el);

        int local_vid = -1;
        {
            const int *vv = el->GetVertices();
            for (int k = 0; k < el->GetNVertices(); ++k)
            {
                if (vv[k] == vi) { local_vid = k; break; }
            }
        }
        MFEM_VERIFY(local_vid >= 0, "Internal error locating local vertex id.");

        const mfem::IntegrationRule *ir = mfem::Geometries.GetVertices(base_geom);
        MFEM_VERIFY(ir, "Geometries.GetVertices returned null.");
        MFEM_VERIFY(local_vid < ir->GetNPoints(), "local_vid out of range for reference vertices.");

        const mfem::IntegrationPoint &ip = ir->IntPoint(local_vid);
        T->SetIntPoint(&ip);

        E_gf.GetVectorValue(*T, ip, eval);

        double tuple[3] = { eval[0], eval[1], (dim == 3) ? eval[2] : 0.0 };
        E->SetTuple(vi, tuple);
    }

    ug->GetPointData()->AddArray(E);
    ug->GetPointData()->SetActiveVectors("E");

    return ug;
}


struct VTKTraceConfig
{
    bool follow_negative_E = true;

    double max_propagation = 1.0;
    double initial_step    = 1e-3;
    double min_step        = 1e-6;
    double max_step        = 1e-2;
    vtkIdType max_steps    = 20000;

    bool record_tracks     = true; // <-- your requested option
};

// Helper: build VTK seed polydata from your seeds, dropping seeds outside cells.
static vtkSmartPointer<vtkPolyData> BuildSeedPolyData(
    vtkDataSet* dataset,
    const std::vector<mfem::Vector>& seed_pos,
    vtkStaticCellLocator* locator,
    std::vector<int>& kept_seed_ids,
    bool axisymmetric)
{
    auto seedPts  = vtkSmartPointer<vtkPoints>::New();
    auto seedVert = vtkSmartPointer<vtkCellArray>::New();

    kept_seed_ids.clear();
    kept_seed_ids.reserve(seed_pos.size());

    for (int i = 0; i < (int)seed_pos.size(); ++i)
    {
        const mfem::Vector& s = seed_pos[i];

        double p[3] = {0,0,0};

        if (s.Size() == 3)
        {
            p[0] = s[0]; p[1] = s[1]; p[2] = s[2];
        }
        else
        {
            // same embedding rule as mesh conversion
            if (axisymmetric)
            {
                p[0] = s[0]; // r
                p[1] = s[1]; // z
                p[2] = 0.0;
            }
            else
            {
                p[0] = s[0];
                p[1] = s[1];
                p[2] = 0.0;
            }
        }

        vtkIdType cellId = locator->FindCell(p);
        if (cellId < 0) continue;

        vtkIdType id = seedPts->InsertNextPoint(p);
        seedVert->InsertNextCell(1);
        seedVert->InsertCellPoint(id);

        kept_seed_ids.push_back(i);
    }

    auto seedPoly = vtkSmartPointer<vtkPolyData>::New();
    seedPoly->SetPoints(seedPts);
    seedPoly->SetVerts(seedVert);
    return seedPoly;
}
// -----------------------------------------------------------------------------
// VTK PREPARED: does NOT build ug/locator. Uses cached ug+locator.
// -----------------------------------------------------------------------------
static void TraceElectronFieldLinesVTKPrepared(
    vtkUnstructuredGrid* ug,                     // cached
    vtkStaticCellLocator* locator,               // cached
    const Seeds &seeds,
    std::vector<ElectronTraceResult>& out_results,
    const VTKTraceConfig& tcfg,
    bool axisymmetric,
    const Config &cfg)
{
    MFEM_VERIFY(ug != nullptr, "VTKPrepared: ug is null");
    MFEM_VERIFY(locator != nullptr, "VTKPrepared: locator is null");

    // Validate vector field
    vtkDataArray* eArr = ug->GetPointData()->GetArray("E");
    MFEM_VERIFY(eArr && eArr->GetNumberOfComponents() == 3,
                "VTK dataset missing point-vector array 'E' (3 comps).");

    // Seeds: build polydata each call (can be cached too if seeds don't change)
    std::vector<int> kept_seed_ids;
    auto seedPoly = BuildSeedPolyData(
        ug,
        seeds.positions,
        locator,
        kept_seed_ids,
        axisymmetric);

    const std::size_t n_requested = seeds.positions.size();
    const std::size_t n_kept      = kept_seed_ids.size();

    if (n_kept != n_requested)
    {
        throw std::logic_error("VTK: Not all requested seeds accepted, particle seeding is invalid");
    }

    out_results.assign(n_requested, ElectronTraceResult{});

    // Integrator choice from cfg.tracing_params.method
    vtkSmartPointer<vtkInitialValueProblemSolver> method;
    if (cfg.tracing_params.method == "RK4")
        method = vtkSmartPointer<vtkRungeKutta4>::New();
    else if (cfg.tracing_params.method == "Euler-Cauchy")
        method = vtkSmartPointer<vtkRungeKutta2>::New();
    else
        method = vtkSmartPointer<vtkRungeKutta4>::New();

    auto tracer = vtkSmartPointer<vtkStreamTracer>::New();
    tracer->SetIntegrator(method);

    tracer->SetInputData(ug);
    tracer->SetSourceData(seedPoly);

    tracer->SetComputeVorticity(false);
    tracer->SetIntegrationDirectionToBackward();

    tracer->SetMaximumPropagation(tcfg.max_propagation);
    tracer->SetInitialIntegrationStep(tcfg.initial_step);
    tracer->SetMinimumIntegrationStep(tcfg.min_step);
    tracer->SetMaximumIntegrationStep(tcfg.max_step);
    tracer->SetMaximumNumberOfSteps(tcfg.max_steps);

    tracer->SetInputArrayToProcess(
        0, 0, 0,
        vtkDataObject::FIELD_ASSOCIATION_POINTS,
        "E");

    tracer->Update();

    vtkPolyData* out = tracer->GetOutput();
    if (!out || out->GetNumberOfCells() == 0 || out->GetNumberOfPoints() == 0)
    {
        // no streamlines produced
        return;
    }

    const vtkIdType nLines = out->GetNumberOfCells();
    const vtkIdType nSeeds = static_cast<vtkIdType>(kept_seed_ids.size());
    const vtkIdType nUse   = std::min(nLines, nSeeds);

    for (vtkIdType li = 0; li < nUse; ++li)
    {
        const int seed_idx = kept_seed_ids[static_cast<int>(li)];

        vtkCell* cell = out->GetCell(li);
        if (!cell) continue;

        vtkIdList* ids = cell->GetPointIds();
        if (!ids || ids->GetNumberOfIds() == 0) continue;

        if (tcfg.record_tracks)
        {
            auto &pts = out_results[seed_idx].points;
            pts.clear();
            pts.reserve(static_cast<std::size_t>(ids->GetNumberOfIds()));

            for (vtkIdType k = 0; k < ids->GetNumberOfIds(); ++k)
            {
                double pk[3];
                out->GetPoint(ids->GetId(k), pk);

                mfem::Vector v(3);
                v[0] = pk[0];
                v[1] = pk[1];
                v[2] = pk[2];
                pts.push_back(std::move(v));
            }
        }
        else
        {
            // Store only last point
            vtkIdType lastPid = ids->GetId(ids->GetNumberOfIds() - 1);
            double pk[3];
            out->GetPoint(lastPid, pk);

            mfem::Vector v(3);
            v[0] = pk[0];
            v[1] = pk[1];
            v[2] = pk[2];

            out_results[seed_idx].points.clear();
            out_results[seed_idx].points.push_back(std::move(v));
            out_results[seed_idx].points.shrink_to_fit();
        }
    }

    // Exit classification from last point (same as your original)
    TpcGeometry geom(cfg.tracing_params);
    for (vtkIdType li = 0; li < nUse; ++li)
    {
        const int seed_idx = kept_seed_ids[static_cast<int>(li)];
        auto &res = out_results[seed_idx];

        if (res.points.empty())
        {
            res.exit_code = ElectronExitCode::None;
            continue;
        }

        const mfem::Vector &last = res.points.back();
        res.exit_code = geom.ClassifyBoundary(last[0], last[1], cfg.tracing_params.geom_tol);
    }
}
static void TraceElectronFieldLinesVTK(
    mfem::ParMesh& pmesh,
    const mfem::GridFunction& E_gf_serial,
    const Seeds &seeds,
    std::vector<ElectronTraceResult>& out_results,
    const VTKTraceConfig& tcfg,
    bool axisymmetric,
    Config cfg)
{
    auto ug = BuildVTKGridFromMFEM(pmesh, E_gf_serial, axisymmetric);

    auto locator = vtkSmartPointer<vtkStaticCellLocator>::New();
    locator->SetDataSet(ug);
    locator->BuildLocator();

    TraceElectronFieldLinesVTKPrepared(ug, locator, seeds, out_results, tcfg, axisymmetric, cfg);
}

// --------------------------- Tracing Object holders
// 2) Context build: identical mesh prep + FindPointsGSLIB setup as BOOST
void MPITraceContext::Build(mfem::ParMesh &mesh, bool debug)
{
    h_ref = ComputeGlobalMinEdgeLength(mesh, debug);

    // Ensure nodes exist (same as BOOST)
    if (mesh.GetNodes() == nullptr)
    {
        const int order = 1;
        const bool discont = false;
        const int space_dim = mesh.SpaceDimension();
        const int ordering = mfem::Ordering::byNODES;
        mesh.SetCurvature(order, discont, space_dim, ordering);
    }

    finder = std::make_unique<mfem::FindPointsGSLIB>(mesh.GetComm());
    finder->Setup(mesh);
    finder->SetDefaultInterpolationValue(std::numeric_limits<double>::quiet_NaN());
}
void BoostTraceContext::Build(mfem::ParMesh &mesh, bool debug)
{
    const int sdim = mesh.SpaceDimension();

    h_ref = ComputeGlobalMinEdgeLength(mesh, debug);

    // Ensure nodes exist
    if (mesh.GetNodes() == nullptr)
    {
        const int order = 1;
        const bool discont = false;
        const int space_dim = mesh.SpaceDimension();
        const int ordering = mfem::Ordering::byNODES;
        mesh.SetCurvature(order, discont, space_dim, ordering);
    }

    finder = std::make_unique<mfem::FindPointsGSLIB>(mesh.GetComm());
    finder->Setup(mesh);
    finder->SetDefaultInterpolationValue(std::numeric_limits<double>::quiet_NaN());

    pos.SetSize(sdim);
    Eout.SetSize(sdim);
}
void VTKTraceContext::Build(mfem::ParMesh &mesh,
                        const mfem::ParGridFunction &E_gf_serial,
                        bool axisymmetric,
                        bool debug)
{
    h_ref = ComputeGlobalMinEdgeLength(mesh, debug);

    grid = BuildVTKGridFromMFEM(mesh, E_gf_serial, axisymmetric);

    locator = vtkSmartPointer<vtkStaticCellLocator>::New();
    locator->SetDataSet(grid);
    locator->BuildLocator();
}
bool VTKTraceContext::Ready() const { return (grid != nullptr && locator != nullptr); }

// ElectronFieldLineTracer methods 
void ElectronFieldLineTracer::Reset()
{
    pmesh = nullptr;

    fec_vec.reset();
    fes_vec.reset();
    E_gf_owned.reset();
    E_gf = nullptr;

    boost.reset();
    vtk.reset();
    trace_fn = nullptr;

    mpitracer.reset();
}
void ElectronFieldLineTracer::Setup(const SimulationResult &result,
                                    const ElectronTraceParams &tp,
                                    const Config &cfg_,
                                    bool axisymmetric)
{
    Reset();

    // Attach mesh
    pmesh = result.mesh.get();
    MFEM_VERIFY(pmesh, "SimulationResult.mesh is null");

    // Attach potential
    auto *phi = result.V.get();
    MFEM_VERIFY(phi, "SimulationResult.V is null");

    const int dim = pmesh->SpaceDimension();

    const mfem::FiniteElementSpace *fes_phi = phi->FESpace();
    MFEM_VERIFY(fes_phi, "phi has no FESpace");
    const int order = fes_phi->GetMaxElementOrder();

    // Build vector FE space for E
    fec_vec = std::make_unique<mfem::H1_FECollection>(order, dim);
    fes_vec = std::make_unique<mfem::ParFiniteElementSpace>(
        pmesh, fec_vec.get(), dim, mfem::Ordering::byVDIM);

    E_gf_owned = std::make_unique<mfem::ParGridFunction>(fes_vec.get());
    *E_gf_owned = 0.0;

    ElectricFieldCoeff E_coeff(*phi, -1.0);
    E_gf_owned->ProjectCoefficient(E_coeff);

    E_gf = E_gf_owned.get();

    params = tp;
    cfg    = cfg_;
    if (params.provider == "MPITracer")
    {
        // MPITracer is intended for multi-rank; allow size==1 too (for debugging).
        BuildMPITracer_(cfg.debug.debug);
        SelectProvider_(params.provider, axisymmetric);
    }
    else if (params.provider == "BOOST")
    {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) // TODO All MPI SPecifics in this method gotta go 
        {
            BuildBOOST_(cfg.debug.debug);
            SelectProvider_(params.provider, axisymmetric);
        }
    }
    else if (params.provider == "VTK")
    {
        int comm_size = 1;
        MPI_Comm_size(pmesh->GetComm(), &comm_size);
        if (comm_size > 1)
        {
            throw std::logic_error("VTK tracing does not support multiple MPI ranks.");
        }
        BuildVTK_(axisymmetric, cfg.debug.debug);
        SelectProvider_(params.provider, axisymmetric);
    }
    else
    {
        throw std::invalid_argument("Tracing provider '" + params.provider + "' not implemented.");
    }
}
void ElectronFieldLineTracer::Trace(const Seeds &seeds,
            std::vector<ElectronTraceResult> &out_results,
            bool axisymmetric,
            bool save_paths,
            const double *z_max_overrides) const
{
    MFEM_VERIFY((bool)trace_fn, "SelectProvider(...) must be called before Trace().");
    trace_fn(seeds, out_results, axisymmetric, save_paths, z_max_overrides);
    // Misc things
    PrintExitConditionSummary(cfg, seeds, out_results);
    DumpElectronPathsCSV(cfg, out_results);
    ValidateElectronExitCodesAllSeeds(out_results);
}
void ElectronFieldLineTracer::BuildBOOST_(bool debug)
{
    if (!pmesh || !E_gf)
    {
        throw std::logic_error("BuildBOOST_: mesh/E_gf not attached.");
    }
    boost.emplace();
    boost->Build(*pmesh, debug);
}
void ElectronFieldLineTracer::BuildVTK_(bool axisymmetric, bool debug)
{
    vtk.emplace();
    vtk->Build(*pmesh, *E_gf, axisymmetric, debug);
}
void ElectronFieldLineTracer::BuildMPITracer_(bool debug)
{
    if (!pmesh || !E_gf)
    {
        throw std::logic_error("BuildMPITracer_: mesh/E_gf not attached.");
    }
    mpitracer.emplace();
    mpitracer->Build(*pmesh, debug);
}
void ElectronFieldLineTracer::SelectProvider_(const std::string &provider, bool axisymmetric)
{   
    if (provider == "MPITracer")
    {
        if (!mpitracer || !mpitracer->Ready())
        {
            throw std::logic_error("SelectProvider_: MPITracer context not ready.");
        }

        trace_fn = [this](const Seeds &seeds,
                          std::vector<ElectronTraceResult> &out_results,
                          bool axisymmetric,
                          bool save_paths,
                          const double * /*z_max_overrides*/)
        {
            // per your current constraint: no save_paths in multi-rank (enforced earlier),
            // but keep a local guard for safety.
            int comm_size = 1;
            MPI_Comm_size(pmesh->GetComm(), &comm_size);

            // z_max_overrides not supported in this provider yet (miniapp-style particles)
            // If you need it later, incorporate per-particle z_max into tags/fields.

            mpitracing::TraceDistributedEuler(
                *pmesh,
                *mpitracer->finder,
                *E_gf,
                params,
                seeds,
                out_results,
                mpitracer->h_ref,
                axisymmetric,
                /*redistribution_every=*/cfg.tracing_params.redistribution_every, // add to Config
                /*debug=*/cfg.debug.debug,
                /*debug_every=*/0,
                cfg.debug.dumpdata);
        };
    }
    else if (provider == "BOOST")
    {
        if (!boost || !boost->finder)
            throw std::logic_error("SelectProvider_: BOOST context not ready.");

        trace_fn = [this](const Seeds &seeds,
                                        std::vector<ElectronTraceResult> &out_results,
                                        bool axisymmetric,
                                        bool save_paths,
                                        const double *z_max_overrides)
        { TraceElectronFieldLinesInnerPreparedBOOST(*pmesh,
                                                    *boost->finder,
                                                    *E_gf,
                                                    params,
                                                    cfg,
                                                    seeds,
                                                    out_results,
                                                    axisymmetric,
                                                    save_paths,
                                                    z_max_overrides,
                                                    boost->pos,
                                                    boost->Eout,
                                                    boost->h_ref); };
    }
    else if (provider == "VTK")
    {
        if (!vtk || !vtk->Ready()) { throw std::logic_error("SelectProvider_: VTK context not ready."); }

        trace_fn = [this](const Seeds &seeds,
                            std::vector<ElectronTraceResult> &out_results,
                            bool axisymmetric,
                            bool save_paths,
                            const double * /*z_max_overrides*/)
        {
            VTKTraceConfig tcfg;
            tcfg.follow_negative_E = true;

            tcfg.max_propagation   = params.max_traversals;
            tcfg.initial_step      = params.c_step * vtk->h_ref;
            tcfg.min_step          = 0.1 * tcfg.initial_step;
            tcfg.max_step          = 10.0 * tcfg.initial_step;
            tcfg.max_steps         = static_cast<vtkIdType>(params.max_traversals / (0.001 * tcfg.initial_step));
            tcfg.record_tracks     = save_paths;

            TraceElectronFieldLinesVTKPrepared(vtk->grid, vtk->locator,
                                            seeds, out_results, tcfg, axisymmetric, cfg);
        };

    }
    else { throw std::invalid_argument("SelectProvider_: unknown provider '" + provider + "'."); }
}

