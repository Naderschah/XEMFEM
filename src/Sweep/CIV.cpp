#include <mfem.hpp>
#include <vector>
#include <cmath>
#include <omp.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "CIV.h"


// Return true for any form of charge insensitivity (cathode or wall)
static bool IsChargeInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitWall);
}

// Build volume-weighted seeds for CIV / optimization / tracing
// based on the mesh and cfg.tracing_params.
static CivSeeds ExtractCivSeeds(const Config &cfg, const SimulationResult &result)
{
    mfem::ParMesh &pmesh = *result.mesh;
    const int dim  = pmesh.Dimension();

    TpcGeometry geom(cfg.tracing_params);

    const int ir_order = cfg.civ_params.ir_order;

    mfem::IntegrationRules irs;
    CivSeeds seeds;

    const int ne = pmesh.GetNE();
    seeds.positions.reserve(ne);
    seeds.elements.reserve(ne);
    seeds.ips.reserve(ne);

    std::vector<double> element_volume(ne, 0.0);
    std::vector<int>    seeds_per_element(ne, 0);

    // ------------------ STEP 1: find geometry-eligible elements ------------------
    std::vector<int> eligible_elems;
    eligible_elems.reserve(ne);

    for (int e = 0; e < ne; ++e)
    {
        mfem::ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "ExtractCivSeeds: ElementTransformation is null");

        mfem::Geometry::Type geom_type  = T->GetGeometryType();
        const mfem::IntegrationRule &ir = irs.Get(geom_type, ir_order);
        const int n_ip                  = ir.GetNPoints();

        bool any_inside = false;

        for (int i = 0; i < n_ip; ++i)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            mfem::Vector x(dim);
            T->Transform(ip, x);
            const double r = x(0);
            const double z = x(1);

            if (geom.Inside(r, z))
            {
                any_inside = true;
                break;
            }
        }

        if (any_inside)
        {
            eligible_elems.push_back(e);
        }
    }
    // ---------------------------------------------------------------------------

    // ------------------ STEP 2: randomly select from eligible elements ---------
    int num_seed_elements = cfg.civ_params.num_seed_elements;
    const int n_eligible  = static_cast<int>(eligible_elems.size());

    if (num_seed_elements <= 0 || num_seed_elements > n_eligible)
    {
        num_seed_elements = n_eligible; // fall back to all eligible
    }

    std::mt19937 rng;
    if (cfg.civ_params.rng_seed != 0)
    {
        rng.seed(static_cast<std::mt19937::result_type>(cfg.civ_params.rng_seed));
    }
    else
    {
        std::random_device rd;
        rng.seed(rd());
    }

    std::vector<int> selected_elems;
    selected_elems.reserve(num_seed_elements);

    std::sample(eligible_elems.begin(), eligible_elems.end(),
                std::back_inserter(selected_elems),
                num_seed_elements,
                rng);
    // ---------------------------------------------------------------------------

    // ------------------ STEP 3: extract seeds only from selected elements ------
    for (int idx = 0; idx < num_seed_elements; ++idx)
    {
        int e = selected_elems[idx];

        mfem::ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "ExtractCivSeeds: ElementTransformation is null");

        mfem::Geometry::Type geom_type  = T->GetGeometryType();
        const mfem::IntegrationRule &ir = irs.Get(geom_type, ir_order);
        const int n_ip                  = ir.GetNPoints();

        double V_e = 0.0;

        for (int i = 0; i < n_ip; ++i)
        {
            const mfem::IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            mfem::Vector x(dim);
            T->Transform(ip, x);
            const double r = x(0);
            const double z = x(1);

            // Contribution to element volume (axisymmetric without 2π)
            const double detJ  = T->Weight();   // |det(J)| in (r,z)
            const double w_ref = ip.weight;     // reference weight
            const double dV    = r * detJ * w_ref;

            V_e += dV;

            if (!geom.Inside(r, z))
            {
                continue;
            }

            seeds.positions.push_back(x);
            seeds.elements.push_back(e);

            mfem::IntegrationPoint ip_copy = ip;
            seeds.ips.push_back(ip_copy);

            seeds_per_element[e]++;
        }

        element_volume[e] = V_e;
    }
    // ---------------------------------------------------------------------------

    // ------------------ STEP 4: same volume partition as before ----------------
    const int n_seeds = static_cast<int>(seeds.positions.size());
    seeds.volumes.resize(n_seeds);

    for (int j = 0; j < n_seeds; ++j)
    {
        const int e         = seeds.elements[j];
        const int n_in_elem = seeds_per_element[e];

        MFEM_VERIFY(n_in_elem > 0,
                    "ExtractCivSeeds: element has seeds reference but zero seeds_per_element");

        seeds.volumes[j] = element_volume[e] / static_cast<double>(n_in_elem);
    }

    return seeds;
}



double ComputeCIV_RandomSample(const Config           &cfg,
                               const SimulationResult &result)
{
    using namespace mfem;

    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (RandomSample)" << std::endl;
    }

    ParMesh &pmesh = *result.mesh;

    CivSeeds seeds = ExtractCivSeeds(cfg, result);

    // If no seeds throw error
    if (seeds.positions.empty())
    {
        throw std::runtime_error("No seeds to compute CIV for");
    }

    std::vector<ElectronTraceResult> trace_results;

    TraceElectronFieldLines(result, cfg, seeds, trace_results);

    double V_total = 0.0; // active TPC volume
    double V_civ   = 0.0; // Volume where electrons do not reach liquid gas interface

    for (std::size_t i = 0; i < seeds.positions.size(); ++i)
    {
        const double                dV  = seeds.volumes[i];
        const ElectronTraceResult  &res = trace_results[i];

        V_total += dV;

        // charge-insensitive if it does NOT reach liquid-gas interface
        if (IsChargeInsensitive(res))
        {
            V_civ += dV;
        }
    }

    if (V_total <= 0.0)
    {
        return 0.0;
    }

    // Return volume fraction
    return V_civ / V_total;
}


// ============================================================================
// CIV by binary search 
// ============================================================================
static bool FindElementForPointGlobal(
    mfem::ParMesh          &mesh,
    const mfem::Vector     &x_phys,
    int                    &out_elem,
    mfem::IntegrationPoint &out_ip)
{
    using namespace mfem;

    const int ne = mesh.GetNE();

    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = mesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "FindElementForPointGlobal: ElementTransformation is null");

        InverseElementTransformation invT(T);

        IntegrationPoint ip_trial;
        int res = invT.Transform(x_phys, ip_trial);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = e;
            out_ip   = ip_trial;
            return true;
        }
    }

    // No element contains x_phys
    return false;
}

// TODO: Whats the point of this again
struct TracerContext
{
    // Per-thread-owned MFEM objects
    mfem::ParMesh               mesh;        // local copy of global mesh
    mfem::ParFiniteElementSpace fes_phi;     // H1 space for potential
    mfem::GridFunction          phi;         // potential on local mesh
    ElectricFieldCoeff          E_coeff;     // E = -grad(phi)

    // Shared read-only data (references or pointers)
    const ElementAdjacency              &adj;

    // Physics / configuration
    TpcGeometry         tpc;
    ElectronTraceParams params;
    bool                axisymmetric;
    Config              cfg;

    TracerContext(mfem::ParMesh                     &global_mesh,
                  const mfem::FiniteElementCollection *fec_phi,
                  int                                  ordering,
                  const mfem::GridFunction           &global_phi,
                  const ElementAdjacency             &adj_ref,
                  const Config                       &cfg)
        : mesh(global_mesh),
          fes_phi(&mesh, fec_phi, /*vdim=*/1, ordering),
          phi(&fes_phi),
          E_coeff(phi, -1.0),
          adj(adj_ref),
          tpc(cfg.tracing_params),
          params(cfg.tracing_params),
          axisymmetric(cfg.solver.axisymmetric),
          cfg(cfg)
    {
        phi = global_phi; // copy potential DOFs into local mesh
    }

    TracerContext(const TracerContext&) = delete;
    TracerContext& operator=(const TracerContext&) = delete;
    TracerContext(TracerContext&&) = default;
    TracerContext& operator=(TracerContext&&) = default;
};
// Trace a single electron starting from physical (r,z).
// Returns true if the point was inside the mesh and we traced it.
// On success, fills 'res' and 'start_elem'.
static bool TraceSingleCIProbe(TracerContext       &ctx,
                               double               r,
                               double               z,
                               int                 &elem_guess,
                               ElectronTraceResult &res)
{
    using namespace mfem;

    // Reject points outside the TPC region
    if (!ctx.tpc.Inside(r, z))
    {
        return false;
    }

    // Physical coordinate vector
    const int dim = ctx.mesh.SpaceDimension();

    Vector x_phys(dim);
    x_phys[0] = r;
    x_phys[1] = z;

    int              elem = -1;
    IntegrationPoint ip;

    // 1) try local search from elem_guess if valid
    if (elem_guess >= 0 &&
        elem_guess < ctx.mesh.GetNE() &&
        FindElementForPointLocal(ctx.mesh, ctx.adj, elem_guess, x_phys, elem, ip))
    {
        elem_guess = elem;
    }
    // 2) otherwise global search
    else if (!FindElementForPointGlobal(ctx.mesh, x_phys, elem, ip))
    {
        return false; // not in mesh
    }
    else
    {
        elem_guess = elem;
    }

    // Now trace from (elem, ip) using existing RK4 integrator
    res = TraceSingleElectronLine(ctx.mesh,
                                  ctx.E_coeff,
                                  ctx.adj,
                                  ctx.tpc,
                                  ctx.params,
                                  elem,
                                  ip,
                                  ctx.axisymmetric,
                                  /*save_pathlines=*/false);

    return true;
}


// Used to do simple downwards integration of the resulting boundaries instead of 
// summing mesh elements
static double ComputeCIV_FromInterface(const Config                  &cfg,
                                       const std::vector<InterfacePoint> &iface)
{
    if (iface.empty())
    {
        return 0.0;
    }

    TpcGeometry tpc(cfg.tracing_params);

    const double r_min = tpc.r_min;
    const double r_max = tpc.r_max;
    const double z_min = tpc.z_min;

    // Use the same "top of active drift region" as in the old element-based
    // ComputeCIVFromInterface (z_check_max = cfg.tracing_params.z_max).
    const double z_top = cfg.tracing_params.z_max;  // or tpc.z_max if you prefer

    if (r_max <= r_min || z_top <= z_min)
    {
        return 0.0;
    }

    std::vector<InterfacePoint> pts = iface;

    auto clamp_z = [&](double z) -> double {
        if (z < z_min) return z_min;
        if (z > z_top) return z_top;
        return z;
    };

    double A_ci = 0.0; // CI area in (r,z)

    // Segment from r_min up to first interface point: use z_ci = pts.front().z_ci
    const double r0 = pts.front().r;
    if (r_min < r0)
    {
        double zc = clamp_z(pts.front().z_ci);
        double h  = zc - z_min;
        if (h > 0.0)
        {
            A_ci += h * (r0 - r_min);
        }
    }

    // Interior segments: piecewise linear between successive interface points
    const int n = static_cast<int>(pts.size());
    for (int i = 0; i + 1 < n; ++i)
    {
        const double ra = pts[i].r;
        const double rb = pts[i+1].r;
        if (rb <= ra) { continue; }

        double za = clamp_z(pts[i].z_ci);
        double zb = clamp_z(pts[i+1].z_ci);

        double ha = za - z_min;
        double hb = zb - z_min;

        // Clip negative heights (if z_ci dips below z_min numerically)
        if (ha < 0.0) ha = 0.0;
        if (hb < 0.0) hb = 0.0;

        if (ha == 0.0 && hb == 0.0) { continue; }

        // Trapezoidal area under a linear segment: (average height) * width
        const double dr    = rb - ra;
        const double h_avg = 0.5 * (ha + hb);
        A_ci += h_avg * dr;
    }

    // Segment from last interface point up to r_max: use z_ci = pts.back().z_ci
    const double rn = pts.back().r;
    if (rn < r_max)
    {
        double zc = clamp_z(pts.back().z_ci);
        double h  = zc - z_min;
        if (h > 0.0)
        {
            A_ci += h * (r_max - rn);
        }
    }

    // Total area of the active rectangle in (r,z)
    const double A_tot = (r_max - r_min) * (z_top - z_min);

    return A_ci / A_tot;
}
// ============================================================================
// CIV by radial column sweeps (fixed r-slices, z-search per slice)
// ============================================================================
double ComputeCIV_ColumnSweep(const SimulationResult &sim,
                              const Config           &cfg)
{
    using namespace mfem;

    MFEM_VERIFY(sim.mesh, "ComputeCIV_ColumnSweep: mesh is null");
    ParMesh &global_mesh = *sim.mesh;

    MFEM_VERIFY(sim.V, "ComputeCIV_ColumnSweep: potential V is null");
    GridFunction &global_phi = *sim.V;

    const FiniteElementSpace      *fes_phi = global_phi.FESpace();
    const FiniteElementCollection *fec     = fes_phi->FEColl();
    const int                      ordering = fes_phi->GetOrdering();

    // Shared geometry/connectivity
    ElementAdjacency              adj         = BuildAdjacency(global_mesh);

    // Base tracer context
    TracerContext base_ctx(global_mesh,
                           fec,
                           ordering,
                           global_phi,
                           adj,
                           cfg);
    // Speed up by not requiring to traverse the entire TPC 

    // Geometry limits and params
    TpcGeometry tpc(cfg.tracing_params);
    const double r_min    = tpc.r_min;
    const double r_max    = tpc.r_max;
    const double z_min    = tpc.z_min;
    const double z_max    = tpc.z_max;
    const double geom_tol = cfg.tracing_params.geom_tol;
    const bool   debug    = cfg.debug.debug;

    if (debug)
    {
        std::cout << "[DEBUG:CIV] r_min=" << r_min
                  << " r_max=" << r_max
                  << " z_min=" << z_min
                  << " z_max=" << z_max
                  << " geom_tol=" << geom_tol << "\n";
    }

    // Number of radial slices
    const int n_r_slices = cfg.civ_params.n_slices; 

    // Radial range we actually sample (avoid exact boundaries)
    const double r_lo = r_min + geom_tol;
    const double r_hi = r_max - geom_tol;
    if (r_hi <= r_lo)
    {
        if (debug)
        {
            std::cout << "[DEBUG:CIV] Invalid radial range in ColumnSweep: r_hi <= r_lo\n";
        }
        return 0.0;
    }

    // Max number of traces per column
    int n_threads = cfg.compute.threads.num;
    if (n_threads <= 0) { n_threads = 1; }
    const int max_probes_per_column = 50; // TOOD Config hook

    if (debug)
    {
        std::cout << "[DEBUG:CIV] n_slices=" << n_r_slices
                  << " threads=" << n_threads
                  << " max_probes_per_column=" << max_probes_per_column << "\n";
    }

    // Build per-thread contexts
    int n_used_threads = std::max(1, n_threads);
    std::vector<std::unique_ptr<TracerContext>> extra;
    extra.reserve(n_used_threads - 1);

    std::vector<TracerContext*> ctx_pool(n_used_threads);
    ctx_pool[0] = &base_ctx;

    for (int t = 1; t < n_used_threads; ++t)
    {
        extra.emplace_back(std::make_unique<TracerContext>(
            base_ctx.mesh, fec, ordering, base_ctx.phi,
            adj, cfg));
        ctx_pool[t] = extra.back().get();
    }

    // Storage for interface points (per slice) and flags
    std::vector<InterfacePoint> iface_tmp(n_r_slices);
    std::vector<bool>           col_found(n_r_slices, false);

    // Parallel over radial slices; each thread uses its own TracerContext
    #pragma omp parallel for num_threads(n_used_threads)
    for (int k = 0; k < n_r_slices; ++k)
    {
        int tid = omp_get_thread_num();
        if (tid >= n_used_threads) { tid %= n_used_threads; }

        TracerContext &ctx_k = *ctx_pool[tid];

        // Radial location for this column
        const double t   = (k + 0.5) / static_cast<double>(n_r_slices);
        const double r_k = r_lo + t * (r_hi - r_lo);

        // -----------------------------------------
        // REUSE elem_guess across all z probes in this column
        // -----------------------------------------
        int elem_guess = -1;

        // Bracket top/bottom
        double z_top = -1.5;//z_max - geom_tol; //TODO add config entry
        double z_bot = z_min + geom_tol; 

        // Top probe
        bool ok_top = false, ci_top = false;
        {
            double rr = r_k, zz = z_top;
            ElectronTraceResult res;
            ok_top = TraceSingleCIProbe(ctx_k, rr, zz, elem_guess, res);
            ci_top = ok_top && IsChargeInsensitive(res);
        }

        // Bottom probe (reuses elem_guess)
        bool ok_bot = false, ci_bot = false;
        {
            double rr = r_k, zz = z_bot;
            ElectronTraceResult res;
            ok_bot = TraceSingleCIProbe(ctx_k, rr, zz, elem_guess, res);
            ci_bot = ok_bot && IsChargeInsensitive(res);
        }

        bool   found_ci = false;
        double z_ci     = 0.0;

        if (!ok_bot || !ci_bot)
        {
            // No CI even at bottom -> this column contributes nothing
            found_ci = false;
        }
        else if (ci_top)
        {
            // CI already at top -> whole column CI
            found_ci = true;
            z_ci     = z_top;
        }
        else
        {
            // Bisection between z_top (non-CI) and z_bot (CI)
            double z_hi = z_top;
            double z_lo = z_bot;
            z_ci        = z_lo;

            int remaining = max_probes_per_column - 2;

            while (remaining > 0 && (z_hi - z_lo) > 5.0 * geom_tol)
            {
                double z_mid = 0.5 * (z_hi + z_lo);
                double rr = r_k, zz = z_mid;

                ElectronTraceResult res;
                bool ok_mid = TraceSingleCIProbe(ctx_k, rr, zz, elem_guess, res);
                bool ci_mid = ok_mid && IsChargeInsensitive(res);

                if (ci_mid)
                {
                    z_ci = z_mid;
                    z_lo = z_mid;   // CI side moves up
                }
                else
                {
                    z_hi = z_mid;   // non-CI side moves down
                }

                --remaining;
            }

            found_ci = true;
        }

        if (found_ci)
        {
            col_found[k] = true;
            InterfacePoint ip;
            ip.r          = r_k;
            ip.z_ci       = z_ci;
            ip.elem       = elem_guess;  // optional; not needed for volume
            ip.cathode_ci = false;
            ip.wall_ci    = false;
            iface_tmp[k]  = ip;
        }
        else
        {
            col_found[k] = false;
        }
    }    // end parallel over columns

    // Compact interface points to only those columns where we found CI
    std::vector<InterfacePoint> iface;
    iface.reserve(n_r_slices);
    for (int k = 0; k < n_r_slices; ++k)
    {
        if (col_found[k])
        {
            iface.push_back(iface_tmp[k]);
        }
    }

    if (iface.empty())
    {
        if (debug)
        {
            std::cout << "[DEBUG:CIV] ColumnSweep: no interface points found; CIV = 0.\n";
        }
        return 0.0;
    }

    // Sort by radius so EvaluateInterfaceZ works correctly
    std::sort(iface.begin(), iface.end(),
              [](const InterfacePoint &a, const InterfacePoint &b)
              {
                  return a.r < b.r;
              });

    if (cfg.civ_params.dump_civ_boundary)
    {
        namespace fs = std::filesystem;
        fs::path outdir(cfg.save_path);
        fs::path fname = outdir / "civ_interface_boundary_columnsweep.csv";
        std::ofstream ofs(fname);
        ofs << "# r, z_ci\n";
        for (const auto &ip : iface)
            ofs << ip.r << "," << ip.z_ci << "\n";

        if (debug)
        {
            std::cout << "[DEBUG:CIV] Dumped ColumnSweep CIV boundary to: "
                      << fname.string() << "\n";
        }
    }
    double V_CI = ComputeCIV_FromInterface(cfg, iface);

    double V_total = (r_max - r_min) * (z_max - z_min ); //vres.V_total;
    //double V_CI    = vres.V_CI;
    double frac    = (V_total > 0.0) ? (V_CI / V_total) : 0.0;

    if (debug)
    {
        std::cout << "\n----- Charge Insensitive Volume (ColumnSweep) -----\n";
        std::cout << "Total Volume considered: " << V_total << "\n";
        std::cout << "CIV Volume:              " << V_CI    << "\n";
        std::cout << "CIV Fraction:            " << frac    << "\n";
        std::cout << "---------------------------------------------------\n";
    }

    return frac;
}
// ============================================================================
// Informed Sweep 
// ============================================================================
struct ElementGeometryData
{
    std::vector<double> r_centroid;
    std::vector<double> z_centroid;
    std::vector<double> volume;
};
static ElementGeometryData ComputeElementGeometry(mfem::ParMesh &mesh)
{
    using namespace mfem;

    const int ne  = mesh.GetNE();
    const int dim = mesh.SpaceDimension();

    ElementGeometryData geom;
    geom.r_centroid.resize(ne);
    geom.z_centroid.resize(ne);
    geom.volume.resize(ne);

    IntegrationRules irs;

    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = mesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "ComputeElementGeometry: ElementTransformation is null");

        Geometry::Type geom_type = T->GetGeometryType();
        // Order 2 is enough for centroid and volume here
        const IntegrationRule &ir = irs.Get(geom_type, 2);
        const int n_ip            = ir.GetNPoints();

        double V_e   = 0.0;
        double r_acc = 0.0;
        double z_acc = 0.0;

        for (int i = 0; i < n_ip; ++i)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            Vector x(dim);
            T->Transform(ip, x);
            const double r = x(0);
            const double z = x(1);

            const double detJ  = T->Weight();
            const double w_ref = ip.weight;
            const double dV    = r * detJ * w_ref; // axisymmetric, no 2π

            V_e   += dV;
            r_acc += r * dV;
            z_acc += z * dV;
        }

        if (V_e > 0.0)
        {
            geom.r_centroid[e] = r_acc / V_e;
            geom.z_centroid[e] = z_acc / V_e;
        }
        else
        {
            geom.r_centroid[e] = 0.0;
            geom.z_centroid[e] = 0.0;
        }

        geom.volume[e] = V_e;
    }

    return geom;
}
static const char* ToString(ElectronExitCode c)
{
    switch (c)
    {
        case ElectronExitCode::None:               return "None";
        case ElectronExitCode::HitCathode:         return "HitCathode";
        case ElectronExitCode::HitLiquidGas:       return "HitLiquidGas";
        case ElectronExitCode::HitWall:             return "HitWall";
        case ElectronExitCode::LeftVolume:         return "LeftVolume";
        case ElectronExitCode::MaxSteps:           return "MaxSteps";
        case ElectronExitCode::DegenerateTimeStep: return "DegenerateTimeStep";
        case ElectronExitCode::HitAxis:             return "HitAxis";
        default:                                   return "Unknown";
    }
}
// Do initial column sweep
static void FindWallColumnSeedCI(TracerContext &ctx, 
                                const ElementGeometryData &elem_geom, 
                                bool &seed_found, 
                                int &seed_elem, 
                                double &seed_r, 
                                double &seed_z_boundary,
                                bool debug)
{
    using namespace mfem;

    seed_found      = false;
    seed_elem       = -1;
    seed_r          = 0.0;
    seed_z_boundary = 0.0;

    ParMesh      &global_mesh = ctx.mesh;
    GridFunction &global_phi  = ctx.phi;

    const ElementAdjacency &adj = ctx.adj;
    const Config           &cfg = ctx.cfg;

    // Geometry bounds (used only to generate probe heights)
    const TpcGeometry &tpc = ctx.tpc;
    const double r_max = tpc.r_max;
    const double z_min = tpc.z_min;
    const double z_max = tpc.z_max;

    const double geom_tol    = ctx.params.geom_tol;
    const double min_col_pos = cfg.civ_params.min_col_pos;

    // Wall-near radius (bounds checks handled inside tracer) FIXME : DEBUG
    const double r0 = r_max - 2*geom_tol;

    // Start at 1/4 height
    const double z_start = z_min + 0.5 * (z_max - z_min);

    // Probe schedule: (K-1) halving heights + explicit near-bottom probe
    int K = std::max(2, cfg.compute.threads.num);

    std::vector<double> z_probes;
    z_probes.reserve((std::size_t)K);

    // Halving toward z_min: z = z_min + (z_start - z_min) * 2^{-k}
    // k=0 -> z_start, k=1 -> halfway toward z_min, ...
    for (int k = 0; k < K - 1; ++k)
    {
        const double frac = std::pow(0.5, k);
        double z = z_min + (z_start - z_min) * frac;
        if (z < z_min) z = z_min;
        if (z > z_max) z = z_max;
        z_probes.push_back(z);
    }

    const double z_bottom = std::max(z_min, std::min(z_max, z_min + min_col_pos));
    z_probes.push_back(z_bottom);

    MFEM_VERIFY(!z_probes.empty(), "FindWallColumnSeedCI: z_probes unexpectedly empty");

    // Verify ordering assumption; if violated, warn and recover by sorting desc
    const double last_halving = z_probes[z_probes.size() - 2];

    if (last_halving < z_bottom)
    {
        std::cout << "\033[33m"
                    << "[WARN:CIV] FindWallColumnSeedCI: halving schedule not descending as expected.\n"
                    << "          z_min=" << z_min
                    << " last_halving=" << last_halving
                    << " z_bottom=" << z_bottom
                    << " min_col_pos=" << min_col_pos
                    << " -> sorting probes (desc) to recover.\n"
                    << "          -> min_col_pos likely needs adjusting.\n"
                    << "\033[0m\n";

        std::sort(z_probes.begin(), z_probes.end(), std::greater<double>());
    }
    // ------------------------------------------------------------
    // Map probes to elements/IPs and trace in one batch
    // ------------------------------------------------------------
    CivSeeds seeds;
    seeds.positions.reserve(z_probes.size());
    seeds.elements.reserve(z_probes.size());
    seeds.ips.reserve(z_probes.size());
    seeds.volumes.reserve(z_probes.size()); // unused here

    std::vector<double> z_eval;
    z_eval.reserve(z_probes.size());

    const int dim = global_mesh.SpaceDimension();
    if (dim != 2) { return; }

    // TODO Pragma parallel
    for (double z : z_probes)
    {
        Vector x(dim);
        x[0] = r0;
        x[1] = z;

        int elem = -1;
        IntegrationPoint ip;
        if (!FindElementForPointGlobal(global_mesh, x, elem, ip))
        {
            continue;
        }

        seeds.positions.push_back(x);
        seeds.elements.push_back(elem);
        seeds.ips.push_back(ip);
        seeds.volumes.push_back(1.0);
        z_eval.push_back(z);
    }

    if (seeds.positions.empty())
    {
        throw std::runtime_error("[CRITICAL_ERROR:CIV] FindWallColumnSeedCI: Seed positions empty in column Sweep -> Something is wrong with the mesh (maybe?)");
    }

    std::vector<ElectronTraceResult> results(seeds.positions.size());

    const FiniteElementCollection *fec_phi      = ctx.fes_phi.FEColl();
    const int                      ordering_phi = ctx.fes_phi.GetOrdering();

    const ElectronTraceParams params = ctx.params;

    const bool axisymmetric = cfg.solver.axisymmetric;
    const bool save_paths   = cfg.debug.dumpdata;

    // Per-seed z_max overrides:
    //  - first seed traces to cfg.optimize.z_max
    //  - each subsequent seed traces only up to the previous seed's start height
    std::vector<double> zmax_overrides;
    zmax_overrides.reserve(seeds.positions.size());

    for (std::size_t i = 0; i < z_eval.size(); ++i)
    {
        if (i == 0) { zmax_overrides.push_back(cfg.optimize.z_max); }
        else        { zmax_overrides.push_back(z_eval[i - 1]);     }
    }
    // In case the order is wrong
    std::sort(zmax_overrides.begin(), zmax_overrides.end(), std::greater<double>());

    TraceElectronFieldLinesInner(global_mesh,
                                 global_phi,
                                 adj,
                                 fec_phi,
                                 ordering_phi,
                                 params,
                                 cfg,
                                 seeds,
                                 results,
                                 axisymmetric,
                                 save_paths,
                                 zmax_overrides.data());

    // ------------------------------------------------------------
    // Find first CI (assumes descending z order) and bracket
    // ------------------------------------------------------------
    std::cout << "Entering Loop "<< std::endl;
    int i_first_ci = -1;
    for (int i = 0; i < (int)results.size(); ++i)
    {
        std::cout << "Index " << i << std::boolalpha << IsChargeInsensitive(results[i]) << std::endl;
        if (IsChargeInsensitive(results[i]))
        {
            i_first_ci = i;
            break;
        }
    }

    if (i_first_ci < 0)
    {
        if (debug) { std::cout << "[DEBUG:CIV] FindWallColumnSeedCI: No CI element found in initial column sweep" << std::endl; }
        return; // no CIV on this wall-near column
    }

    // Previous non-CI above the first CI (preferred boundary estimate)
    // TODO Remove the loop when know this works
    int i_prev_non_ci = i_first_ci - 1;


    for (std::size_t i = 0; i < results.size(); ++i) {
        std::cout
            << "Index " << i
            << " r " << seeds.positions[i][0]
            << " z " << seeds.positions[i][1]
            << " exit Condition " << ToString(results[i].exit_code)
            << " charge insensitive? "
            << std::boolalpha
            << IsChargeInsensitive(results[i])
            << std::endl;
    }
    std::cout << "First charge insensitive index " <<i_first_ci <<std::endl;


    if (!IsChargeInsensitive(results[i_first_ci - 1]))
    {
        throw std::runtime_error("[CRITICAL_ERROR:CIV] FindWallColumnSeedCI: Initial CI element search failed -> Implementation Error");
    }

    double z_hi = (i_prev_non_ci >= 0) ? z_eval[i_prev_non_ci] : z_eval[i_first_ci]; // non-CI (preferred)
    double z_lo = z_eval[i_first_ci];                                                // CI

    int elem_lo = seeds.elements[i_first_ci];

    // ------------------------------------------------------------
    // Refine bracket between z_hi (non-CI) and z_lo (CI), if we have both
    // Output boundary = highest non-CI (z_hi).
    // ------------------------------------------------------------
    if (debug) {
    std::cout
        << "[DEBUG:CIV] FindWallColumnSeedCI: Found z bound on first try : "
        << std::boolalpha
        << !(i_prev_non_ci >= 0 && (z_hi - z_lo) > 0.0)
        << std::endl;
    }
    // TODO Change to nthreads samples at a time 
    if (i_prev_non_ci >= 0 && (z_hi - z_lo) > 0.0)
    {
        const double z_tol = std::max(geom_tol, 1e-12);

        CivSeeds one;
        one.positions.resize(1);
        one.elements.resize(1);
        one.ips.resize(1);
        one.volumes.resize(1, 1.0);

        std::vector<ElectronTraceResult> one_res(1);

        // In refinement, keep the same "previous-start" cap: trace up to z_hi (non-CI side)
        double one_zmax = z_hi;

        for (int it = 0; it < 64 && (z_hi - z_lo) > z_tol; ++it)
        {
            const double z_mid = 0.5 * (z_hi + z_lo);

            if (debug)
            {
                std::cout << "[DEBUG:CIV] FindWallColumnSeedCI: testing z_lo : " << z_lo << " z_max : " << z_max << std::endl;
            }

            Vector x(dim);
            x[0] = r0;
            x[1] = z_mid;

            int elem_mid = -1;
            IntegrationPoint ip_mid;
            if (!FindElementForPointGlobal(global_mesh, x, elem_mid, ip_mid))
            {
                z_lo = z_mid;
                continue;
            }

            one.positions[0] = x;
            one.elements[0]  = elem_mid;
            one.ips[0]       = ip_mid;

            TraceElectronFieldLinesInner(global_mesh,
                                         global_phi,
                                         adj,
                                         fec_phi,
                                         ordering_phi,
                                         params,
                                         cfg,
                                         one,
                                         one_res,
                                         axisymmetric,
                                         save_paths,
                                         &one_zmax);

            if (IsChargeInsensitive(one_res[0]))
            {
                z_lo   = z_mid;
                elem_lo = one.elements[0];
            }
            else
            {
                z_hi = z_mid;
                one_zmax = z_hi; // keep cap consistent with current non-CI height
            }
        }
    }

    seed_found      = true;
    seed_elem       = elem_lo; // CI element at/just below boundary
    seed_r          = r0;
    seed_z_boundary = z_hi;    // last known non-CI height ("previous to CI")
}
static void GetFaceNeighbors_LeftEligible(mfem::ParMesh          &mesh,
                                  const ElementAdjacency &adj,
                                  int                     e,
                                  double                  r_ref,
                                  double                  eps_r,
                                  std::vector<int>       &out,
                                  mfem::Array<int>       &verts,
                                  const std::vector<uint8_t> &visited)
{
    out.clear();

    const int ne = mesh.GetNE();
    const auto &face = adj.neighbors[e];

    out.reserve(face.size());

    for (int nb : face)
    {
        if (nb < 0 || nb >= ne) continue;
        if (visited[(std::size_t)nb]) continue;

        verts.SetSize(0);
        mesh.GetElementVertices(nb, verts);

        bool left_ok = false;
        for (int i = 0; i < verts.Size(); ++i)
        {
            const auto *X = mesh.GetVertex(verts[i]);
            if (X[0] <= r_ref + eps_r) { left_ok = true; break; }
        }
        if (left_ok) out.push_back(nb);
    }
}
static void GetVertexNeighbors_LeftEligible(mfem::ParMesh                        &mesh,
                                    const std::vector<std::vector<int>>  &unified_adj,
                                    int                                   e,
                                    double                                r_ref,
                                    double                                eps_r,
                                    std::vector<int>                     &out,
                                    mfem::Array<int>                     &verts,
                                    const std::vector<uint8_t>           &visited)
{
    out.clear();

    const int ne = mesh.GetNE();
    const auto &nbrs = unified_adj[e];

    out.reserve(nbrs.size());

    for (int nb : nbrs)
    {
        if (nb < 0 || nb >= ne) continue;
        if (visited[(std::size_t)nb]) continue;

        verts.SetSize(0);
        mesh.GetElementVertices(nb, verts);

        bool left_ok = false;
        for (int i = 0; i < verts.Size(); ++i)
        {
            const auto *X = mesh.GetVertex(verts[i]);
            if (X[0] <= r_ref + eps_r) { left_ok = true; break; }
        }
        if (left_ok) out.push_back(nb);
    }
}

static void GetNeighbors_DownCandidates(mfem::ParMesh                       &mesh,
                                const std::vector<std::vector<int>> &unified_adj,
                                int                                  e,
                                double                               z_ref,
                                double                               eps_z,
                                std::vector<int>                    &out,
                                mfem::Array<int>                    &verts,
                                const std::vector<uint8_t>          &visited)
{
    out.clear();

    const int ne = mesh.GetNE();
    const auto &nbrs = unified_adj[e];

    out.reserve(nbrs.size());

    for (int nb : nbrs)
    {
        if (nb < 0 || nb >= ne) continue;
        if (visited[(std::size_t)nb]) continue;

        verts.SetSize(0);
        mesh.GetElementVertices(nb, verts);

        bool down_ok = false;
        for (int i = 0; i < verts.Size(); ++i)
        {
            const auto *X = mesh.GetVertex(verts[i]);
            if (X[1] < z_ref - eps_z) { down_ok = true; break; }
        }
        if (down_ok) out.push_back(nb);
    }
}

static bool IsAtBottomBoundary(mfem::ParMesh &mesh,
                        int            e,
                        double         z_min,
                        double         eps_z,
                        mfem::Array<int> &verts)
{
    verts.SetSize(0);
    mesh.GetElementVertices(e, verts);

    for (int i = 0; i < verts.Size(); ++i)
    {
        const auto *X = mesh.GetVertex(verts[i]);
        if (X[1] <= z_min + eps_z) { return true; }
    }
    return false;
}
static std::vector<std::vector<int>> BuildVertexAdjacency(mfem::ParMesh &mesh)
{
    using namespace mfem;

    const int ne = mesh.GetNE();
    const int nv = mesh.GetNV();

    // vertex -> elements containing that vertex
    std::vector<std::vector<int>> vertex_to_elems(nv);

    Array<int> verts;
    for (int e = 0; e < ne; ++e)
    {
        mesh.GetElementVertices(e, verts);
        for (int k = 0; k < verts.Size(); ++k)
        {
            int v = verts[k];
            MFEM_ASSERT(v >= 0 && v < nv, "Invalid vertex index");
            vertex_to_elems[v].push_back(e);
        }
    }

    // element -> vertex neighbors (all elems sharing any vertex)
    std::vector<std::vector<int>> elem_vertex_neighbors(ne);

    for (int e = 0; e < ne; ++e)
    {
        mesh.GetElementVertices(e, verts);
        auto &nbrs = elem_vertex_neighbors[e];

        for (int k = 0; k < verts.Size(); ++k)
        {
            int v = verts[k];
            for (int ee : vertex_to_elems[v])
            {
                if (ee != e)
                {
                    nbrs.push_back(ee);
                }
            }
        }

        // optional: deduplicate
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    return elem_vertex_neighbors;
}
static std::vector<std::vector<int>>
BuildUnifiedAdjacency(const ElementAdjacency &adj,
                      const std::vector<std::vector<int>> &vertex_adj,
                      int ne)
{
    std::vector<std::vector<int>> unified(ne);

    for (int e = 0; e < ne; ++e)
    {
        auto &u = unified[e];

        // insert face neighbors
        const auto &fn = adj.neighbors[e];
        u.insert(u.end(), fn.begin(), fn.end());

        // insert vertex neighbors
        const auto &vn = vertex_adj[e];
        u.insert(u.end(), vn.begin(), vn.end());

        // remove duplicates & self
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), u.end()), u.end());

        // safety: remove the element itself if present
        u.erase(std::remove(u.begin(), u.end(), e), u.end());
    }

    return unified;
}
static void BuildInterfaceSamples_LeftOnly(TracerContext &ctx, 
                                           const ElementGeometryData &elem_geom, 
                                           const std::vector<std::vector<int>> &unified_adj, 
                                           int seed_elem, 
                                           double seed_r, 
                                           double seed_z_boundary, 
                                           std::vector<double> &iface_r, 
                                           std::vector<double> &iface_z)
{
    using namespace mfem;

    iface_r.clear();
    iface_z.clear();

    // Seed sample (wall column boundary from Step 1)
    iface_r.push_back(seed_r);
    iface_z.push_back(seed_z_boundary);

    ParMesh      &mesh = ctx.mesh;
    GridFunction &phi  = ctx.phi;

    const ElementAdjacency &adj = ctx.adj;
    const Config           &cfg = ctx.cfg;

    const ElectronTraceParams params = ctx.params;

    const bool axisymmetric = cfg.solver.axisymmetric;
    const bool save_paths   = cfg.debug.dumpdata;

    const int ne = mesh.GetNE();
    if (seed_elem < 0 || seed_elem >= ne) { return; }

    // Visited to avoid cycles when allowing <= for "left"
    std::vector<uint8_t> visited((std::size_t)ne, 0);
    visited[(std::size_t)seed_elem] = 1;

    const double eps_r = std::max(1e-12, ctx.params.geom_tol);
    const double eps_z = std::max(1e-12, ctx.params.geom_tol);

    Array<int> verts; // scratch

    // Candidate lists (reused to avoid realloc each iteration)
    std::vector<int> face_left;
    std::vector<int> vert_left;
    std::vector<int> down_cands;

    // Mapping helper: local first, fallback global
    auto map_point_to_seed = [&](int elem_hint,
                                 const Vector &x,
                                 int &out_elem,
                                 IntegrationPoint &out_ip) -> bool
    {
        out_elem = -1;
        if (FindElementForPointLocal(mesh, adj, elem_hint, x, out_elem, out_ip)) { return true; }
        return FindElementForPointGlobal(mesh, x, out_elem, out_ip);
    };

    // Trace candidates in-order; return first CI element.
    auto trace_first_ci = [&](const std::vector<int> &cands,
                              int &out_elem_ci) -> bool
    {
        out_elem_ci = -1;
        if (cands.empty()) return false;

        CivSeeds seeds;
        seeds.positions.reserve(cands.size());
        seeds.elements.reserve(cands.size());
        seeds.ips.reserve(cands.size());
        seeds.volumes.reserve(cands.size()); // unused

        std::vector<double> z_start;
        z_start.reserve(cands.size());

        const int dim = mesh.SpaceDimension();

        for (int nb : cands)
        {
            // (visited already filtered by neighbor builders)
            Vector x(dim);
            x[0] = elem_geom.r_centroid[nb];
            x[1] = elem_geom.z_centroid[nb];

            int elem = -1;
            IntegrationPoint ip;
            if (!map_point_to_seed(nb, x, elem, ip)) { continue; }

            seeds.positions.push_back(x);
            seeds.elements.push_back(elem);
            seeds.ips.push_back(ip);
            seeds.volumes.push_back(1.0);
            z_start.push_back(x[1]);
        }

        if (seeds.positions.empty()) return false;

        std::vector<ElectronTraceResult> res(seeds.positions.size());

        const FiniteElementCollection *fec_phi      = ctx.fes_phi.FEColl();
        const int                      ordering_phi = ctx.fes_phi.GetOrdering();

        std::vector<double> zmax_overrides;
        zmax_overrides.reserve(z_start.size());
        for (std::size_t i = 0; i < z_start.size(); ++i)
        {
            if (i == 0) zmax_overrides.push_back(cfg.optimize.z_max);
            else        zmax_overrides.push_back(z_start[i - 1]);
        }

        TraceElectronFieldLinesInner(mesh,
                                     phi,
                                     adj,
                                     fec_phi,
                                     ordering_phi,
                                     params,
                                     cfg,
                                     seeds,
                                     res,
                                     axisymmetric,
                                     save_paths,
                                     zmax_overrides.data());

        for (std::size_t i = 0; i < res.size(); ++i)
        {
            if (IsChargeInsensitive(res[i]))
            {
                out_elem_ci = seeds.elements[i];
                return true;
            }
        }
        return false;
    };

    int e = seed_elem;

    while (true)
    {
        // Current reference vertex: "top-left-most" (maximize z; tie minimize r)
        mesh.GetElementVertices(e, verts);

        double r_ref = std::numeric_limits<double>::infinity();
        double z_ref = -std::numeric_limits<double>::infinity();

        for (int i = 0; i < verts.Size(); ++i)
        {
            const auto *X = mesh.GetVertex(verts[i]);
            const double rv = X[0];
            const double zv = X[1];

            if (zv > z_ref + eps_z)
            {
                z_ref = zv;
                r_ref = rv;
            }
            else if (std::fabs(zv - z_ref) <= eps_z && rv < r_ref)
            {
                r_ref = rv;
            }
        }

        // Build candidates (helpers filter visited and do vertex tests)
        GetFaceNeighbors_LeftEligible(mesh, adj, e, r_ref, eps_r, face_left, verts, visited);
        GetVertexNeighbors_LeftEligible(mesh, unified_adj, e, r_ref, eps_r, vert_left, verts, visited);

        int next_elem = -1;

        // 1) Face-left first CI
        if (!trace_first_ci(face_left, next_elem))
        {
            // 2) Vertex-left first CI
            if (!trace_first_ci(vert_left, next_elem))
            {
                // 3) Down fallback candidates
                GetNeighbors_DownCandidates(mesh, unified_adj, e, z_ref, eps_z, down_cands, verts, visited);

                if (down_cands.empty())
                {
                    if (!IsAtBottomBoundary(mesh, e, ctx.tpc.z_min, eps_z, verts))
                    {
                        std::cout << "\033[33m"
                                  << "[WARN:CIV] Marching stuck before bottom boundary. "
                                  << "Last z_ref=" << z_ref << "\n"
                                  << "\033[0m";
                    }
                    break;
                }

                // Trace all down candidates in one batch, then choose:
                //   smallest z_down_metric (min z among vertices with r <= r_ref+eps_r),
                //   tie: smallest r_down_metric.
                CivSeeds seeds;
                seeds.positions.reserve(down_cands.size());
                seeds.elements.reserve(down_cands.size());
                seeds.ips.reserve(down_cands.size());
                seeds.volumes.reserve(down_cands.size());

                std::vector<double> z_start;
                z_start.reserve(down_cands.size());

                std::vector<double> z_down_metric;
                std::vector<double> r_down_metric;
                z_down_metric.reserve(down_cands.size());
                r_down_metric.reserve(down_cands.size());

                const int dim = mesh.SpaceDimension();

                for (int nb : down_cands)
                {
                    Vector x(dim);
                    x[0] = elem_geom.r_centroid[nb];
                    x[1] = elem_geom.z_centroid[nb];

                    int elem = -1;
                    IntegrationPoint ip;
                    if (!map_point_to_seed(nb, x, elem, ip)) { continue; }

                    seeds.positions.push_back(x);
                    seeds.elements.push_back(elem);
                    seeds.ips.push_back(ip);
                    seeds.volumes.push_back(1.0);
                    z_start.push_back(x[1]);

                    // Metrics for down selection (reuse verts scratch)
                    mesh.GetElementVertices(elem, verts);

                    double z_min_left = std::numeric_limits<double>::infinity();
                    double r_min_left = std::numeric_limits<double>::infinity();

                    for (int vi = 0; vi < verts.Size(); ++vi)
                    {
                        const auto *Xv = mesh.GetVertex(verts[vi]);
                        const double rv = Xv[0];
                        const double zv = Xv[1];

                        if (rv <= r_ref + eps_r)
                        {
                            if (zv < z_min_left) z_min_left = zv;
                            if (rv < r_min_left) r_min_left = rv;
                        }
                    }

                    z_down_metric.push_back(z_min_left);
                    r_down_metric.push_back(r_min_left);
                }

                if (seeds.positions.empty())
                {
                    if (!IsAtBottomBoundary(mesh, e, ctx.tpc.z_min, eps_z, verts))
                    {
                        std::cout << "\033[33m"
                                  << "[WARN:CIV] Down fallback produced no mappable seeds. "
                                  << "Last z_ref=" << z_ref << "\n"
                                  << "\033[0m";
                    }
                    break;
                }

                std::vector<ElectronTraceResult> res(seeds.positions.size());

                const FiniteElementCollection *fec_phi      = ctx.fes_phi.FEColl();
                const int                      ordering_phi = ctx.fes_phi.GetOrdering();

                std::vector<double> zmax_overrides;
                zmax_overrides.reserve(z_start.size());
                for (std::size_t i = 0; i < z_start.size(); ++i)
                {
                    if (i == 0) zmax_overrides.push_back(cfg.optimize.z_max);
                    else        zmax_overrides.push_back(z_start[i - 1]);
                }

                TraceElectronFieldLinesInner(mesh,
                                             phi,
                                             adj,
                                             fec_phi,
                                             ordering_phi,
                                             params,
                                             cfg,
                                             seeds,
                                             res,
                                             axisymmetric,
                                             save_paths,
                                             zmax_overrides.data());

                bool   found_down_ci = false;
                double best_zm       = std::numeric_limits<double>::infinity();
                double best_rm       = std::numeric_limits<double>::infinity();
                int    best_e        = -1;

                for (std::size_t i = 0; i < res.size(); ++i)
                {
                    if (!IsChargeInsensitive(res[i])) continue;

                    const double zm = z_down_metric[i];
                    const double rm = r_down_metric[i];

                    if (!found_down_ci || zm < best_zm - eps_z ||
                        (std::fabs(zm - best_zm) <= eps_z && rm < best_rm))
                    {
                        found_down_ci = true;
                        best_zm       = zm;
                        best_rm       = rm;
                        best_e        = seeds.elements[i];
                    }
                }

                if (!found_down_ci || best_e < 0)
                {
                    if (!IsAtBottomBoundary(mesh, e, ctx.tpc.z_min, eps_z, verts))
                    {
                        std::cout << "\033[33m"
                                  << "[WARN:CIV] Marching terminated: no CI found in down fallback. "
                                  << "Last z_ref=" << z_ref << "\n"
                                  << "\033[0m";
                    }
                    break;
                }

                next_elem = best_e;
            }
        }

        if (next_elem < 0 || next_elem >= ne) break;
        if (visited[(std::size_t)next_elem])  break;

        // Advance
        visited[(std::size_t)next_elem] = 1;
        e = next_elem;

        // Record interface sample for this element: top-left-most vertex
        mesh.GetElementVertices(e, verts);

        double r_s = std::numeric_limits<double>::infinity();
        double z_s = -std::numeric_limits<double>::infinity();

        for (int i = 0; i < verts.Size(); ++i)
        {
            const auto *X = mesh.GetVertex(verts[i]);
            const double rv = X[0];
            const double zv = X[1];

            if (zv > z_s + eps_z)
            {
                z_s = zv;
                r_s = rv;
            }
            else if (std::fabs(zv - z_s) <= eps_z && rv < r_s)
            {
                r_s = rv;
            }
        }

        iface_r.push_back(r_s);
        iface_z.push_back(z_s);
    }
}

static double ComputeCIVFractionFromInterface(const TracerContext &ctx, 
                                              const ElementGeometryData &elem_geom, 
                                              const std::vector<double> &r_sorted, 
                                              const std::vector<double> &z_sorted)
{
    using namespace mfem;

    const ParMesh &mesh = ctx.mesh;
    const int ne = mesh.GetNE();

    if (ne <= 0) { return 0.0; }
    if (r_sorted.empty() || z_sorted.empty()) { return 0.0; }

    const std::size_t n = std::min(r_sorted.size(), z_sorted.size());
    if (n == 0) { return 0.0; }

    // --- Evaluate boundary z(r) by piecewise-linear interpolation (clamped) ---
    auto eval_z_boundary = [&](double r) -> double
    {
        if (n == 1) { return z_sorted[0]; }

        if (r <= r_sorted.front()) { return z_sorted.front(); }
        if (r >= r_sorted.back())  { return z_sorted.back(); }

        // Find segment [i, i+1] with r_sorted[i] <= r <= r_sorted[i+1]
        // (Linear scan is fine for small n; switch to lower_bound if n grows.)
        for (std::size_t i = 0; i + 1 < n; ++i)
        {
            const double r0 = r_sorted[i];
            const double r1 = r_sorted[i + 1];
            if (r >= r0 && r <= r1)
            {
                const double z0 = z_sorted[i];
                const double z1 = z_sorted[i + 1];

                const double dr = (r1 - r0);
                if (dr <= 0.0) { return std::max(z0, z1); }

                const double t = (r - r0) / dr;
                return (1.0 - t) * z0 + t * z1;
            }
        }

        return z_sorted.back();
    };

    // --- Integrate volumes by centroid classification ---
    double V_total = 0.0;
    double V_ci    = 0.0;

    const TpcGeometry &tpc = ctx.tpc;
    const ElementGeometryData &geom = elem_geom;

    // Optional region restriction (if you want full domain, remove these)
    const double z_check_max = ctx.cfg.tracing_params.tracing_z_max;

    for (int e = 0; e < ne; ++e)
    {
        const double V_e = geom.volume[e];
        if (V_e <= 0.0) { continue; }

        const double r_c = geom.r_centroid[e];
        const double z_c = geom.z_centroid[e];

        if (!tpc.Inside(r_c, z_c)) { continue; }
        if (z_c > z_check_max)     { continue; }

        V_total += V_e;

        const double z_b = eval_z_boundary(r_c);
        if (z_c <= z_b)
        {
            V_ci += V_e;
        }
    }

    if (V_total <= 0.0) { return 0.0; }
    return V_ci / V_total;
}


static void SortAndCompactInterfaceByR(const std::vector<double> &iface_r, 
                                       const std::vector<double> &iface_z, 
                                       double r_wall, 
                                       double z_wall, 
                                       double r_radial_min, 
                                       double eps_r, 
                                       double eps_z, 
                                       std::vector<double> &r_sorted, 
                                       std::vector<double> &z_sorted)
{
    r_sorted.clear();
    z_sorted.clear();

    const std::size_t n = std::min(iface_r.size(), iface_z.size());
    if (n == 0) { return; }

    struct P { double r, z; };
    std::vector<P> pts;
    pts.reserve(n + 2);

    for (std::size_t i = 0; i < n; ++i)
    {
        pts.push_back({iface_r[i], iface_z[i]});
    }

    // Ensure the wall point from Step 1 exists (or is at least represented).
    {
        bool have_wall = false;
        for (const auto &p : pts)
        {
            if (std::fabs(p.r - r_wall) <= eps_r && std::fabs(p.z - z_wall) <= eps_z)
            {
                have_wall = true;
                break;
            }
        }
        if (!have_wall)
        {
            pts.push_back({r_wall, z_wall});
        }
    }

    // Sort by radius ascending (required by later interpolation).
    std::sort(pts.begin(), pts.end(),
              [](const P &a, const P &b)
              {
                  if (a.r != b.r) return a.r < b.r;
                  return a.z < b.z;
              });

    // Compact near-duplicate radii; keep the maximum z at that radius.
    r_sorted.reserve(pts.size());
    z_sorted.reserve(pts.size());

    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        const double r = pts[i].r;
        const double z = pts[i].z;

        if (r_sorted.empty())
        {
            r_sorted.push_back(r);
            z_sorted.push_back(z);
            continue;
        }

        const double r_prev = r_sorted.back();
        const double z_prev = z_sorted.back();

        if (std::fabs(r - r_prev) <= eps_r)
        {
            // Warn on duplicates (show both points)
            std::cout << "\033[33m"
                      << "[WARN:CIV] Interface duplicate/near-duplicate radius detected.\n"
                      << "          kept:   (r=" << r_prev << ", z=" << z_prev << ")\n"
                      << "          merged: (r=" << r      << ", z=" << z      << ")\n"
                      << "\033[0m";

            // Keep the larger z as the boundary at that radius.
            if (z > z_prev) { z_sorted.back() = z; }
        }
        else
        {
            r_sorted.push_back(r);
            z_sorted.push_back(z);
        }
    }

    if (r_sorted.empty()) { return; }

    // Ensure we connect horizontally to the radial boundary at the same Z.
    // (If the marched samples do not reach r_radial_min, add (r_radial_min, z_at_min_r)).
    if (r_sorted.front() > r_radial_min + eps_r)
    {
        const double z_min_r = z_sorted.front();
        r_sorted.insert(r_sorted.begin(), r_radial_min);
        z_sorted.insert(z_sorted.begin(), z_min_r);
    }

    // Ensure the curve includes the exact wall radius point (in case sorting/compaction altered it).
    // If the last point isn't at r_wall (within eps_r), append a horizontal closure to r_wall at z_wall.
    if (std::fabs(r_sorted.back() - r_wall) > eps_r)
    {
        r_sorted.push_back(r_wall);
        z_sorted.push_back(z_wall);
    }
}


double ComputeCIV_Marching(const SimulationResult &sim, const Config &cfg)
{
    using namespace mfem;
    // ------------- Set Up
    ParMesh      &mesh = *sim.mesh;
    GridFunction &phi  = *sim.V;
    const FiniteElementSpace      *fes      = phi.FESpace();
    const FiniteElementCollection *fec      = fes->FEColl();
    const int                      ordering = fes->GetOrdering();
    // Precompute objects
    ElementAdjacency               adj         = BuildAdjacency(mesh);
    ElementGeometryData            geom        = ComputeElementGeometry(mesh);
    std::vector<std::vector<int>>  vertex_adj  = BuildVertexAdjacency(mesh);
    std::vector<std::vector<int>>  unified_adj = BuildUnifiedAdjacency(adj, vertex_adj, mesh.GetNE());
    TracerContext ctx(mesh, fec, ordering, phi, adj, cfg);

    // Search along TPC wall halving steps identifying the right most CI to non CI boundary
    // Highest CI element will be placed in seed_elem 
    bool   seed_found      = false;
    int    seed_elem       = -1;
    double seed_r          = 0.0;
    double seed_z_boundary = 0.0;
    FindWallColumnSeedCI(ctx, geom, seed_found, seed_elem, seed_r, seed_z_boundary, cfg.debug.debug);

    if (!seed_found || seed_elem < 0) { return 0.0; }

    // Starting at the last CI element we find the next left CI non CI boundary and repeat this
    // Until there are no more to be found in the general left down direction 
    // This will produce a list of CI element boundaries from which we integrate the CIV 
    std::vector<double> iface_r;
    std::vector<double> iface_z;
    BuildInterfaceSamples_LeftOnly(ctx,
                                   geom,    
                                   unified_adj,
                                   seed_elem,
                                   seed_r,
                                   seed_z_boundary,
                                   iface_r,
                                   iface_z);

    if (iface_r.empty()) { return 0.0; }

    // Step 3: sort + compact interface samples by radius
    TpcGeometry tpc(cfg.tracing_params);
    std::vector<double> r_sorted;
    std::vector<double> z_sorted;
    SortAndCompactInterfaceByR(iface_r, iface_z,
                           /*r_wall=*/seed_r,
                           /*z_wall=*/seed_z_boundary,
                           /*r_radial_min=*/tpc.r_min,
                           /*eps_r=*/cfg.tracing_params.geom_tol,
                           /*eps_z=*/cfg.tracing_params.geom_tol,
                           r_sorted, z_sorted);

    if (r_sorted.empty()) { return 0.0; }

    // Step 4: integrate CIV by classifying elements below z_boundary(r)
    // Returns fraction V_CI / V_total in the region of interest.
    return ComputeCIVFractionFromInterface(ctx, geom, r_sorted, z_sorted);
}




// ============================================================================
// Public API: tracing
// ============================================================================

double compute_civ(const Config &cfg, const SimulationResult &result)
/*
Compute Charge Insensitive Volume on the mesh
Ie 1-(Volume where electrons reach liquid gas interface) / (Total Volume)
*/
{
    auto t_start = std::chrono::steady_clock::now();
    if (cfg.debug.debug)
        std::cout << "[DEBUG:CIV] Timing: start CIV" << std::endl;

    const std::string &method = cfg.civ_params.method; // "InformedSweep" or "RandomSample"

    if (method == "RandomSample")
    {
        return ComputeCIV_RandomSample(cfg, result);
    }
    else if (method == "InformedSweep")
    {
        return ComputeCIV_Marching(result, cfg);
    }
    else if (method == "ColumnSweep")
    {
        if (cfg.debug.debug)
        {
            std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (Column Sweep)" << std::endl;
        }
        double civ = ComputeCIV_ColumnSweep(result, cfg);
        if (cfg.debug.debug)
        {
            auto t_end = std::chrono::steady_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            std::cout << "[DEBUG:CIV] Timing: ColumnSweep took " << ms << " ms" << std::endl;
        }
        return civ;
    }
    else
    {
        // Unknown method: default to random sampling, but warn
        std::cerr << "[WARN:OPTIMIZATION] Unknown CIV method '" << method << "',\n";
        return 1.0;
    }
}

