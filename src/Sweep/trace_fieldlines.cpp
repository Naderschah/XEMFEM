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
#include <random>
#include <deque>
#include <chrono>

#include "trace_fieldlines.h"

using namespace mfem;

// ============================================================================
// CIV by random sample 
// Also contains basic tracing logic
// ============================================================================

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
        V *= sign; // sign = -1.0 gives physical E = -grad(phi)
    }
};


// Build adjacency and characteristic size h per element
static ElementAdjacency BuildAdjacency(mfem::ParMesh &mesh)
{
    ElementAdjacency adj;

    const int ne  = mesh.GetNE();
    const int dim = mesh.SpaceDimension();

    adj.neighbors.resize(ne);
    adj.h.resize(ne);

    mfem::Array<int> vert_ids;

    const mfem::Table &el_to_el = mesh.ElementToElementTable();

    for (int e = 0; e < ne; ++e)
    {
        // neighbors
        {
            mfem::Array<int> nbrs;
            el_to_el.GetRow(e, nbrs);
            adj.neighbors[e].assign(nbrs.begin(), nbrs.end());
        }

        // element "size": shortest edge length
        mesh.GetElementVertices(e, vert_ids);
        const int nv = vert_ids.Size();

        double h_min2 = std::numeric_limits<double>::infinity();

        for (int i = 0; i < nv; ++i)
        {
            const double *xi = mesh.GetVertex(vert_ids[i]);
            for (int j = i+1; j < nv; ++j)
            {
                const double *xj = mesh.GetVertex(vert_ids[j]);
                double dx2 = 0.0;
                for (int d = 0; d < dim; ++d)
                {
                    double diff = xi[d] - xj[d];
                    dx2 += diff * diff;
                }
                if (dx2 < h_min2) { h_min2 = dx2; }
            }
        }

        if (h_min2 == std::numeric_limits<double>::infinity())
        {
            h_min2 = 0.0;
        }

        adj.h[e] = std::sqrt(h_min2);
    }

    return adj;
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


static bool FindElementForPointLocal(
    mfem::ParMesh            &mesh,
    const ElementAdjacency   &adj,
    int                       current_elem,
    const mfem::Vector       &x_new,
    int                      &out_elem,
    mfem::IntegrationPoint   &out_ip)
{
    using namespace mfem;

    // --- 1) Try current element ---
    {
        ElementTransformation *T = mesh.GetElementTransformation(current_elem);
        InverseElementTransformation invT(T);

        IntegrationPoint ip;
        int res = invT.Transform(x_new, ip);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = current_elem;
            out_ip   = ip;
            return true;
        }
    }

    // --- 2) Try neighbors ---
    const auto &nbrs = adj.neighbors[current_elem];
    for (int nb : nbrs)
    {
        if (nb < 0 || nb >= mesh.GetNE()) { continue; }

        ElementTransformation *T = mesh.GetElementTransformation(nb);
        InverseElementTransformation invT(T);

        IntegrationPoint ip;
        int res = invT.Transform(x_new, ip);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = nb;
            out_ip   = ip;
            return true;
        }
    }

    return false;
}

static bool FindElementForPointRobust(
    mfem::ParMesh            &mesh,
    const ElementAdjacency   &adj,
    int                       current_elem,
    const mfem::Vector       &x_new,
    int                      &out_elem,
    mfem::IntegrationPoint   &out_ip)
{
    using namespace mfem;
    // Need this in case we hit a vertex

    // First try cheap local search (current element + neighbors)
    if (FindElementForPointLocal(mesh, adj, current_elem, x_new, out_elem, out_ip))
    {
        return true;
    }

    // Fallback: global search over all elements.
    const int ne = mesh.GetNE();
    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = mesh.GetElementTransformation(e);
        InverseElementTransformation invT(T);

        IntegrationPoint ip;
        int res = invT.Transform(x_new, ip);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = e;
            out_ip   = ip;
            return true;
        }
    }

    // No element contains x_new anywhere in the mesh
    return false;
}

// Single electron tracing (internal helper)
namespace
{
ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh                &mesh,
    ElectricFieldCoeff           &E_coeff,
    const ElementAdjacency       &adj,
    const TpcGeometry            &geom,
    const ElectronTraceParams    &params,
    int                           start_elem,
    const mfem::IntegrationPoint &start_ip,
    bool                          axisymmetric,
    bool                          save_pathlines)
{
    using namespace mfem;

    ElectronTraceResult result;
    result.exit_code    = ElectronExitCode::None;
    result.exit_element = start_elem;

    const int dim = mesh.SpaceDimension();

    int              elem_id = start_elem;
    IntegrationPoint ip      = start_ip;

    Vector x(dim);

    // Initial point
    {
        ElementTransformation *T0 = mesh.GetElementTransformation(elem_id);
        T0->SetIntPoint(&ip);
        T0->Transform(ip, x);

        if (axisymmetric && x[0] <= params.geom_tol)
        {
            x[0] = params.geom_tol;
        }

        if (axisymmetric && params.terminate_on_axis &&
            x[0] <= params.geom_tol)
        {
            result.exit_code    = ElectronExitCode::HitAxis;
            result.exit_element = elem_id;
            if (save_pathlines) { result.points.push_back(x); }
            return result;
        }

        if (save_pathlines) { result.points.push_back(x); }
    }

    int    step       = 0;
    double ds_current = 0.0; // adaptive step length for this electron

    while (step < params.max_steps)
    {
        // Local element size
        const double hK = adj.h[elem_id];
        if (hK <= 0.0)
        {
            result.exit_code    = ElectronExitCode::DegenerateTimeStep;
            result.exit_element = elem_id;
            if (save_pathlines) { result.points.push_back(x); }
            break;
        }

        // Per-element bounds and initial ds
        const double ds_min = params.ds_min_factor * hK;
        const double ds_max = params.ds_max_factor * hK;

        if (ds_current <= 0.0)
        {
            ds_current = params.c_step * hK;
        }
        ds_current = std::min(std::max(ds_current, ds_min), ds_max);

        bool accepted_step = false;

        while (!accepted_step)
        {
            const double ds = ds_current;

            // ------------------ RK4 trial step ------------------
            int              e0  = elem_id;
            IntegrationPoint ip0 = ip;

            // Stage 1
            ElementTransformation *T1 = mesh.GetElementTransformation(e0);
            T1->SetIntPoint(&ip0);

            Vector x0(dim);
            T1->Transform(ip0, x0);

            if (axisymmetric && x0[0] <= params.geom_tol)
            {
                x0[0] = params.geom_tol;
            }

            Vector E1(dim);
            E_coeff.Eval(E1, *T1, ip0);
            double Enorm1 = E1.Norml2();
            if (Enorm1 <= 0.0)
            {
                result.exit_code    = ElectronExitCode::DegenerateTimeStep;
                result.exit_element = e0;
                if (save_pathlines) { result.points.push_back(x0); }
                return result;
            }

            Vector v1(E1);
            v1 *= -1.0 / Enorm1;  // unit direction

            Vector k1(dim), k2(dim), k3(dim), k4(dim);
            k1 = v1;

            bool             need_shrink = false;
            bool             terminated  = false;
            ElectronExitCode term_code   = ElectronExitCode::None;

            Vector          x_new(dim);
            int             elem_new = e0;
            IntegrationPoint ip_new  = ip0;

            // Helper for stages 2–4
            auto rk_stage = [&](double step_factor,
                                const Vector &k_prev,
                                int           search_elem,
                                Vector       &x_stage,
                                int          &elem_stage,
                                IntegrationPoint &ip_stage,
                                Vector       &k_stage)
            {
                if (terminated || need_shrink) { return; }

                x_stage = x0;
                x_stage.Add(step_factor * ds, k_prev);
                if (axisymmetric && x_stage[0] <= params.geom_tol)
                {
                    x_stage[0] = params.geom_tol;
                }

                int              e_tmp  = -1;
                IntegrationPoint ip_tmp;
                if (!FindElementForPointRobust(mesh, adj, search_elem, x_stage, e_tmp, ip_tmp) ||
                        e_tmp < 0 || e_tmp >= mesh.GetNE())
                {
                    need_shrink = true;
                    return;
                }

                ElementTransformation *T = mesh.GetElementTransformation(e_tmp);
                T->SetIntPoint(&ip_tmp);
                Vector E(dim);
                E_coeff.Eval(E, *T, ip_tmp);
                double Enorm = E.Norml2();
                if (Enorm <= 0.0)
                {
                    terminated  = true;
                    term_code   = ElectronExitCode::DegenerateTimeStep;
                    x_new       = x_stage;
                    elem_stage  = e_tmp;
                    ip_stage    = ip_tmp;
                    return;
                }

                Vector v(E);
                v *= -1.0 / Enorm;
                k_stage   = v;
                elem_stage = e_tmp;
                ip_stage   = ip_tmp;
            };

            Vector x2(dim), x3(dim), x4(dim);

            // Stage 2
            rk_stage(0.5, k1, e0, x2, elem_new, ip_new, k2);

            // Stage 3
            rk_stage(0.5, k2, elem_new, x3, elem_new, ip_new, k3);

            // Stage 4
            if (!terminated && !need_shrink)
            {
                rk_stage(1.0, k3, elem_new, x4, elem_new, ip_new, k4);

                if (!terminated && !need_shrink)
                {
                    // Build final x_new from RK4 combination
                    x_new = x0;
                    x_new.Add(ds / 6.0, k1);
                    x_new.Add(ds / 3.0, k2);
                    x_new.Add(ds / 3.0, k3);
                    x_new.Add(ds / 6.0, k4);

                    // Axis termination
                    if (axisymmetric && params.terminate_on_axis &&
                        x_new[0] <= params.geom_tol)
                    {
                        terminated = true;
                        term_code  = ElectronExitCode::HitAxis;
                    }
                    else
                    {
                        // Geometry-based boundary classification
                        const double r_new = x_new[0];
                        const double z_new = x_new[1];

                        const ElectronExitCode bc =
                            geom.ClassifyBoundary(r_new, z_new, params.geom_tol);

                        if (bc != ElectronExitCode::None)
                        {
                            terminated = true;
                            term_code  = bc;
                        }
                    }
                }
            }

            // Handle termination or need_shrink
            if (terminated)
            {
                result.exit_code    = term_code;
                result.exit_element = elem_id; // last known valid in-volume element
                if (save_pathlines) { result.points.push_back(x_new); }
                return result;
            }

            if (need_shrink)
            {
                // Step too large w.r.t. element topology; shrink and retry
                double ds_new = ds_current * params.adapt_shrink;
                if (ds_new < ds_min)
                {
                    // Cannot shrink further: treat as leaving volume
                    result.exit_code    = ElectronExitCode::LeftVolume;
                    result.exit_element = elem_id;
                    if (save_pathlines) { result.points.push_back(x); }
                    return result;
                }
                ds_current = ds_new;
                continue; // retry trial with smaller ds
            }

            // ------------------ Curvature-based "error" ------------------
            Vector dv(dim);
            dv = k4;
            dv -= k1;
            double err_dir = dv.Norml2(); // 0: straight, larger: more bending

            const double tol = params.tol_rel;
            if (err_dir > tol && ds_current > ds_min)
            {
                // Too much bending for current step length → shrink
                double ds_new = ds_current * params.adapt_shrink;
                if (ds_new < ds_min) { ds_new = ds_min; }
                ds_current = ds_new;
                continue; // retry trial with smaller ds
            }

            // Accept step
            elem_id = elem_new;
            ip      = ip_new;
            x       = x_new;

            if (save_pathlines) { result.points.push_back(x_new); }
            result.exit_element = elem_id;
            accepted_step = true;

            // Optional step growth in smooth regions
            if (err_dir < 0.25 * tol)
            {
                double ds_new = ds_current * params.adapt_grow;
                if (ds_new > ds_max) { ds_new = ds_max; }
                ds_current = ds_new;
            }
        } // end adaptive trial loop

        ++step;
    }

    if (result.exit_code == ElectronExitCode::None)
    {
        result.exit_code = ElectronExitCode::MaxSteps;
    }

    return result;
}
} // anonymous namespace


// ============================================================================
// Public API: seed extraction
// ============================================================================

CivSeeds ExtractCivSeeds(const Config &cfg, const SimulationResult &result)
{
    mfem::ParMesh &pmesh = *result.mesh;
    const int dim  = pmesh.Dimension();

    TpcGeometry geom(cfg);

    const int ir_order = cfg.tracing_params.ir_order;

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
    int num_seed_elements = cfg.tracing_params.num_seed_elements;
    const int n_eligible  = static_cast<int>(eligible_elems.size());

    if (num_seed_elements <= 0 || num_seed_elements > n_eligible)
    {
        num_seed_elements = n_eligible; // fall back to all eligible
    }

    std::mt19937 rng;
    if (cfg.tracing_params.rng_seed != 0)
    {
        rng.seed(static_cast<std::mt19937::result_type>(cfg.tracing_params.rng_seed));
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

// Old random-sampling CIV computation, factored out so we can dispatch by method.
static double ComputeCIV_RandomSample(const Config           &cfg,
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
        if (res.exit_code != ElectronExitCode::HitLiquidGas)
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

void TraceElectronFieldLines(const SimulationResult           &sim,
                             const Config                     &cfg,
                             const CivSeeds                   &seeds,
                             std::vector<ElectronTraceResult> &out_results)
{
    using namespace mfem;

    // Copy tracing params so we can tweak them locally (e.g. in debug mode)
    ElectronTraceParams params = cfg.tracing_params;

    ParMesh &global_mesh = *sim.mesh;
    const int dim = global_mesh.SpaceDimension();

    // H1 potential used to compute E = -grad(phi)
    GridFunction &global_phi = *sim.V;

    const std::size_t n_seeds = seeds.positions.size();
    MFEM_VERIFY(seeds.elements.size() == n_seeds &&
                seeds.ips.size()      == n_seeds,
                "TraceElectronFieldLines: CivSeeds size mismatch");

    // Prepare output container for parallel write by index
    out_results.clear();
    out_results.resize(n_seeds);

    // Build adjacency ONCE from the global mesh
    ElementAdjacency adj = BuildAdjacency(global_mesh);

    TpcGeometry geom(cfg);

    // FE info for the potential phi
    const FiniteElementSpace        *fes_phi   = global_phi.FESpace();
    const FiniteElementCollection   *fec_phi   = fes_phi->FEColl();
    const int                        ordering_phi = fes_phi->GetOrdering();

    const bool axisymmetric = cfg.solver.axisymmetric;
    const bool save_paths   = cfg.debug.dumpdata;

    // -------------------------------------------------------------------------
    // Single-seed debug mode
    // -------------------------------------------------------------------------
    if (cfg.debug.debug_single_seed)
    {
        // 1) Construct centre point in physical coordinates
        const double r_c = 0.5 * (params.r_min + params.r_max);
        const double z_c = 0.5 * (params.z_min + params.z_max);

        mfem::Vector x_center(dim);
        x_center[0] = r_c;
        x_center[1] = z_c;

        // 2) Find element and reference IntegrationPoint that contain x_center
        int                   elem_center = -1;
        mfem::IntegrationPoint ip_center;

        for (int e = 0; e < global_mesh.GetNE(); ++e)
        {
            mfem::ElementTransformation *T = global_mesh.GetElementTransformation(e);
            mfem::InverseElementTransformation invT(T);

            mfem::IntegrationPoint ip_trial;
            int res = invT.Transform(x_center, ip_trial);
            if (res == InverseElementTransformation::Inside)
            {
                elem_center = e;
                ip_center   = ip_trial;
                break;
            }
        }

        mfem::ElementTransformation *Tcheck =
            global_mesh.GetElementTransformation(elem_center);
        Tcheck->SetIntPoint(&ip_center);

        mfem::Vector x_check(dim);
        Tcheck->Transform(ip_center, x_check);

        std::cout << "[DEBUG:TRACE] Centre from config : r=" << r_c
                  << ", z=" << z_c << "\n";
        std::cout << "[DEBUG:TRACE] Mapped back centre : r=" << x_check[0]
                  << ", z=" << x_check[1] << "\n";

        if (cfg.debug.debug)
        {
            std::cout << "[DEBUG:TRACE] Single-seed debug mode active.\n"
                      << "               Centre point: r=" << r_c
                      << ", z=" << z_c << "\n"
                      << "               Element id : " << elem_center << "\n";
        }

        // 3) Local mesh and potential
        mfem::ParMesh local_mesh(global_mesh);

        mfem::ParFiniteElementSpace local_fes_phi(
            &local_mesh, fec_phi, /*vdim=*/1, ordering_phi);

        mfem::GridFunction local_phi(&local_fes_phi);
        local_phi = global_phi; // copy potential DOFs

        // Local E coefficient: E = -grad(phi)
        ElectricFieldCoeff local_E_coeff(local_phi, -1.0);

        const ElementAdjacency &local_adj = adj;

        // 4) Local params with optional c_step override
        ElectronTraceParams local_params = params;
        if (cfg.debug.debug_c_step_override > 0.0)
        {
            local_params.c_step = cfg.debug.debug_c_step_override;
        }

        // 5) Resize outputs and trace exactly one electron
        out_results.clear();
        out_results.resize(1);

        out_results[0] = TraceSingleElectronLine(local_mesh,
                                                 local_E_coeff,
                                                 local_adj,
                                                 geom,
                                                 local_params,
                                                 elem_center,
                                                 ip_center,
                                                 axisymmetric,
                                                 save_paths);
    }
    else
    {
        // ---------------------------------------------------------------------
        // Normal mode: trace all seeds in parallel
        // ---------------------------------------------------------------------

        #pragma omp parallel
        {
            // Local copy of the mesh (geometry + topology)
            ParMesh local_mesh(global_mesh);

            // Local FE space and potential on local mesh
            ParFiniteElementSpace local_fes_phi(&local_mesh, fec_phi, /*vdim=*/1, ordering_phi);
            GridFunction local_phi(&local_fes_phi);
            local_phi = global_phi; // copy potential DOFs

            // Local E coefficient: physical field E = -grad(phi)
            ElectricFieldCoeff local_E_coeff(local_phi, -1.0);

            // Adjacency can be shared
            const ElementAdjacency &local_adj = adj;

            // Parallelize over seeds within this thread's context
            #pragma omp for schedule(static)
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(n_seeds); ++i)
            {
                const int start_elem = seeds.elements[i];
                const IntegrationPoint &start_ip = seeds.ips[i];

                ElectronTraceParams local_params = params;

                ElectronTraceResult res = TraceSingleElectronLine(local_mesh,
                                                                  local_E_coeff,
                                                                  local_adj,
                                                                  geom,
                                                                  local_params,
                                                                  start_elem,
                                                                  start_ip,
                                                                  axisymmetric,
                                                                  save_paths);

                out_results[static_cast<std::size_t>(i)] = std::move(res);
            }
        } // end parallel
    }

    // Print exit condition summary (serial)
    if (cfg.debug.debug)
    {
        std::size_t count_hit_lgi      = 0;
        std::size_t count_hit_cathode  = 0;
        std::size_t count_hit_wall     = 0;
        std::size_t count_left_volume  = 0;
        std::size_t count_max_steps    = 0;
        std::size_t count_deg_dt       = 0;
        std::size_t count_none         = 0;
        std::size_t count_hit_axis     = 0; 

        for (const auto &res : out_results)
        {
            switch (res.exit_code)
            {
                case ElectronExitCode::HitLiquidGas:       count_hit_lgi++;     break;
                case ElectronExitCode::HitCathode:         count_hit_cathode++; break;
                case ElectronExitCode::HitWall:            count_hit_wall++;    break;
                case ElectronExitCode::LeftVolume:         count_left_volume++; break;
                case ElectronExitCode::MaxSteps:           count_max_steps++;   break;
                case ElectronExitCode::DegenerateTimeStep: count_deg_dt++;      break;
                case ElectronExitCode::HitAxis:            count_hit_axis++;    break;
                case ElectronExitCode::None:               count_none++;        break;
            }
        }

        const double N = static_cast<double>(out_results.size());
        auto frac = [&](std::size_t c) -> double {
            return (N > 0.0 ? static_cast<double>(c) / N : 0.0);
        };

        double V_total        = 0.0;
        double V_hit_lgi      = 0.0;
        double V_hit_cathode  = 0.0;
        double V_hit_wall     = 0.0;
        double V_left_volume  = 0.0;
        double V_max_steps    = 0.0;
        double V_deg_dt       = 0.0;
        double V_none         = 0.0;
        double V_hit_axis     = 0.0;

        for (std::size_t i = 0; i < seeds.positions.size(); ++i)
        {
            const double dV = seeds.volumes[i];
            const auto  &res = out_results[i];

            V_total += dV;

            switch (res.exit_code)
            {
                case ElectronExitCode::HitLiquidGas:       V_hit_lgi     += dV; break;
                case ElectronExitCode::HitCathode:         V_hit_cathode += dV; break;
                case ElectronExitCode::HitWall:            V_hit_wall    += dV; break;
                case ElectronExitCode::LeftVolume:         V_left_volume += dV; break;
                case ElectronExitCode::MaxSteps:           V_max_steps   += dV; break;
                case ElectronExitCode::DegenerateTimeStep: V_deg_dt      += dV; break;
                case ElectronExitCode::HitAxis:            V_hit_axis    += dV; break;
                case ElectronExitCode::None:               V_none        += dV; break;
            }
        }

        auto vfrac = [&](double V) -> double {
            return (V_total > 0.0 ? V / V_total : 0.0);
        };

        std::cout << "\n---------------- Electron Tracing Summary ----------------\n";
        std::cout << "Total seeds traced: " << out_results.size() << "\n\n";

        std::cout << std::left
                  << std::setw(22) << "Exit Type"
                  << std::setw(18) << "Seed Fraction"
                  << std::setw(18) << "Volume Fraction\n";

        std::cout << "-----------------------------------------------------------------\n";

        auto row = [&](const char *label, double f_seed, double f_vol) {
            std::cout << std::left
                      << std::setw(22) << label
                      << std::setw(18) << f_seed
                      << std::setw(18) << f_vol
                      << "\n";
        };

        row("Hit Liquid-Gas",
            frac(count_hit_lgi),
            vfrac(V_hit_lgi));

        row("Hit Cathode",
            frac(count_hit_cathode),
            vfrac(V_hit_cathode));

        row("Hit Wall",
            frac(count_hit_wall),
            vfrac(V_hit_wall));

        row("Hit Axis",
            frac(count_hit_axis),
            vfrac(V_hit_axis));

        row("Left Volume",
            frac(count_left_volume),
            vfrac(V_left_volume));

        row("Max Steps",
            frac(count_max_steps),
            vfrac(V_max_steps));

        row("Degenerate dt",
            frac(count_deg_dt),
            vfrac(V_deg_dt));

        row("None (unexpected)",
            frac(count_none),
            vfrac(V_none));

        std::cout << "-----------------------------------------------------------------\n";
    }

    if (cfg.debug.dumpdata)
    {
        std::cout << "[DEBUG] Dumping Electron Paths in CIV" << std::endl;
        namespace fs = std::filesystem;
        fs::path outdir(cfg.save_path);
        fs::path outfile = outdir / "electron_paths_debug.csv";
        std::ofstream ofs(outfile);
        ofs << "# id, step, r, z, exit_code\n";

        for (std::size_t i = 0; i < out_results.size(); ++i)
        {
            const auto &res = out_results[i];
            const auto &pts = res.points;

            for (std::size_t k = 0; k < pts.size(); ++k)
            {
                const mfem::Vector &x = pts[k];
                const double r = x(0);
                const double z = x(1);

                ofs << i << "," << k << "," << r << "," << z << ","
                    << static_cast<int>(res.exit_code) << "\n";
            }
        }
    }
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
// Return true if this exit is cathode-based charge insensitivity
static bool IsCathodeInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitCathode);
}

// Return true if this exit is wall-based charge insensitivity
static bool IsWallInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitWall);
}

// Return true for any form of charge insensitivity (cathode or wall)
static bool IsChargeInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitCathode ||
            res.exit_code == ElectronExitCode::HitWall);
}

// Axisymmetric element centroids and volumes (no 2π factor)
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
struct TracerContext
{
    // Per-thread-owned MFEM objects
    mfem::ParMesh               mesh;        // local copy of global mesh
    mfem::ParFiniteElementSpace fes_phi;     // H1 space for potential
    mfem::GridFunction          phi;         // potential on local mesh
    ElectricFieldCoeff          E_coeff;     // E = -grad(phi)

    // Shared read-only data (references or pointers)
    const ElementAdjacency              &adj;
    const ElementGeometryData           &geom;
    const std::vector<std::vector<int>> &unified_adj;

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
                  const ElementGeometryData          &geom_ref,
                  const std::vector<std::vector<int>> &unified_adj_ref,
                  const Config                       &cfg)
        : mesh(global_mesh),
          fes_phi(&mesh, fec_phi, /*vdim=*/1, ordering),
          phi(&fes_phi),
          E_coeff(phi, -1.0),
          adj(adj_ref),
          geom(geom_ref),
          unified_adj(unified_adj_ref),
          tpc(cfg),
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
// Parallel interface
static void ProbeBatch(const std::vector<ProbeRequest>       &requests,
                       const std::vector<TracerContext*>     &tctx_pool,
                       std::vector<ProbeResult>              &out_results)
{
    const int N = (int)requests.size();
    out_results.clear();
    out_results.resize(N);

    if (N == 0 || tctx_pool.empty())
    {
        return;
    }

    const int n_ctx = (int)tctx_pool.size();
    const int max_threads = std::max(1, n_ctx);

    #pragma omp parallel for num_threads(max_threads)
    for (int i = 0; i < N; ++i)
    {
        int tid = omp_get_thread_num();
        // Safety in case OpenMP chooses more threads than tctx_pool size
        if (tid >= n_ctx) tid = tid % n_ctx;

        TracerContext &ctx = *tctx_pool[tid];

        const ProbeRequest &req = requests[i];
        ProbeResult        &out = out_results[i];

        int elem_guess = req.elem_guess;
        ElectronTraceResult res;

        bool ok = TraceSingleCIProbe(ctx,
                                     req.r,
                                     req.z,
                                     elem_guess,
                                     res);

        out.ok         = ok;
        out.r          = req.r;
        out.z          = req.z;
        out.elem       = elem_guess;
        out.exit_code  = ok ? int(res.exit_code) : -1;
        out.is_ci      = ok && IsChargeInsensitive(res);
        out.is_cathode = ok && IsCathodeInsensitive(res);
        out.is_wall    = ok && IsWallInsensitive(res);
    }
}

// Build CivSeeds from a list of physical points (r,z)
// Any point that cannot be mapped to an element is skipped.
static CivSeeds BuildSeedsFromPhysicalPoints(mfem::ParMesh                         &mesh,
                                             const std::vector<mfem::Vector>       &points)
{
    using namespace mfem;

    CivSeeds seeds;

    const int dim = mesh.SpaceDimension();
    MFEM_VERIFY(dim == 2, "BuildSeedsFromPhysicalPoints: only 2D axisymmetric meshes supported.");

    const std::size_t n_pts = points.size();
    seeds.positions.reserve(n_pts);
    seeds.elements.reserve(n_pts);
    seeds.ips.reserve(n_pts);
    seeds.volumes.reserve(n_pts); // will typically be dummy in interface search

    for (const auto &x_phys : points)
    {
        MFEM_VERIFY(x_phys.Size() == dim,
                    "BuildSeedsFromPhysicalPoints: point dimension mismatch.");

        int elem = -1;
        IntegrationPoint ip;

        if (!FindElementForPointGlobal(mesh, x_phys, elem, ip))
        {
            // Point lies outside the mesh; skip
            continue;
        }

        seeds.positions.push_back(x_phys);
        seeds.elements.push_back(elem);
        seeds.ips.push_back(ip);

        // In the interface search we don't use per-seed volumes for CIV;
        // set a dummy value (e.g. 1.0) or 0.0 as a placeholder.
        seeds.volumes.push_back(1.0);
    }

    return seeds;
}
// Trace a batch of electrons starting from physical (r,z) positions.
// Points that cannot be mapped to an element are silently skipped;
// out_results will have the same length as the resulting seeds.
static void TraceBatchFromPhysicalPoints(const SimulationResult         &sim,
                                         const Config                   &cfg,
                                         const std::vector<mfem::Vector> &points,
                                         std::vector<ElectronTraceResult> &out_results)
{
    using namespace mfem;

    ParMesh &mesh = *sim.mesh;

    // Build seeds (this will skip points outside the mesh)
    CivSeeds seeds = BuildSeedsFromPhysicalPoints(mesh, points);

    out_results.clear();

    const std::size_t n_seeds = seeds.positions.size();
    if (n_seeds == 0)
    {
        return; // nothing to trace
    }

    out_results.resize(n_seeds);

    if (cfg.debug.debug) {std::cout << "Tracing Batch of electorns n_elems=" << n_seeds << std::endl;}

    TraceElectronFieldLines(sim, cfg, seeds, out_results);
}
// Search along a vertical column at fixed r for charge-insensitive exits.
// We generate a geometric sequence of z-values between z_mid and z_min
// (towards the cathode) and trace all of them in one batch.
static ColumnCISearchResult SearchColumnForCI(const SimulationResult &sim,
                                              const Config           &cfg,
                                              double                  r_fixed,
                                              double                  z_min,
                                              double                  z_max)
{
    using namespace mfem;

    ColumnCISearchResult col_res;
    col_res.found         = false;
    col_res.z_ci          = z_min;
    col_res.is_cathode_ci = false;
    col_res.is_wall_ci    = false;

    const int dim = sim.mesh->SpaceDimension();
    MFEM_VERIFY(dim == 2, "SearchColumnForCI: only 2D axisymmetric meshes supported.");

    const double z_mid = 0.5 * (z_min + z_max);

    // Build a sequence of test heights: z_mid, then halfway towards z_min,
    // then halfway again, etc. We keep them all and trace in one batch.
    const int max_levels = 8; // enough to get close to cathode

    std::vector<double> z_tests;
    z_tests.reserve(max_levels);

    double z_hi = z_mid;
    double z_lo = z_min;

    // z0 = mid, then mid - (mid - z_min)/2, etc
    for (int k = 0; k < max_levels; ++k)
    {
        double zk;
        if (k == 0)
        {
            zk = z_mid;
        }
        else
        {
            double frac = std::pow(0.5, k); // 1/2, 1/4, 1/8, ...
            zk = z_mid - (z_mid - z_min) * frac;
        }
        if (zk < z_min) { zk = z_min; }
        if (zk > z_max) { zk = z_max; }

        // Avoid duplicates if rounding collapses
        if (!z_tests.empty() && std::fabs(zk - z_tests.back()) < 1e-12)
        {
            continue;
        }

        z_tests.push_back(zk);
    }

    if (z_tests.empty())
    {
        return col_res;
    }

    // Build physical points along this column
    std::vector<Vector> points;
    points.reserve(z_tests.size());

    for (double z : z_tests)
    {
        Vector x_phys(dim);
        x_phys[0] = r_fixed;
        x_phys[1] = z;
        points.push_back(x_phys);
    }

    // Trace all in one batch
    std::vector<ElectronTraceResult> results;
    TraceBatchFromPhysicalPoints(sim, cfg, points, results);

    // Note: some points may have been skipped by BuildSeedsFromPhysicalPoints
    // if they were outside the mesh. So we need to be careful: we assume
    // that if a point could not be mapped, it simply won't appear in results.
    // For a well-chosen r_fixed,z-range this should be rare; if it happens,
    // we effectively ignore that z in the search.

    // Walk through results to find highest z that is CI
    double z_ci_max         = z_min;
    bool   found_any_ci     = false;
    bool   any_cathode_ci   = false;
    bool   any_wall_ci      = false;

    // We assume results are in the same order as 'points'
    const std::size_t n_res = results.size();
    const std::size_t n_pts = points.size();
    const std::size_t n_eval = std::min(n_res, n_pts);

    for (std::size_t i = 0; i < n_eval; ++i)
    {
        const ElectronTraceResult &res = results[i];
        const Vector              &x0  = points[i];
        const double z              = x0[1];

        if (IsChargeInsensitive(res))
        {
            found_any_ci = true;
            if (z > z_ci_max) { z_ci_max = z; }

            if (IsCathodeInsensitive(res)) { any_cathode_ci = true; }
            if (IsWallInsensitive(res))    { any_wall_ci    = true; }
        }
    }

    if (!found_any_ci)
    {
        return col_res; // no CI in this column
    }

    col_res.found         = true;
    col_res.z_ci          = z_ci_max;
    col_res.is_cathode_ci = any_cathode_ci;
    col_res.is_wall_ci    = any_wall_ci;

    return col_res;
}
static ColumnRefineResult RefineColumnCI(const SimulationResult &sim,
                                         const Config           &cfg,
                                         double                  r_fixed,
                                         double                  z_lo,
                                         double                  z_hi,
                                         int                     n_samples)
{
    using namespace mfem;

    ColumnRefineResult ref;
    ref.found         = false;
    ref.z_ci          = z_lo;
    ref.is_cathode_ci = false;
    ref.is_wall_ci    = false;

    const int dim = sim.mesh->SpaceDimension();
    MFEM_VERIFY(dim == 2, "RefineColumnCI: only 2D axisymmetric meshes supported.");

    if (n_samples <= 0 || z_hi <= z_lo)
    {
        return ref;
    }

    std::vector<Vector> points;
    points.reserve(n_samples);

    const double dz = (z_hi - z_lo) / static_cast<double>(n_samples + 1);

    for (int k = 1; k <= n_samples; ++k)
    {
        double z = z_lo + k * dz;
        Vector x_phys(dim);
        x_phys[0] = r_fixed;
        x_phys[1] = z;
        points.push_back(x_phys);
    }

    std::vector<ElectronTraceResult> results;
    TraceBatchFromPhysicalPoints(sim, cfg, points, results);

    const std::size_t n_res = results.size();
    const std::size_t n_pts = points.size();
    const std::size_t n_eval = std::min(n_res, n_pts);

    double z_ci_max       = z_lo;
    bool   found_any_ci   = false;
    bool   any_cathode_ci = false;
    bool   any_wall_ci    = false;

    for (std::size_t i = 0; i < n_eval; ++i)
    {
        const ElectronTraceResult &res = results[i];
        const Vector              &x0  = points[i];
        const double z              = x0[1];

        if (IsChargeInsensitive(res))
        {
            found_any_ci = true;
            if (z > z_ci_max) { z_ci_max = z; }

            if (IsCathodeInsensitive(res)) { any_cathode_ci = true; }
            if (IsWallInsensitive(res))    { any_wall_ci    = true; }
        }
    }

    if (!found_any_ci)
    {
        return ref;
    }

    ref.found         = true;
    ref.z_ci          = z_ci_max;
    ref.is_cathode_ci = any_cathode_ci;
    ref.is_wall_ci    = any_wall_ci;

    return ref;
}


static InitialColumnCIResult FindInitialColumnCI(TracerContext &base_ctx)
{
    using namespace mfem;

    InitialColumnCIResult out;
    out.found      = false;
    out.r          = 0.0;
    out.z_ci       = 0.0;
    out.elem       = -1;
    out.cathode_ci = false;
    out.wall_ci    = false;

    const TpcGeometry        &tpc    = base_ctx.tpc;
    const ElectronTraceParams &params = base_ctx.params;
    const Config &cfg = base_ctx.cfg;
    int n_threads = cfg.compute.threads.num;
    bool debug = cfg.debug.debug;

    const double r_min    = tpc.r_min;
    const double r_max    = tpc.r_max;
    const double z_min    = tpc.z_min;
    const double z_max    = tpc.z_max;
    const double geom_tol = params.geom_tol;

    const double r0 = r_max - geom_tol;
    if (r0 <= r_min + geom_tol) {
        return out;
    }

    const double z_mid = 0.5 * (z_min + z_max);

    // 1) Build a batch of candidate z’s between z_mid and z_min
    const int max_samples = std::max(4, n_threads);  // at least a few, at most a few dozen
    std::vector<ProbeRequest> reqs;
    reqs.reserve(max_samples);

    for (int k = 0; k < max_samples; ++k) {
        // Geometric-ish spacing from mid downwards
        double frac = std::pow(0.5, k); // 1, 1/2, 1/4, ...
        double z    = z_mid - (z_mid - z_min) * frac;
        if (z < z_min) z = z_min;

        ProbeRequest pr;
        pr.r         = r0;
        pr.z         = z;
        pr.elem_guess = -1;
        reqs.push_back(pr);
    }

    // 2) Build per-thread contexts as in 1.3
    int n_used_threads = std::max(1, n_threads);
    std::vector<std::unique_ptr<TracerContext>> extra;
    extra.reserve(n_used_threads - 1);

    std::vector<TracerContext*> ctx_pool(n_used_threads);
    ctx_pool[0] = &base_ctx;

    const FiniteElementCollection *fec = base_ctx.fes_phi.FEColl();
    const int                      ord = base_ctx.fes_phi.GetOrdering();

    for (int t = 1; t < n_used_threads; ++t) {
        extra.emplace_back(std::make_unique<TracerContext>(
            base_ctx.mesh, fec, ord, base_ctx.phi,
            base_ctx.adj, base_ctx.geom, base_ctx.unified_adj, cfg));
        ctx_pool[t] = extra.back().get();
    }

    // 3) Run a single parallel ProbeBatch
    std::vector<ProbeResult> results;
    ProbeBatch(reqs, ctx_pool, results);
    

    // 4) Select the highest CI z
    bool   found_ci    = false;
    double best_z      = z_min;
    int    best_elem   = -1;
    bool   saw_cathode = false;
    bool   saw_wall    = false;

    for (std::size_t i = 0; i < results.size(); ++i) {
        const ProbeResult &pr = results[i];
        if (!pr.ok || !pr.is_ci) continue;

        if (!found_ci || pr.z > best_z) {
            found_ci  = true;
            best_z    = pr.z;
            best_elem = pr.elem;
        }

        saw_cathode = saw_cathode || pr.is_cathode;
        saw_wall    = saw_wall    || pr.is_wall;
    }

    if (!found_ci) {
        return out;
    }

    // Optional: small serial refinement around best_z (reuse old bisection if you like),
    // but now centered near best_z and with reduced z-interval.

    out.found      = true;
    out.r          = r0;
    out.z_ci       = best_z;
    out.elem       = best_elem;
    out.cathode_ci = saw_cathode;
    out.wall_ci    = saw_wall;

    if (debug) {
        std::cout << "[DEBUG:MARCH] InitialColumnCI_Parallel: r=" << out.r
                  << " z_ci=" << out.z_ci
                  << " elem=" << out.elem
                  << " cath=" << out.cathode_ci
                  << " wall=" << out.wall_ci << "\n";
    }

    return out;
}


// Evaluate CI boundary height z_b at radius r using linear interpolation
// over the given interface points (assumed sorted by r).
// If r is outside the sampled range, we extrapolate using the nearest point.
static double EvaluateInterfaceZ(const std::vector<InterfacePoint> &iface,
                                 double                             r)
{
    const std::size_t n = iface.size();
    if (n == 0)
    {
        return 0.0;
    }

    if (n == 1)
    {
        return iface[0].z_ci;
    }

    // If r is below the smallest sampled radius, clamp to first point
    if (r <= iface.front().r)
    {
        return iface.front().z_ci;
    }

    // If r is above the largest sampled radius, clamp to last point
    if (r >= iface.back().r)
    {
        return iface.back().z_ci;
    }

    // Find segment [i, i+1] s.t. iface[i].r <= r <= iface[i+1].r
    for (std::size_t i = 0; i + 1 < n; ++i)
    {
        const double r0 = iface[i].r;
        const double r1 = iface[i+1].r;
        if (r >= r0 && r <= r1)
        {
            const double z0 = iface[i].z_ci;
            const double z1 = iface[i+1].z_ci;

            const double t = (r - r0) / (r1 - r0);
            return (1.0 - t) * z0 + t * z1;
        }
    }

    // Fallback (should not reach here if above logic is correct)
    return iface.back().z_ci;
}
// Classify elements as charge-insensitive by comparing their centroids
// to the interface curve, then sum volumes.
static CIVolumeResult ComputeCIVFromInterface(mfem::ParMesh                 &mesh,
                                              const Config                  &cfg,
                                              const ElementGeometryData     &geom,
                                              const std::vector<InterfacePoint> &iface)
{
    using namespace mfem;

    CIVolumeResult res;
    res.V_total = 0.0;
    res.V_CI    = 0.0;

    const int ne = mesh.GetNE();
    if (ne <= 0 || iface.empty())
    {
        return res;
    }

    TpcGeometry tpc(cfg);
    const double z_check_max = cfg.tracing_params.z_max; // or a dedicated z_check_max if you add it

    for (int e = 0; e < ne; ++e)
    {
        const double r_c = geom.r_centroid[e];
        const double z_c = geom.z_centroid[e];
        const double V_e = geom.volume[e];

        if (V_e <= 0.0)
        {
            continue;
        }

        // Optionally restrict CIV computation to lower part of TPC
        if (z_c > z_check_max)
        {
            continue;
        }

        // Ignore elements outside the TPC cylinder
        if (!tpc.Inside(r_c, z_c))
        {
            continue;
        }

        res.V_total += V_e;

        // Interface-based CI decision: below or on the curve => CI
        const double z_boundary = EvaluateInterfaceZ(iface, r_c);
        if (z_c <= z_boundary)
        {
            res.V_CI += V_e;
        }
    }

    return res;
}

static std::vector<InterfacePoint>
BuildInterfaceCurve_Marching(TracerContext &ctx)
{
    using namespace mfem;

    ParMesh                   &mesh    = ctx.mesh;
    const ElementAdjacency    &adj     = ctx.adj;
    const ElementGeometryData &geom    = ctx.geom;
    const TpcGeometry         &tpc     = ctx.tpc;
    const ElectronTraceParams &params  = ctx.params;
    const auto                &unified_adj = ctx.unified_adj; // <--- use unified adjacency
    bool debug                         = ctx.cfg.debug.debug;

    if (debug)
    {
        std::cout << "[DEBUG:MARCH] BuildInterfaceCurve_Marching\n";
        std::cout << "[DEBUG:MARCH] Mesh NE=" << mesh.GetNE() << "\n";
        std::cout << "[DEBUG:MARCH] r_min=" << tpc.r_min
                  << " r_max=" << tpc.r_max
                  << " z_min=" << tpc.z_min
                  << " z_max=" << tpc.z_max
                  << " geom_tol=" << params.geom_tol << "\n";
    }

    std::vector<InterfacePoint> iface;

    const double r_min    = tpc.r_min;
    const double r_max    = tpc.r_max;
    const double z_min    = tpc.z_min;
    const double z_max    = tpc.z_max;
    const double geom_tol = params.geom_tol;

    // ---------------------------------------------------------------------
    // 1) Initial column CI (at r ~ r_max)
    // ---------------------------------------------------------------------
    InitialColumnCIResult init = FindInitialColumnCI(ctx);

    if (debug)
    {
        std::cout << "[DEBUG:MARCH] InitialColumnCI: found=" << init.found
                  << " r=" << init.r
                  << " z_ci=" << init.z_ci
                  << " elem=" << init.elem
                  << " cathode_ci=" << init.cathode_ci
                  << " wall_ci=" << init.wall_ci << "\n";
    }

    if (!init.found || init.elem < 0)
    {
        if (debug)
            std::cout << "[DEBUG:MARCH] No initial CI found. Exiting.\n";
        return iface;
    }

    // First interface point from initial column search
    iface.push_back({init.r, init.z_ci, init.elem,
                     init.cathode_ci, init.wall_ci});

    // Marching state
    int    elem_i = init.elem;
    double r_i    = init.r;
    double z_i    = init.z_ci;

    double smallest_r_reached = r_i;

    // Cathode detection state
    bool   cathode_seen   = init.cathode_ci;
    double z_first_cath   = init.cathode_ci ? init.z_ci
                                            : std::numeric_limits<double>::infinity();

    const int   max_steps = 200;
    const double eps_r    = 1e-9;
    const double eps_z    = 1e-9;

    // ---------------------------------------------------------------------
    // 2) Main adjacency-based march: move to CI neighbors that are down-left.
    //    Use unified adjacency (face + vertex neighbors).
    //    Stop as soon as we see the first cathode hit.
    // ---------------------------------------------------------------------
    if (!cathode_seen) // if initial column already saw cathode, skip marching
    {
        for (int step = 0; step < max_steps; ++step)
        {
            if (debug)
            {
                std::cout << "\n[DEBUG:MARCH] --- Step " << step
                          << " at elem=" << elem_i
                          << " r_i=" << r_i
                          << " z_i=" << z_i << "\n";
            }

            // Stop near inner radius for wall-branch
            if (r_i <= r_min + geom_tol)
            {
                if (debug)
                    std::cout << "[DEBUG:MARCH] Reached inner radius stop (wall marching).\n";
                break;
            }

            // Use unified adjacency instead of face-only neighbors
            const auto &nbrs = unified_adj[elem_i];

            if (debug)
            {
                std::cout << "[DEBUG:MARCH] Unified neighbors of elem " << elem_i
                          << ": " << nbrs.size() << "\n";
            }

            bool   found_ci_here   = false;
            double best_z_ci       = -std::numeric_limits<double>::infinity();
            int    best_elem       = -1;
            bool   cathode_here    = false;
            bool   wall_here       = false;

            bool   trigger_cathode = false; // becomes true on first cathode hit

            // Scan neighbors; pick those whose centroids lie "down-left"
            for (int nb : nbrs)
            {
                if (nb < 0 || nb >= mesh.GetNE())
                    continue;

                const double r_n = geom.r_centroid[nb];
                const double z_n = geom.z_centroid[nb];

                // down-left filter (allow tiny tolerance)
                if (r_n > r_i + eps_r) continue;
                if (z_n > z_i + eps_z) continue;

                if (debug)
                {
                    std::cout << "[DEBUG:MARCH]   neighbor elem=" << nb
                              << " r_n=" << r_n
                              << " z_n=" << z_n << " (candidate)\n";
                }

                double rr = r_n;
                double zz = z_n;

                ElectronTraceResult res;
                int elem_guess = nb;

                bool ok = TraceSingleCIProbe(ctx, rr, zz, elem_guess, res);

                if (debug)
                {
                    std::cout << "[DEBUG:MARCH]      probe ok=" << ok
                              << " exit=" << int(res.exit_code)
                              << " cathode=" << IsCathodeInsensitive(res)
                              << " wall=" << IsWallInsensitive(res)
                              << "\n";
                }

                if (!ok) continue;
                if (!IsChargeInsensitive(res)) continue;

                // First time we see a cathode hit anywhere → trigger cathode phase
                if (IsCathodeInsensitive(res) && !cathode_seen)
                {
                    cathode_seen   = true;
                    z_first_cath   = z_n;
                    trigger_cathode = true;

                    if (debug)
                    {
                        std::cout << "[DEBUG:MARCH] First cathode hit detected at "
                                  << "r=" << r_n << " z=" << z_n
                                  << " → stopping wall marching and entering cathode phase.\n";
                    }
                    break; // break neighbor loop
                }

                // Wall-based CI candidate
                // Keep highest z among CI neighbors
                if (!found_ci_here || z_n > best_z_ci)
                {
                    found_ci_here = true;
                    best_z_ci     = z_n;
                    best_elem     = elem_guess;
                    cathode_here  = IsCathodeInsensitive(res);
                    wall_here     = IsWallInsensitive(res);
                }
                else
                {
                    // OR-aggregate CI type flags
                    cathode_here = cathode_here || IsCathodeInsensitive(res);
                    wall_here    = wall_here    || IsWallInsensitive(res);
                }
            } // neighbor loop

            if (trigger_cathode)
            {
                // We just detected the first cathode hit; stop wall marching.
                break;
            }

            if (!found_ci_here || best_elem < 0)
            {
                if (debug)
                    std::cout << "[DEBUG:MARCH] No CI neighbor found → terminating wall march\n";
                break;
            }

            // Accept new wall interface point
            elem_i = best_elem;
            r_i    = geom.r_centroid[best_elem];
            z_i    = best_z_ci;

            iface.push_back({r_i, z_i, elem_i, cathode_here, wall_here});

            smallest_r_reached = std::min(smallest_r_reached, r_i);
        } // wall-marching loop
    } // if (!cathode_seen)

    // ---------------------------------------------------------------------
    // 3) Cathode plane construction (triggered as soon as a cathode hit is seen)
    // ---------------------------------------------------------------------
    if (cathode_seen)
    {
        if (debug)
        {
            std::cout << "[DEBUG:MARCH] Entering cathode plane construction.\n";
            std::cout << "[DEBUG:MARCH] z_first_cath (initial guess)=" << z_first_cath << "\n";
        }

        double z_cath = z_first_cath;
        if (!std::isfinite(z_cath))
        {
            // Fallback: if for some reason we didn't record z_first_cath,
            // take the last interface z or fall back to z_min.
            z_cath = (iface.empty() ? z_min : iface.back().z_ci);
        }

        double r_plane_start = smallest_r_reached;

        if (!iface.empty())
        {
            // Find index of smallest-r point
            int idx_min = 0;
            double r_min_iface = iface[0].r;
            for (int i = 1; i < (int)iface.size(); ++i)
            {
                if (iface[i].r < r_min_iface)
                {
                    r_min_iface = iface[i].r;
                    idx_min     = i;
                }
            }

            // Update that point to lie on the cathode plane
            iface[idx_min].z_ci       = z_cath;
            iface[idx_min].cathode_ci = true;
            iface[idx_min].wall_ci    = false;

            // Add an additional plane point at inner radius (if distinct in r)
            const double r_inner = r_min + geom_tol;
            if (r_inner < r_min_iface - 1e-12)
            {
                InterfacePoint p;
                p.r          = r_inner;
                p.z_ci       = z_cath;
                p.elem       = iface[idx_min].elem;
                p.cathode_ci = true;
                p.wall_ci    = false;
                iface.push_back(p);
            }
        }
        else
        {
            // If somehow the interface has no points, add just a cathode plane segment
            const double r_inner = r_min + geom_tol;
            InterfacePoint p1{r_plane_start, z_cath, -1, true, false};
            InterfacePoint p2{r_inner,       z_cath, -1, true, false};
            iface.push_back(p1);
            iface.push_back(p2);
        }
    }

    // ---------------------------------------------------------------------
    // 4) Sort by radius and print summary
    // ---------------------------------------------------------------------
    std::sort(iface.begin(), iface.end(),
              [](const InterfacePoint &a, const InterfacePoint &b)
              {
                  return a.r < b.r;
              });

    if (debug)
    {
        std::cout << "[DEBUG:MARCH] Completed marching. Points = "
                  << iface.size() << "\n";
        for (const auto &p : iface)
        {
            std::cout << "  r=" << p.r
                      << " z_ci=" << p.z_ci
                      << " elem=" << p.elem
                      << " cath=" << p.cathode_ci
                      << " wall=" << p.wall_ci << "\n";
        }
    }

    return iface;
}


double ComputeCIV_Marching(const SimulationResult &sim,
                           const Config           &cfg)
{
    using namespace mfem;

    MFEM_VERIFY(sim.mesh, "ComputeCIV_Marching: mesh is null");
    ParMesh &global_mesh = *sim.mesh;

    MFEM_VERIFY(sim.V, "ComputeCIV_Marching: potential V is null");
    GridFunction &global_phi = *sim.V;

    const FiniteElementSpace      *fes_phi = global_phi.FESpace();
    const FiniteElementCollection *fec     = fes_phi->FEColl();
    const int                      ordering = fes_phi->GetOrdering();

    // 1) Build shared geometry/connectivity ONCE
    ElementAdjacency           adj         = BuildAdjacency(global_mesh);
    ElementGeometryData        geom        = ComputeElementGeometry(global_mesh);
    std::vector<std::vector<int>> vertex_adj = BuildVertexAdjacency(global_mesh);
    std::vector<std::vector<int>> unified_adj =
        BuildUnifiedAdjacency(adj, vertex_adj, global_mesh.GetNE());

    // 2) Build a single base TracerContext
    TracerContext base_ctx(global_mesh,
                           fec,
                           ordering,
                           global_phi,
                           adj,
                           geom,
                           unified_adj,
                           cfg);

    // 3) Build interface using base_ctx (it will internally clone per-thread)
    std::vector<InterfacePoint> iface =
        BuildInterfaceCurve_Marching(base_ctx);

    if (iface.empty()) {
        if (cfg.debug.debug)
            std::cout << "[DEBUG:CIV] No interface points found; CIV = 0.\n";
        return 0.0;
    }

    std::sort(iface.begin(), iface.end(),
              [](const InterfacePoint &a, const InterfacePoint &b)
              { return a.r < b.r; });

    if (cfg.tracing_params.dump_civ_boundary) {
        namespace fs = std::filesystem;
        fs::path outdir(cfg.save_path);
        fs::path fname = outdir / "civ_interface_boundary.csv";
        std::ofstream ofs(fname);
        ofs << "# r, z_ci\n";
        for (const auto &ip : iface)
            ofs << ip.r << "," << ip.z_ci << "\n";
        if (cfg.debug.debug)
            std::cout << "[DEBUG:CIV] Dumped CIV boundary to: "
                      << fname.string() << "\n";
    }

    // 4) Volume classification uses the same geom we already have
    CIVolumeResult vres = ComputeCIVFromInterface(global_mesh, cfg, geom, iface);
    double V_total = vres.V_total;
    double V_CI    = vres.V_CI;
    double frac    = (V_total > 0.0) ? (V_CI / V_total) : 0.0;

    if (cfg.debug.debug) {
        std::cout << "\n----- Charge Insensitive Volume (Marching) -----\n";
        std::cout << "Total Volume considered: " << V_total << "\n";
        std::cout << "CIV Volume:              " << V_CI    << "\n";
        std::cout << "CIV Fraction:            " << frac    << "\n";
        std::cout << "------------------------------------------------\n";
    }

    return frac;
}
// Used to do simple downwards integration of the resulting boundaries instead of 
// summing mesh elements
static double ComputeCIV_FromInterfaceAnalytic(const Config                  &cfg,
                                               const std::vector<InterfacePoint> &iface)
{
    if (iface.empty())
    {
        return 0.0;
    }

    TpcGeometry tpc(cfg);

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
    if (A_tot <= 0.0)
    {
        return 0.0;
    }

    const double frac = A_ci / A_tot;
    return frac;
}
// ============================================================================
// CIV by radial column sweeps (fixed r-slices, z-search per slice)
// ============================================================================
static double ComputeCIV_ColumnSweep(const SimulationResult &sim,
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
    ElementGeometryData           geom        = ComputeElementGeometry(global_mesh);
    std::vector<std::vector<int>> vertex_adj  = BuildVertexAdjacency(global_mesh);
    std::vector<std::vector<int>> unified_adj =
        BuildUnifiedAdjacency(adj, vertex_adj, global_mesh.GetNE());

    // Base tracer context
    TracerContext base_ctx(global_mesh,
                           fec,
                           ordering,
                           global_phi,
                           adj,
                           geom,
                           unified_adj,
                           cfg);

    // Speed up by not requiring to traverse the entire TPC 

    // Geometry limits and params
    TpcGeometry tpc(cfg);
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
    const int n_r_slices = 25; // TODO Config entry 

    // Radial range we actually sample (avoid exact boundaries)
    const double r_lo = r_min + geom_tol;
    // Really dont need to start high
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
        std::cout << "[DEBUG:CIV] n_r_slices=" << n_r_slices
                  << " threads=" << n_threads
                  << " max_probes_per_column=" << max_probes_per_column << "\n";
    }

    // Build per-thread contexts (same pattern as in FindInitialColumnCI)
    int n_used_threads = std::max(1, n_threads);
    std::vector<std::unique_ptr<TracerContext>> extra;
    extra.reserve(n_used_threads - 1);

    std::vector<TracerContext*> ctx_pool(n_used_threads);
    ctx_pool[0] = &base_ctx;

    for (int t = 1; t < n_used_threads; ++t)
    {
        extra.emplace_back(std::make_unique<TracerContext>(
            base_ctx.mesh, fec, ordering, base_ctx.phi,
            adj, geom, unified_adj, cfg));
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

    if (cfg.tracing_params.dump_civ_boundary)
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
    // THis is too slow i dont want to parallelize it we are doing simple integration
    // Use the same volume classification as with the marching interface
    // CIVolumeResult vres = ComputeCIVFromInterface(global_mesh, cfg, geom, iface);
    double V_CI = ComputeCIV_FromInterfaceAnalytic(cfg, iface);

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

    const std::string &method = cfg.tracing_params.method; // "InformedSweep" or "RandomSample"

    if (method == "RandomSample")
    {
        return ComputeCIV_RandomSample(cfg, result);
    }
    else if (method == "InformedSweep")
    {
        // Marching-interface CIV
        if (cfg.debug.debug)
        {
            std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (InformedSweep / marching)" << std::endl;
        }
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

