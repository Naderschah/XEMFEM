/*
We have 2 options
- Interpolate a grid
- Integrate on the mesh <- Doing this 

Little GPT summary on the error propagation here

Field-Line Tracing Error Model
------------------------------

We trace trajectories of the discrete electric field E_h defined on the FEM mesh.
The exact field is E, the discrete field is E_h, and the numerical trajectory uses
a time step Δt.

1. Field Approximation Error
   Let x(t) solve ẋ = −E(x), and let x_h(t) solve ẋ = −E_h(x).
   If  ||E − E_h||_{L∞} ≤ C h^m  and E_h is Lipschitz with constant L, then

      ||x_h(t) − x(t)|| ≤ C_T ||E − E_h||_{L∞} ≤ C_T' h^m ,

   so field lines of E_h converge to those of E as the mesh is refined.

2. Time-Discretization Error
   Let x̃_h^n be the numerical solution of order p with step Δt. Then

      ||x̃_h(t_n) − x_h(t_n)|| ≤ C_T Δt^p ,

   giving total trajectory error

      ||x̃_h(t_n) − x(t_n)|| ≤ C_1 h^m + C_2 Δt^p .

3. Step-Size Condition
   We enforce

      Δt ≤ c h_K / ||v_h(x)|| ,

   where h_K is the size of the current element. Then

      ||x_{n+1} − x_n|| ≤ c h_K ,

   which guarantees x_{n+1} lies in the current element K or one of its
   immediate neighbors (shape-regularity). This ensures neighbor-only
   element lookup is exact and introduces no additional modeling error.

4. Comparison with Grid Interpolation
   If the field is additionally interpolated onto a structured grid
   (spacing H, interpolation order q), producing E_g, then

      ||x_g(t) − x(t)|| ≤ C_T ( h^m + H^q ),

   so grid-based tracing adds an extra spatial error term H^q that the
   direct FEM-based tracer does not incur.


------------------------------
Functions and logic 

class: ElectricFieldCoeff
- Just a wrapper to get the vector valued field at the point 



*/


#include <mfem.hpp>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "trace_fieldlines.h"

// Wrapper class to simplify calls 
ElectricFieldCoeff::ElectricFieldCoeff(const mfem::GridFunction &E)
    : mfem::VectorCoefficient(E.FESpace()->GetMesh()->SpaceDimension())
    , E_(E)
{ }

void ElectricFieldCoeff::Eval(mfem::Vector &V,
                              mfem::ElementTransformation &T,
                              const mfem::IntegrationPoint &ip)
{
    E_.GetVectorValue(T.ElementNo, ip, V);
}

// Get adjacency of each mesh element - pre iter constructor
static ElementAdjacency BuildAdjacencyAndActiveMask(
    const SimulationResult &sim,
    const Config &cfg)
{
    using namespace mfem;

    ElementAdjacency adj;

    ParMesh &mesh = *sim.mesh;
    const int ne  = mesh.GetNE();
    const int dim = mesh.SpaceDimension();

    adj.neighbors.resize(ne);
    adj.h.resize(ne);
    adj.active.assign(ne, false);

    const double r_min = 0.0;
    const double r_max = cfg.mesh.tpc_r;
    const double z_min = cfg.mesh.z_cathode;
    const double z_max = cfg.mesh.z_liquidgas;

    Array<int> vert_ids;

    // Build element-to-element connectivity once
    const Table &el_to_el = mesh.ElementToElementTable();

    for (int e = 0; e < ne; ++e)
    {
        // ---------------- Neighbors from ElementToElementTable ----------------
        {
            Array<int> nbrs;
            el_to_el.GetRow(e, nbrs);
            adj.neighbors[e].assign(nbrs.begin(), nbrs.end());
        }

        // ---------------- Element size via bounding box ----------------
        mesh.GetElementVertices(e, vert_ids);
        const int nv = vert_ids.Size();

        Vector x_min(dim), x_max(dim);
        x_min = 0.0;
        x_max = 0.0;

        for (int j = 0; j < nv; ++j)
        {
            const int vid = vert_ids[j];
            const double *coord = mesh.GetVertex(vid);

            for (int d = 0; d < dim; ++d)
            {
                const double val = coord[d];
                if (j == 0)
                {
                    x_min[d] = x_max[d] = val;
                }
                else
                {
                    if (val < x_min[d]) x_min[d] = val;
                    if (val > x_max[d]) x_max[d] = val;
                }
            }
        }

        double dx2 = 0.0;
        for (int d = 0; d < dim; ++d)
        {
            const double diff = x_max[d] - x_min[d];
            dx2 += diff * diff;
        }
        adj.h[e] = std::sqrt(dx2);

        // ---------------- Active flag (center inside TPC volume) ----------------
        Vector ctr(dim);
        for (int d = 0; d < dim; ++d)
        {
            ctr[d] = 0.5 * (x_min[d] + x_max[d]);
        }

        const double r = ctr[0];
        const double z = (dim > 1) ? ctr[1] : 0.0;

        if (r >= r_min && r <= r_max &&
            z >= z_min && z <= z_max)
        {
            adj.active[e] = true;
        }
    }

    return adj;
}


// Look through current element and its neighbors
static bool FindElementForPointLocal(
    mfem::ParMesh      &mesh,
    const ElementAdjacency   &adj,
    int                       current_elem,
    const mfem::Vector       &x_new,
    int                      &out_elem,
    mfem::IntegrationPoint   &out_ip)
{
    // Try current element
    {
        mfem::ElementTransformation *T = mesh.GetElementTransformation(current_elem);
        mfem::IntegrationPoint ip;
        const int success = T->TransformBack(x_new, ip);
        if (success)
        {
            out_elem = current_elem;
            out_ip   = ip;
            return true;
        }
    }

    // Try neighbors.
    const auto &nbrs = adj.neighbors[current_elem];
    for (int nb : nbrs)
    {
        if (nb < 0 || nb >= mesh.GetNE()) { continue; }

        mfem::ElementTransformation *T = mesh.GetElementTransformation(nb);
        mfem::IntegrationPoint ip;
        const int success = T->TransformBack(x_new, ip);
        if (success)
        {
            out_elem = nb;
            out_ip   = ip;
            return true;
        }
    }

    return false;
}

static ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh        &mesh,
    ElectricFieldCoeff   &E_coeff,
    const ElementAdjacency     &adj,
    const Config               &cfg,
    int                         start_elem,
    const mfem::IntegrationPoint &start_ip,
    const ElectronTraceParams  &params)
{
    using namespace mfem;

    ElectronTraceResult result;

    const double r_min = 0.0;
    const double r_max = cfg.mesh.tpc_r;
    const double z_min = cfg.mesh.z_cathode;
    const double z_max = cfg.mesh.z_liquidgas;

    const int dim = mesh.SpaceDimension();

    int elem_id = start_elem;
    IntegrationPoint ip = start_ip;

    double t = 0.0;
    int step = 0;

    while (step < params.max_steps && t < params.max_time)
    {
        ElementTransformation *T = mesh.GetElementTransformation(elem_id);
        T->SetIntPoint(&ip);
        Vector x(dim);
        T->Transform(ip, x);  // x = (r,z)

        Vector E(dim);
        E_coeff.Eval(E, *T, ip);

        const double Enorm = E.Norml2();
        if (Enorm < params.min_field_norm)
        {
            result.exit_code = ElectronExitCode::WeakField;
            result.points.push_back(x);
            result.exit_element = elem_id;
            break;
        }

        Vector v(E);
        v *= -1.0;  // electrons drift opposite to E

        const double hK = adj.h[elem_id];
        double dt = params.c_step * hK / Enorm;

        if (dt <= 0.0)
        {
            result.exit_code = ElectronExitCode::WeakField;
            result.points.push_back(x);
            result.exit_element = elem_id;
            break;
        }

        if (t + dt > params.max_time)
        {
            dt = params.max_time - t;
        }

        Vector x_new(dim);
        x_new = x;
        x_new.Add(dt, v);

        const double r = x_new[0];
        const double z = x_new[1];

        // physical boundaries
        if (r <= params.geom_tol)
        {
            result.exit_code = ElectronExitCode::HitAxis;
            result.points.push_back(x_new);
            result.exit_element = elem_id;
            break;
        }
        if (r >= r_max - params.geom_tol)
        {
            result.exit_code = ElectronExitCode::HitWall;
            result.points.push_back(x_new);
            result.exit_element = elem_id;
            break;
        }
        if (z <= z_min + params.geom_tol)
        {
            result.exit_code = ElectronExitCode::HitCathode;
            result.points.push_back(x_new);
            result.exit_element = elem_id;
            break;
        }
        if (z >= z_max - params.geom_tol)
        {
            result.exit_code = ElectronExitCode::HitLiquidGas;
            result.points.push_back(x_new);
            result.exit_element = elem_id;
            break;
        }

        int    elem_new = -1;
        IntegrationPoint ip_new;
        bool found = FindElementForPointLocal(mesh, adj, elem_id, x_new,
                                              elem_new, ip_new);

        if (!found || elem_new < 0 || elem_new >= mesh.GetNE() || !adj.active[elem_new])
        {
            result.exit_code = ElectronExitCode::LeftVolume;
            result.points.push_back(x_new);
            result.exit_element = elem_new;
            break;
        }

        result.points.push_back(x_new);
        elem_id = elem_new;
        ip      = ip_new;
        result.exit_element = elem_id;

        t += dt;
        ++step;
    }

    if (result.exit_code == ElectronExitCode::None)
    {
        result.exit_code = ElectronExitCode::MaxSteps;
    }

    return result;
}


void TraceElectronFieldLines(
    const SimulationResult                  &sim,
    const Config                            &cfg,
    const std::vector<mfem::Vector>         &seed_points,
    const std::vector<int>                  &seed_elements,
    const std::vector<mfem::IntegrationPoint> &seed_ips,
    const ElectronTraceParams               &params,
    std::vector<ElectronTraceResult>        &out_results)
{
    using namespace mfem;

    out_results.clear();

    ParMesh &mesh = *sim.mesh;
    const GridFunction &E = *sim.E;

    MFEM_VERIFY(seed_points.size() == seed_elements.size() &&
                seed_points.size() == seed_ips.size(),
                "TraceElectronFieldLines: seed sizes mismatch");

    ElementAdjacency adj = BuildAdjacencyAndActiveMask(sim, cfg);
    ElectricFieldCoeff E_coeff(E);

    const double r_min = 0.0;
    const double r_max = cfg.mesh.tpc_r;
    const double z_min = cfg.mesh.z_cathode;
    const double z_max = cfg.mesh.z_liquidgas;

    out_results.reserve(seed_points.size());

    for (std::size_t i = 0; i < seed_points.size(); ++i)
    {
        const Vector &seed = seed_points[i];
        ElectronTraceResult res;

        const double r = seed[0];
        const double z = seed[1];

        if (r < r_min || r > r_max || z < z_min || z > z_max)
        {
            res.exit_code = ElectronExitCode::InvalidSeed;
            res.points.push_back(seed);
            out_results.push_back(std::move(res));
            continue;
        }

        const int start_elem = seed_elements[i];
        const IntegrationPoint &start_ip = seed_ips[i];

        if (start_elem < 0 || start_elem >= mesh.GetNE() || !adj.active[start_elem])
        {
            res.exit_code = ElectronExitCode::LeftVolume;
            res.points.push_back(seed);
            out_results.push_back(std::move(res));
            continue;
        }

        res = TraceSingleElectronLine(mesh, E_coeff, adj, cfg,
                                      start_elem, start_ip, params);

        out_results.push_back(std::move(res));
    }
}

CivSeeds ExtractCivSeeds(const Config &cfg, const SimulationResult &result, int ir_order)
{
    using namespace mfem;

    MFEM_VERIFY(result.mesh, "ExtractCivSeeds: mesh is null");

    ParMesh &pmesh = *result.mesh;
    const int dim  = pmesh.Dimension();
    MFEM_VERIFY(dim == 2, "ExtractCivSeeds: only 2D axisymmetric meshes supported");

    const double r_tpc     = cfg.mesh.tpc_r;
    const double z_cathode = cfg.mesh.z_cathode;
    const double z_lgi     = cfg.mesh.z_liquidgas;

    IntegrationRules irs;

    CivSeeds seeds;

    const int ne = pmesh.GetNE();
    seeds.positions.reserve(ne);
    seeds.volumes.reserve(ne);
    seeds.elements.reserve(ne);
    seeds.ips.reserve(ne);

    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "ElementTransformation is null");

        Geometry::Type geom = T->GetGeometryType();
        const IntegrationRule &ir = irs.Get(geom, ir_order);

        // One representative IP per element (first one)
        for (int i = 0; i < 1; ++i)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            Vector x(dim);
            T->Transform(ip, x);
            const double r = x(0);
            const double z = x(1);

            // Restrict to TPC active volume
            if (r < 0.0 || r > r_tpc || z < z_cathode || z > z_lgi)
            {
                continue;
            }

            const double detJ  = T->Weight();   // |det(J)| in (r,z)
            const double w_ref = ip.weight;     // reference weight
            const double dV    = 2.0 * M_PI * r * detJ * w_ref;

            seeds.positions.push_back(x);
            seeds.volumes.push_back(dV);
            seeds.elements.push_back(e);

            // We must store a copy of the IP (not a reference)
            mfem::IntegrationPoint ip_copy = ip;
            seeds.ips.push_back(ip_copy);
        }
    }

    return seeds;
}