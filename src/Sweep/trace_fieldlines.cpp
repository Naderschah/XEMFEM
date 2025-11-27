#include <mfem.hpp>
#include <vector>
#include <cmath>

#include "trace_fieldlines.h"

using namespace mfem;

// ============================================================================
// Internal helpers
// ============================================================================

// Electric field coefficient wrapper around a vector GridFunction
struct ElectricFieldCoeff : public mfem::VectorCoefficient
{
    explicit ElectricFieldCoeff(const mfem::GridFunction &E_)
        : mfem::VectorCoefficient(E_.FESpace()->GetMesh()->SpaceDimension())
        , E(E_)
    { }

    void Eval(mfem::Vector &V,
              mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override
    {
        E.GetVectorValue(T.ElementNo, ip, V);
    }

    const mfem::GridFunction &E;
};

// Mesh connectivity / size info used by the tracer
struct ElementAdjacency
{
    std::vector<std::vector<int>> neighbors; // face neighbors per element
    std::vector<double>           h;         // characteristic element size
};

// Geometry helper: all TPC-related geometry info and checks
namespace
{
struct TpcGeometry
{
    double r_min;
    double r_max;
    double z_min;
    double z_max;

    explicit TpcGeometry(const Config &cfg)
    {
        const auto &tp = cfg.tracing_params;
        r_min = tp.r_min;
        r_max = tp.r_max;
        z_min = tp.z_min;
        z_max = tp.z_max;
    }

    bool Inside(double r, double z) const
    {
        return (r >= r_min && r <= r_max &&
                z >= z_min && z <= z_max);
    }

    ElectronExitCode ClassifyBoundary(double r, double z, double geom_tol) const
    {
        if (r <= r_min + geom_tol)   { return ElectronExitCode::HitAxis; }
        if (r >= r_max - geom_tol)   { return ElectronExitCode::HitWall; }
        if (z <= z_min + geom_tol)   { return ElectronExitCode::HitCathode; }
        if (z >= z_max - geom_tol)   { return ElectronExitCode::HitLiquidGas; }
        return ElectronExitCode::None;
    }
};
} // anonymous namespace

// Build adjacency and characteristic size h per element
static ElementAdjacency BuildAdjacency(const SimulationResult &sim)
{
    ElementAdjacency adj;

    ParMesh &mesh = *sim.mesh;
    const int ne  = mesh.GetNE();
    const int dim = mesh.SpaceDimension();

    adj.neighbors.resize(ne);
    adj.h.resize(ne);

    Array<int> vert_ids;

    const Table &el_to_el = mesh.ElementToElementTable();

    for (int e = 0; e < ne; ++e)
    {
        // --- neighbors
        {
            Array<int> nbrs;
            el_to_el.GetRow(e, nbrs);
            adj.neighbors[e].assign(nbrs.begin(), nbrs.end());
        }

        // --- element size via bounding box
        mesh.GetElementVertices(e, vert_ids);
        const int nv = vert_ids.Size();

        Vector x_min(dim), x_max(dim);
        x_min = 0.0;
        x_max = 0.0;

        for (int j = 0; j < nv; ++j)
        {
            const int     vid   = vert_ids[j];
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
    }

    return adj;
}

// Local element search: current element + its neighbors
static bool FindElementForPointLocal(
    mfem::ParMesh            &mesh,
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

    // Try neighbors
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

// Single electron tracing (internal helper)
namespace
{
ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh              &mesh,
    ElectricFieldCoeff         &E_coeff,
    const ElementAdjacency     &adj,
    const TpcGeometry          &geom,
    const ElectronTraceParams  &params,
    int                         start_elem,
    const mfem::IntegrationPoint &start_ip)
{
    ElectronTraceResult result;

    const int dim = mesh.SpaceDimension();

    int             elem_id = start_elem;
    IntegrationPoint ip     = start_ip;

    double t   = 0.0;
    int    step = 0;

    while (step < params.max_steps && t < params.max_time)
    {
        ElementTransformation *T = mesh.GetElementTransformation(elem_id);
        T->SetIntPoint(&ip);

        Vector x(dim);
        T->Transform(ip, x);  // x = (r,z)
        const double r = x[0];
        const double z = x[1];

        Vector E(dim);
        E_coeff.Eval(E, *T, ip);

        const double Enorm = E.Norml2();
        if (Enorm < params.min_field_norm)
        {
            result.exit_code    = ElectronExitCode::WeakField;
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
            result.exit_code    = ElectronExitCode::WeakField;
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

        const double r_new = x_new[0];
        const double z_new = x_new[1];

        // Geometry-based classification of boundary hits
        const ElectronExitCode bc =
            geom.ClassifyBoundary(r_new, z_new, params.geom_tol);

        if (bc != ElectronExitCode::None)
        {
            result.exit_code    = bc;
            result.points.push_back(x_new);
            result.exit_element = elem_id;
            break;
        }

        int             elem_new = -1;
        IntegrationPoint ip_new;
        bool found = FindElementForPointLocal(mesh, adj, elem_id, x_new,
                                              elem_new, ip_new);

        // If we cannot locate the new point in current element + neighbors,
        // treat it as "left volume" numerically.
        if (!found || elem_new < 0 || elem_new >= mesh.GetNE())
        {
            result.exit_code    = ElectronExitCode::LeftVolume;
            result.points.push_back(x_new);
            result.exit_element = elem_new;
            break;
        }

        result.points.push_back(x_new);
        elem_id         = elem_new;
        ip              = ip_new;
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
} // anonymous namespace

// ============================================================================
// Public API: seed extraction
// ============================================================================

CivSeeds ExtractCivSeeds(const Config &cfg, const SimulationResult &result)
{
    MFEM_VERIFY(result.mesh, "ExtractCivSeeds: mesh is null");

    ParMesh &pmesh = *result.mesh;
    const int dim  = pmesh.Dimension();
    MFEM_VERIFY(dim == 2, "ExtractCivSeeds: only 2D axisymmetric meshes supported");

    TpcGeometry geom(cfg);

    const int ir_order          = cfg.tracing_params.ir_order;
    const int max_seed_ip_per_e = cfg.tracing_params.integration_points;

    IntegrationRules irs;
    CivSeeds seeds;

    const int ne = pmesh.GetNE();
    seeds.positions.reserve(ne * max_seed_ip_per_e);
    seeds.volumes.reserve(ne * max_seed_ip_per_e);
    seeds.elements.reserve(ne * max_seed_ip_per_e);
    seeds.ips.reserve(ne * max_seed_ip_per_e);

    std::vector<double> element_volume(ne, 0.0);
    std::vector<int>    seeds_per_element(ne, 0);

    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "ExtractCivSeeds: ElementTransformation is null");

        Geometry::Type geom_type   = T->GetGeometryType();
        const IntegrationRule &ir  = irs.Get(geom_type, ir_order);
        const int n_ip             = ir.GetNPoints();
        const int n_seed_ip        = std::min(n_ip, max_seed_ip_per_e);

        double V_e = 0.0;

        for (int i = 0; i < n_ip; ++i)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            Vector x(dim);
            T->Transform(ip, x);
            const double r = x(0);
            const double z = x(1);

            // Contribution to full element volume (axisymmetric 3D)
            const double detJ  = T->Weight();   // |det(J)| in (r,z)
            const double w_ref = ip.weight;     // reference weight
            const double dV    = 2.0 * M_PI * r * detJ * w_ref;

            V_e += dV;

            // Only some IPs (up to n_seed_ip) are used as seeds
            if (i < n_seed_ip)
            {
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
        }

        element_volume[e] = V_e;
    }

    // Partition each element's volume among its seeds
    const int n_seeds = static_cast<int>(seeds.positions.size());
    seeds.volumes.resize(n_seeds);

    for (int j = 0; j < n_seeds; ++j)
    {
        const int e         = seeds.elements[j];
        const int n_in_elem = seeds_per_element[e];

        MFEM_VERIFY(n_in_elem > 0,
                    "ExtractCivSeeds: seed assigned to element with zero seed count");

        seeds.volumes[j] = element_volume[e] / static_cast<double>(n_in_elem);
    }

    return seeds;
}

// ============================================================================
// Public API: tracing
// ============================================================================

void TraceElectronFieldLines(const SimulationResult           &sim,
                             const Config                     &cfg,
                             const CivSeeds                   &seeds,
                             std::vector<ElectronTraceResult> &out_results)
{
    MFEM_VERIFY(sim.mesh, "TraceElectronFieldLines: mesh is null");
    MFEM_VERIFY(sim.E,    "TraceElectronFieldLines: electric field GridFunction is null");

    ParMesh &mesh = *sim.mesh;
    const int dim = mesh.SpaceDimension();
    MFEM_VERIFY(dim == 2, "TraceElectronFieldLines: only 2D axisymmetric meshes supported");

    const GridFunction &E = *sim.E;

    const std::size_t n_seeds = seeds.positions.size();
    MFEM_VERIFY(n_seeds == seeds.elements.size() &&
                n_seeds == seeds.ips.size(),
                "TraceElectronFieldLines: CivSeeds sizes mismatch");

    out_results.clear();
    out_results.reserve(n_seeds);

    ElementAdjacency adj   = BuildAdjacency(sim);
    ElectricFieldCoeff E_coeff(E);
    TpcGeometry geom(cfg);
    const ElectronTraceParams &params = cfg.tracing_params;

    for (std::size_t i = 0; i < n_seeds; ++i)
    {
        ElectronTraceResult res;

        const int start_elem = seeds.elements[i];
        MFEM_VERIFY(start_elem >= 0 && start_elem < mesh.GetNE(),
                    "TraceElectronFieldLines: seed element index out of range");

        const IntegrationPoint &start_ip = seeds.ips[i];

        res = TraceSingleElectronLine(mesh, E_coeff, adj, geom,
                                      params, start_elem, start_ip);

        out_results.push_back(std::move(res));
    }
}
