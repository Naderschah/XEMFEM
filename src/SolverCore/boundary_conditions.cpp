#include "boundary_conditions.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>


// Build a safe boundary marker for the given attribute IDs.
static Array<int> MakeBdrMarker(const Mesh *mesh,
                                std::vector<int> attrs)
{
    mesh->bdr_attributes.HostRead(); // safe in device builds
    int max_id = mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0;
    Array<int> m(max_id); // On threaded size 0 is valid
    m = 0;
    for (int a : attrs)
    {
        if (mesh->bdr_attributes.Find(a) != -1) { m[a - 1] = 1; }
        else { std::cerr << "[BC] Warning: boundary attribute " << a
                         << " not found in mesh.\n"; }
    }
    return m;
}

// Identifies which boundaries should have BC conditions
BoundaryConditionGroups GetBoundaryConditionGroups(const mfem::Mesh *mesh,
                                                   const std::shared_ptr<const Config> &cfg)
{
    bool had_warning = false;
    auto warn = [&](const std::string &msg)
    {
        had_warning = true;
        std::cerr << msg << "\n";
    };

    // Make vector arrays
    std::vector<int> dir_ids;
    std::vector<int> neu_ids;
    for (const auto &[name, bc] : cfg->boundaries)
    {
        if (bc.bdr_id <= 0)
        {
            warn("[BC] Warning: boundary '" + name + "' has invalid bdr_id=" + std::to_string(bc.bdr_id));
            continue;
        }

        if      (bc.type == "dirichlet") { dir_ids.push_back(bc.bdr_id); }
        else if (bc.type == "neumann") { neu_ids.push_back(bc.bdr_id); }
        else if (bc.type == "robin") { throw std::runtime_error("[BC] robin boundary conditions are not implemented yet (boundary '" + name + "')."); }
        else { throw std::runtime_error("[BC] Unknown boundary condition type '" + bc.type + "' for boundary '" + name + "'."); }
    }

    // Check all are unique TODO Should happen on config parsing not here
    auto warn_overlap = [&](const char *a_name, const std::vector<int> &a,
                            const char *b_name, const std::vector<int> &b)
    {
        for (int id : a)
        {
            if (std::find(b.begin(), b.end(), id) != b.end())
            {
                warn(std::string("[BC] Warning: boundary attribute ") + std::to_string(id) +
                     " appears in both " + a_name + " and " + b_name +
                     " (Dirichlet should usually take precedence).");
            }
        }
    };
    warn_overlap("dirichlet", dir_ids, "neumann", neu_ids);

    if (had_warning)
    {
        throw std::runtime_error("[BC] Boundary condition configuration warnings occurred (see stderr).");
    }

    // Make Bdr markers
    BoundaryConditionGroups groups;
    groups.dirichlet_attr = MakeBdrMarker(mesh, dir_ids);
    groups.neumann_attr   = MakeBdrMarker(mesh, neu_ids);

    return groups;
}




void ApplyDirichletValues(GridFunction &V, const Array<int> &dirichlet_attr, const std::shared_ptr<const Config>& cfg)
{
  Mesh *mesh = V.FESpace()->GetMesh();
  if (cfg->debug.printBoundaryConditions) std::cout << "[DEBUG:DirichletBC] \n";
  // Iterate config elements with dirichlet attributes
  for (const auto& [name, bc] :cfg->boundaries)
  {
    if (bc.type != "dirichlet") continue; 

    Array<int> marker = MakeBdrMarker(mesh, {bc.bdr_id});
    if (marker.Max() == 1)
    {
      ConstantCoefficient Vcoef(bc.value);
      V.ProjectBdrCoefficient(Vcoef, marker);
      if (cfg->debug.printBoundaryConditions) { std::cout << "\t\t" << name << " '(bdr_id = " << bc.bdr_id << ")'\n"; }
    }
    else 
    { std::cerr << "\033[33m" << "WARNING: boundary " << name << "'(bdr_id= " << bc.bdr_id << ")' not present in mesh!\n" << "\033[0m\n"; }
  }
}

// Depth dependent charge profile
static std::unique_ptr<mfem::Coefficient> MakeLinearZProfileCoeff(double z_bot, double z_top, double sigma_bot, double sigma_top)
{
    auto f = [=](const Vector &x) -> double
    {
        const double z = x[1];
        const double t = (z - z_bot) / (z_top - z_bot);
        return sigma_bot + t * (sigma_top - sigma_bot);
    };
    return std::make_unique<FunctionCoefficient>(f);
}

static void ValidateAndReportWallChargeBoundary(
    mfem::ParFiniteElementSpace &pfes,
    const Config &cfg,
    const Boundary &bc,
    mfem::Coefficient &axisymmetric_sigma)
{
    mfem::Mesh &mesh = *pfes.GetMesh();
    MPI_Comm comm = pfes.GetComm();
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    long long local_elements = 0;
    double local_r_min = std::numeric_limits<double>::infinity();
    double local_r_max = -std::numeric_limits<double>::infinity();
    double local_z_min = std::numeric_limits<double>::infinity();
    double local_z_max = -std::numeric_limits<double>::infinity();
    double local_integrated_charge = 0.0;

    for(int be = 0; be < mesh.GetNBE(); ++be)
    {
        if(mesh.GetBdrAttribute(be) != bc.bdr_id) continue;
        ++local_elements;

        mfem::Array<int> vertices;
        mesh.GetBdrElementVertices(be, vertices);
        for(int i = 0; i < vertices.Size(); ++i)
        {
            const double *point = mesh.GetVertex(vertices[i]);
            local_r_min = std::min(local_r_min, point[0]);
            local_r_max = std::max(local_r_max, point[0]);
            local_z_min = std::min(local_z_min, point[1]);
            local_z_max = std::max(local_z_max, point[1]);
        }

        mfem::ElementTransformation *T =
            mesh.GetBdrElementTransformation(be);
        MFEM_VERIFY(T != nullptr,
                    "BC_WallCharge has no boundary element transformation.");
        const mfem::Geometry::Type geometry =
            mesh.GetBdrElement(be)->GetGeometryType();
        const mfem::IntegrationRule &rule =
            mfem::IntRules.Get(geometry, 4);
        for(int q = 0; q < rule.GetNPoints(); ++q)
        {
            const mfem::IntegrationPoint &ip = rule.IntPoint(q);
            T->SetIntPoint(&ip);
            local_integrated_charge +=
                ip.weight * T->Weight() * axisymmetric_sigma.Eval(*T, ip);
        }
    }

    long long global_elements = 0;
    double global_r_min = 0.0;
    double global_r_max = 0.0;
    double global_z_min = 0.0;
    double global_z_max = 0.0;
    double global_integrated_charge = 0.0;
    MPI_Allreduce(&local_elements, &global_elements, 1,
                  MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&local_r_min, &global_r_min, 1,
                  MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_r_max, &global_r_max, 1,
                  MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&local_z_min, &global_z_min, 1,
                  MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_z_max, &global_z_max, 1,
                  MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&local_integrated_charge, &global_integrated_charge, 1,
                  MPI_DOUBLE, MPI_SUM, comm);

    MFEM_VERIFY(global_elements > 0,
                "BC_WallCharge has no boundary elements after MFEM mesh loading.");

    const double scale = std::max(
        {1.0, std::abs(global_r_min), std::abs(global_r_max),
         std::abs(global_z_min), std::abs(global_z_max)});
    const double coordinate_tolerance = std::max(1e-8, 1e-7 * scale);
    MFEM_VERIFY(std::abs(global_r_max - global_r_min) <= coordinate_tolerance,
                "BC_WallCharge is not a constant-radius surface after MFEM loading: r=["
                << global_r_min << ", " << global_r_max << "].");
    MFEM_VERIFY(std::abs(global_z_min - bc.z_bot) <= coordinate_tolerance &&
                std::abs(global_z_max - bc.z_top) <= coordinate_tolerance,
                "BC_WallCharge z extent does not match its charge profile: mesh z=["
                << global_z_min << ", " << global_z_max
                << "], configured z=[" << bc.z_bot << ", " << bc.z_top << "].");

    if(cfg.solver.axisymmetric)
    {
        const double radius = 0.5 * (global_r_min + global_r_max);
        const double expected_charge =
            2.0 * M_PI * radius * (bc.z_top - bc.z_bot) *
            0.5 * (bc.value_bot + bc.value_top);
        const double charge_tolerance =
            std::max(1e-15, 1e-9 * std::abs(expected_charge));
        MFEM_VERIFY(
            std::abs(global_integrated_charge - expected_charge) <= charge_tolerance,
            "BC_WallCharge assembled charge does not match the configured linear profile: "
            << "integrated=" << global_integrated_charge
            << " C, expected=" << expected_charge << " C.");

        if(rank == 0)
        {
            std::cout << std::setprecision(17)
                      << "[BC_WallCharge validation] bdr_id=" << bc.bdr_id
                      << " elements=" << global_elements
                      << " r=" << radius
                      << " z=[" << global_z_min << ", " << global_z_max << "]"
                      << " sigma=[" << bc.value_bot << ", " << bc.value_top << "]"
                      << " integrated_charge_C=" << global_integrated_charge
                      << " expected_charge_C=" << expected_charge << "\n";
        }
    }
}

void ApplyNeumannValues(ParLinearForm &b,
                        const Array<int> & /*neumann_attr*/,
                        const std::shared_ptr<const Config> &cfg,
                        mfem::Coefficient &w,
                        std::vector<std::unique_ptr<mfem::Coefficient>> &owned_coeffs,
                        std::vector<mfem::Array<int>> &owned_markers)
{
    ParFiniteElementSpace *pfes = b.ParFESpace();
    MFEM_VERIFY(pfes != nullptr, "ApplyNeumannValues: ParLinearForm has no ParFESpace.");
    Mesh *mesh = pfes->GetMesh();

    if (cfg->debug.printBoundaryConditions) { std::cout << "[DEBUG:NeumannBC]\n"; }

    // ---- Reserve marker storage to keep MFEM's internal marker pointers valid ----
    size_t n_neumann = 0;
    for (const auto &[name, bc] : cfg->boundaries)
    {
        if (bc.type == "neumann") { ++n_neumann; }
    }
    owned_markers.reserve(owned_markers.size() + n_neumann);
    // ---------------------------------------------------------------------------

    for (const auto &[name, bc] : cfg->boundaries)
    {
        if (bc.type != "neumann") { continue; }

        owned_markers.emplace_back(MakeBdrMarker(mesh, {bc.bdr_id}));
        mfem::Array<int> &marker = owned_markers.back();

        const bool local_marker_present =
            marker.Size() > 0 && marker.Max() == 1;
        if (name == "BC_WallCharge")
        {
            // Wall validation contains collectives, so every rank must enter
            // it even when a partition owns no wall segments locally.
            const int local_present = local_marker_present ? 1 : 0;
            int global_present = 0;
            MPI_Allreduce(&local_present, &global_present, 1,
                          MPI_INT, MPI_MAX, pfes->GetComm());
            MFEM_VERIFY(global_present == 1,
                        "BC_WallCharge is absent from the distributed mesh.");
        }
        else if (!local_marker_present)
        {
            continue;
        }

        // Shouldn't be higher
        if (marker.Size() > 0 && marker.Max() > 1)
        {
            std::cerr << "\033[33m"
                      << "WARNING: boundary " << name << "'(bdr_id= " << bc.bdr_id
                      << ")' not present in mesh!\n"
                      << "\033[0m\n";
        }

        if (!bc.depth_dependent)
        {
            auto g  = std::make_unique<ConstantCoefficient>(bc.value);
            auto wg = std::make_unique<ProductCoefficient>(w, *g);  // axisym weight

            b.AddBoundaryIntegrator(new BoundaryLFIntegrator(*wg), marker);

            if (name == "BC_WallCharge")
            {
                Boundary zero_charge = bc;
                zero_charge.value_bot = bc.value;
                zero_charge.value_top = bc.value;
                ValidateAndReportWallChargeBoundary(
                    *pfes, *cfg, zero_charge, *wg);
            }

            owned_coeffs.emplace_back(std::move(g));
            owned_coeffs.emplace_back(std::move(wg));

            if (cfg->debug.printBoundaryConditions)
            {
                std::cout << "\t\t" << name << " '(bdr_id = " << bc.bdr_id
                          << ")' constant value=" << bc.value << "\n";
            }
        }
        else
        {
            const double z_bot     = bc.z_bot;
            const double z_top     = bc.z_top;
            const double sigma_bot = bc.value_bot;
            const double sigma_top = bc.value_top;

            auto sigma_z = [=](const Vector &x) -> double
            {
                const double z = x[1];
                const double t = (z - z_bot) / (z_top - z_bot); // assume validated z_top != z_bot
                return sigma_bot + t * (sigma_top - sigma_bot);
            };

            auto g  = std::make_unique<FunctionCoefficient>(sigma_z);
            auto wg = std::make_unique<ProductCoefficient>(w, *g);

            b.AddBoundaryIntegrator(new BoundaryLFIntegrator(*wg), marker);

            if (name == "BC_WallCharge")
            {
                ValidateAndReportWallChargeBoundary(
                    *pfes, *cfg, bc, *wg);
            }

            owned_coeffs.emplace_back(std::move(g));
            owned_coeffs.emplace_back(std::move(wg));

            if (cfg->debug.printBoundaryConditions)
            {
                std::cout << "\t\t" << name << " '(bdr_id = " << bc.bdr_id
                          << ")' depth-dependent sigma(z): " << sigma_bot
                          << " -> " << sigma_top << " C/m^2"
                          << " over z=[" << z_bot << ", " << z_top << "]\n";
            }
        }
    }
}

