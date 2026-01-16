#include "boundary_conditions.h"


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
void ApplyNeumannValues(ParLinearForm &b, const Array<int> &neumann_attr, const std::shared_ptr<const Config>& cfg, mfem::Coefficient &w, std::vector<std::unique_ptr<mfem::Coefficient>> &owned_coeffs, std::vector<mfem::Array<int>> &owned_markers)
{
  ParFiniteElementSpace *pfes = b.ParFESpace();
  MFEM_VERIFY(pfes != nullptr, "ApplyNeumannValues: ParLinearForm has no ParFESpace.");
  Mesh *mesh = pfes->GetMesh();

  if (cfg->debug.printBoundaryConditions) { std::cout << "[DEBUG:NeumannBC]\n"; }

  // Iterate config elements with neumann attributes
  for (const auto &[name, bc] : cfg->boundaries)
  {
    if (bc.type != "neumann") { continue; }

    owned_markers.emplace_back(MakeBdrMarker(mesh, {bc.bdr_id}));
    mfem::Array<int> &marker = owned_markers.back();

    // If this MPI rank has no boundary attributes, or this attribute is absent, skip.
    if (marker.Size() == 0) { continue; }
    if (marker.Max() == 0)  { continue; }
    // Shouldnt be higher 
    if (marker.Max() != 1) { std::cerr << "\033[33m" << "WARNING: boundary " << name << "'(bdr_id= " << bc.bdr_id << ")' not present in mesh!\n" << "\033[0m\n"; }
    // Default branch - Constant charge density on boundary 
    if (!bc.depth_dependent)
    {
      auto g  = std::make_unique<ConstantCoefficient>(bc.value);
      auto wg = std::make_unique<ProductCoefficient>(w, *g);  // axisym weight

      b.AddBoundaryIntegrator(new BoundaryLFIntegrator(*wg), marker);
      
      owned_coeffs.emplace_back(std::move(g));
      owned_coeffs.emplace_back(std::move(wg));

      if (cfg->debug.printBoundaryConditions)
      {
        std::cout << "\t\t" << name << " '(bdr_id = " << bc.bdr_id
                  << ")' constant value=" << bc.value << "\n";
      }
    }
    else
    { // Depth dependent Neumann Boundary (per-boundary parameters assumed pre-validated)
      const double z_bot     = /* p.z_bot */     0.0;
      const double z_top     = /* p.z_top */     1.0;
      const double sigma_bot = /* p.sigma_bot */ -0.1e-6; // C/m^2
      const double sigma_top = /* p.sigma_top */ -0.5e-6; // C/m^2

      auto sigma_z = [=](const Vector &x) -> double
      {
        const double z = x[1];
        const double t = (z - z_bot) / (z_top - z_bot); // assume validated z_top != z_bot
        return sigma_bot + t * (sigma_top - sigma_bot);
      };

      auto g  = std::make_unique<FunctionCoefficient>(sigma_z);
      auto wg = std::make_unique<ProductCoefficient>(w, *g);
      b.AddBoundaryIntegrator(new BoundaryLFIntegrator(*wg), marker);
      
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

