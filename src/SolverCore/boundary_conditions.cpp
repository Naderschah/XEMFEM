#include "boundary_conditions.h"


// Build a safe boundary marker for the given attribute IDs.
static Array<int> MakeBdrMarker(const Mesh *mesh,
                                std::vector<int> attrs)
{
    int max_id = mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0;
    Array<int> m(std::max(1, max_id)); // at least length 1, never 0
    m = 0;
    for (int a : attrs)
    {
        if (mesh->bdr_attributes.Find(a) != -1) { m[a - 1] = 1; }
        else { std::cerr << "[BC] Warning: boundary attribute " << a
                         << " not found in mesh.\n"; }
    }
    return m;
}

// Identifies which boundaries should have Dirichlet conditions
Array<int> GetDirichletAttributes(mfem::Mesh *mesh,
                                  const std::shared_ptr<const Config>& cfg)
{
    using namespace mfem;

    // Collect all Dirichlet boundary IDs from config
    std::vector<int> dirichlet_ids;
    dirichlet_ids.reserve(cfg->boundaries.size());
    for (const auto& [name, bc] : cfg->boundaries)
    {
        if (bc.type == "dirichlet")
            dirichlet_ids.push_back(bc.bdr_id);
    }
    // Convert to MFEM boundary marker array
    Array<int> ess = MakeBdrMarker(mesh, dirichlet_ids); // size = max bdr attr count

    // --- The ONLY axisymmetric modification: remove the axis boundary ID ---
    if (cfg->solver.axisymmetric)
    {
        int axis_id = cfg->solver.axisymmetric_r0_bd_attribute; // 1-based ID
        if (axis_id > 0 && axis_id <= ess.Size())
        {
            ess[axis_id - 1] = 0;   // ← just kill it
        }
    }

    return ess;
}



void ApplyDirichletValues(GridFunction &V, const Array<int> &dirichlet_attr, const std::shared_ptr<const Config>& cfg)
{
    Mesh *mesh = V.FESpace()->GetMesh();
    if (cfg->debug.debug) std::cout << "[DEBUG:DirichletBC] \n";
    // Iterate config elements with dirichlet attributes
    for (const auto& [name, bc] :cfg->boundaries)
    {
      if (bc.type != "dirichlet") continue; 

      Array<int> marker = MakeBdrMarker(mesh, {bc.bdr_id});
      if (marker.Max() == 1)
      {
        ConstantCoefficient Vcoef(bc.value);
        V.ProjectBdrCoefficient(Vcoef, marker);
        if (cfg->debug.debug) {
          std::cout << "\t\t" << name << " '(bdr_id = " << bc.bdr_id << ")'\n";
        }
      }
      else 
      {
        std::cerr << "\033[33m" << "WARNING: boundary " << name << "'(bdr_id= " << bc.bdr_id << ")' not present in mesh!\n" << "\033[0m\n";
      }
    }
}
