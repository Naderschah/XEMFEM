#include "boundary_conditions.h"
#include "constants.h"


// Build a safe boundary marker for the given attribute IDs.
static Array<int> MakeBdrMarker(const Mesh *mesh,
                                std::initializer_list<int> attrs)
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
Array<int> GetDirichletAttributes(Mesh *mesh)
{
    Array<int> ess = MakeBdrMarker(mesh, {1, 2}); // mark only if present
    return ess;
}

void ApplyDirichletValues(GridFunction &V, const Array<int> &/*unused*/)
{
    Mesh *mesh = V.FESpace()->GetMesh();

    // bottom plate (1) -> 0 V
    {
        Array<int> m = MakeBdrMarker(mesh, {1});
        if (m.Max() == 1) { ConstantCoefficient Vbot(0.0);
                            V.ProjectBdrCoefficient(Vbot, m); }
    }
    // top plate (2) -> V0
    {
        Array<int> m = MakeBdrMarker(mesh, {2});
        if (m.Max() == 1) { ConstantCoefficient Vtop(V0);
                            V.ProjectBdrCoefficient(Vtop, m); }
    }
}
