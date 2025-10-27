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
    Array<int> ess = MakeBdrMarker(mesh, {anode_BC_index, cathode_BC_index}); // mark only if present
    return ess;
}

void ApplyDirichletValues(GridFunction &V, const Array<int> &/*unused*/)
{
    Mesh *mesh = V.FESpace()->GetMesh();

    // TODO Wrap them in if use this

    // Anode Voltage
    {
        Array<int> m = MakeBdrMarker(mesh, {anode_BC_index});
        if (m.Max() == 1) { ConstantCoefficient Vanode(VAnode);
                            V.ProjectBdrCoefficient(Vanode, m); }
    }
    // Gate Voltage
    {
        Array<int> m = MakeBdrMarker(mesh, {gate_BC_index});
        if (m.Max() == 1) { ConstantCoefficient Vcat(VGate);
                            V.ProjectBdrCoefficient(Vcat, m); }
    }
    // Cathode Voltage
    {
        Array<int> m = MakeBdrMarker(mesh, {cathode_BC_index});
        if (m.Max() == 1) { ConstantCoefficient Vcat(VCathode);
                            V.ProjectBdrCoefficient(Vcat, m); }
    }
}
