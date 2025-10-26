#define MFEM_DEBUG
#include "mfem.hpp"
#include "geometry.h"
#include "constants.h"
#include "boundary_conditions.h"
#include "solver.h"
#include "ComputeElectricField.h"

using namespace mfem;

int main(int argc, char *argv[])
{
    // 1. Create the mesh
    Mesh* mesh = CreateSimulationDomain();
    if (debug)
    {
        std::cout << "bdr_attributes Max=" << mesh->bdr_attributes.Max()
            << " list: ";
        for (int i = 0; i < mesh->bdr_attributes.Size(); i++)
            std::cout << mesh->bdr_attributes[i] << " ";
        std::cout << "\nNBE=" << mesh->GetNBE() << "\n";

        int count1=0, count2=0, count3=0;
        for (int be = 0; be < mesh->GetNBE(); be++)
        {
            int a = mesh->GetBdrElement(be)->GetAttribute();
            if (a == 1) count1++;
            if (a == 2) count2++;
            if (a == 3) count3++;
        }
        std::cout << "boundary counts: attr1=" << count1
                << " attr2=" << count2
                << " attr3=" << count3 << "\n";

        std::cout << "Mesh dimensions: " 
                << mesh->GetNE() << " elements, "
                << mesh->GetNBE() << " boundary elements" << std::endl;
    }
    // 2. Create finite element collection and space
    H1_FECollection fec(fe_order, mesh->Dimension());      // Must survive as long as GridFunction
    FiniteElementSpace fespace(mesh, &fec);               // Must survive as long as GridFunction

    // 3. Get Dirichlet boundary attributes
    Array<int> dirichlet_arr = GetDirichletAttributes(mesh);

    // 4. Solve Poisson
    GridFunction V = SolvePoisson(fespace, dirichlet_arr);

    // Vector-valued FE space for E-field
    FiniteElementSpace *vec_fes = new FiniteElementSpace(mesh, &fec, mesh->Dimension());
    // Scalar FE space for |E|
    FiniteElementSpace *scalar_fes = new FiniteElementSpace(mesh, &fec, 1);

    // Allocate GridFunctions
    GridFunction E(vec_fes);     // vector field
    GridFunction Emag(scalar_fes); // magnitude

    // Compute vector E-field and its magnitude
    ComputeElectricField(V, E);      // fills E with -∇V
    ComputeFieldMagnitude(E, Emag);  // computes |E| at nodes

    // 5. Save Data
    mesh->Save("simulation_mesh.mesh");
    V.Save("solution_V.gf");
    Emag.Save("solution_Emag.gf");

    // TODO Update
    //double Enaive = V0 / plate_gap;  // V/m
    //std::cout << "Electric field between plates: " << Enaive*0.01 << " V/cm" << std::endl; 

}

