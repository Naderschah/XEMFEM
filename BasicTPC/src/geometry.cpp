#include "geometry.h"
#include <iostream>

Mesh* CreateSimulationDomain()
{
    std::cout << "[Geometry] Loading mesh from " << msh_file << std::endl;
    Mesh *mesh = new Mesh(msh_file);

    if (mesh->bdr_attributes.Size() == 0)
    {
        std::cerr << "Error: Mesh has no boundary attributes!" << std::endl;
        std::exit(1);
    }

    mesh->EnsureNCMesh();

    std::cout << "[Geometry] Loaded mesh with "
              << mesh->GetNE() << " elements and "
              << mesh->bdr_attributes.Size() << " boundary attributes."
              << std::endl;

    return mesh;
}