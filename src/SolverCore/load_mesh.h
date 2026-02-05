#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "mfem.hpp"

using namespace mfem;


std::unique_ptr<mfem::ParMesh> 
CreateSimulationDomain(const std::string &path,
                       MPI_Comm comm = MPI_COMM_WORLD
);

void CheckAxisymmetricMesh( const mfem::ParMesh &mesh, int radial_coord_index, MPI_Comm comm);
void CheckAxisymmetricMesh( const mfem::ParMesh &mesh, int radial_coord_index);

#endif
