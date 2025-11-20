#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "mfem.hpp"

using namespace mfem;


std::unique_ptr<mfem::Mesh>
CreateSimulationDomain(const std::string &path,
                       bool use_distributed,
#ifdef MFEM_USE_MPI
                       MPI_Comm comm = MPI_COMM_WORLD
#else
                       int comm = 0   // ignored when MPI is off
#endif
);

void CheckAxisymmetricMesh( const mfem::Mesh &mesh,
                            int radial_coord_index,   // 0 = x, 1 = y
                            int axis_bdr_attr);

#endif
