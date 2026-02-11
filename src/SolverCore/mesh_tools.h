#pragma once

#include "mfem.hpp"
#include "Config.h"
#include <iostream>
#include <cmath>
#include <filesystem>

using namespace mfem;


std::unique_ptr<mfem::ParMesh> 
CreateSimulationDomain(const std::string &path,
                       MPI_Comm comm = MPI_COMM_WORLD
);

void CheckAxisymmetricMesh( const mfem::ParMesh &mesh, int radial_coord_index, MPI_Comm comm);
void CheckAxisymmetricMesh( const mfem::ParMesh &mesh, int radial_coord_index);


// --------------------- AMR --------------------------
bool ApplyAMRRefineDerefineStep(mfem::ParMesh &pmesh,
                                mfem::ParFiniteElementSpace &pfes,
                                mfem::ParGridFunction &V,
                                const Config &cfg);
