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

// --------------------- Mesh Writers --------------------------
void SaveParallelMesh(const mfem::ParMesh &pmesh,
                      const std::filesystem::path &out_dir,
                      const std::string &prefix,
                      int precision = 16);

void SaveSerialMesh(const mfem::ParMesh &pmesh,
                    const std::filesystem::path &mesh_path,
                    int precision = 16,
                    const std::string &comments = "");

void SaveParallelVTU(mfem::ParMesh &pmesh,
                     const std::filesystem::path &out_prefix,
                     mfem::VTKFormat format = mfem::VTKFormat::ASCII,
                     bool high_order_output = true,
                     int compression_level = 0,
                     bool bdr_elements = false);

void SaveAMRMeshArtifacts(mfem::ParMesh &pmesh,
                          const std::filesystem::path &mesh_dir,
                          const AMRSettings &amr_io,
                          int precision = 16);

// --------------------- AMR --------------------------
bool ApplyAMRRefineDerefineStep(mfem::ParMesh &pmesh,
                                mfem::ParFiniteElementSpace &pfes,
                                mfem::ParGridFunction &V,
                                const Config &cfg);
