#include "load_mesh.h"
#include "mfem.hpp"
#include <iostream>
#include <cmath>

std::unique_ptr<mfem::ParMesh> CreateSimulationDomain(const std::string &path, MPI_Comm comm)
{
  // Load Mesh 
  // TODO Add flags for these options but I need to understand this  
  auto serial = std::make_unique<mfem::Mesh>(path.c_str(), 0, 0, false);
  if (serial->bdr_attributes.Size() == 0) { 
    std::cerr << "No boundary attributes!\n"; std::exit(1); 
  }
  // If I end up doing Adative Mesh Refinement
  //serial->EnsureNCMesh();
  // Parallelize if required
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, *serial);
  return pmesh;
}

void CheckAxisymmetricMesh(const mfem::Mesh &mesh, int radial_coord_index)
{
    // 1) Must be a 2D meridian mesh
    if (mesh.Dimension() != 2)
    {
        std::cerr << "[Axisym] ERROR: Axisymmetric mode requires a 2D (r,z) mesh.\n";
        std::exit(1);
    }

    // 2) Ensure r >= 0 across all vertices
    double rmin = 1e300, rmax = -1e300;
    for (int v = 0; v < mesh.GetNV(); ++v)
    {
        const double* X = mesh.GetVertex(v);
        const double r  = X[radial_coord_index];

        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
    }

    if (rmin < -1e-12 * std::max(1.0, rmax - rmin))
    {
        std::cerr << "[Axisym] ERROR: Mesh contains vertices with r < 0.\n"
                  << "         Axisymmetric meshes must lie entirely in r >= 0.\n"
                  << "         Minimum r detected: " << rmin << "\n";
        std::exit(1);
    }

    std::cout << "[Axisym] Mesh OK. r-range: [" << rmin << ", " << rmax << "]\n";
}
