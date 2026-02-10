#include "mesh_tools.h"
#include "mfem.hpp"
#include <iostream>
#include <cmath>

std::unique_ptr<mfem::ParMesh> CreateSimulationDomain(const std::string &path, MPI_Comm comm)
{
  // Load Mesh 
  auto serial = std::make_unique<mfem::Mesh>(path.c_str(), 0, 1, false);
  if (serial->bdr_attributes.Size() == 0) { 
    std::cerr << "No boundary attributes!\n"; std::exit(1); 
  }
  // If I end up doing Adative Mesh Refinement
  //serial->EnsureNCMesh();
  // Parallelize if required
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, *serial);
  return pmesh;
}
void CheckAxisymmetricMesh(const mfem::ParMesh &mesh,
                           int radial_coord_index)
{
    MPI_Comm comm = mesh.GetComm();
    CheckAxisymmetricMesh(mesh,
                          radial_coord_index,
                          comm);
}
void CheckAxisymmetricMesh(const mfem::ParMesh &mesh,
                           int radial_coord_index,
                           MPI_Comm comm)
{
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // 1) Must be a 2D meridian mesh (local check is sufficient)
    if (mesh.Dimension() != 2)
    {
        if (rank == 0)
        {
            std::cerr << "[Axisym] ERROR: Axisymmetric mode requires a 2D (r,z) mesh.\n";
        }
        MPI_Abort(comm, 1);
    }

    // 2) Local extrema of r
    double rmin_loc =  1e300;
    double rmax_loc = -1e300;

    for (int v = 0; v < mesh.GetNV(); ++v)
    {
        const double *X = mesh.GetVertex(v);
        const double r  = X[radial_coord_index];
        rmin_loc = std::min(rmin_loc, r);
        rmax_loc = std::max(rmax_loc, r);
    }

    // 3) Global extrema
    double rmin = 0.0, rmax = 0.0;
    MPI_Allreduce(&rmin_loc, &rmin, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&rmax_loc, &rmax, 1, MPI_DOUBLE, MPI_MAX, comm);

    // 4) Global validation
    const double tol = 1e-12 * std::max(1.0, rmax - rmin);
    if (rmin < -tol)
    {
        if (rank == 0)
        {
            std::cerr << "[Axisym] ERROR: Mesh contains vertices with r < 0.\n"
                      << "         Axisymmetric meshes must lie entirely in r >= 0.\n"
                      << "         Global minimum r detected: " << rmin << "\n";
        }
        MPI_Abort(comm, 1);
    }

    // 5) Informational output (rank 0 only)
    if (rank == 0)
    {
        std::cout << "[Axisym] Mesh OK. Global r-range: ["
                  << rmin << ", " << rmax << "]\n";
    }
}