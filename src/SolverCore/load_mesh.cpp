#include "load_mesh.h"
#include "mfem.hpp"
#include <iostream>
#include <cmath>

std::unique_ptr<mfem::Mesh>
CreateSimulationDomain(const std::string &path, bool use_distributed,
#ifdef MFEM_USE_MPI
                       MPI_Comm comm
#else
                       int /*comm*/
#endif
) {
  // Load Mesh  -  generate edges false, refine true, fix_orientation false TODO Add flag for the latter once i understand waht exactly it does
  auto serial = std::make_unique<mfem::Mesh>(path.c_str(), 0, 1, false);
  if (serial->bdr_attributes.Size() == 0) { 
    std::cerr << "No boundary attributes!\n"; std::exit(1); 
  }
  serial->EnsureNCMesh();
  // Parallelize if required
  #ifdef MFEM_USE_MPI
    if (use_distributed) {
      return std::make_unique<mfem::ParMesh>(comm, *serial); // upcasts to Mesh*
    }
  #endif
  return serial;
}

void CheckAxisymmetricMesh(const mfem::Mesh &mesh,
                                  int radial_coord_index,   // 0 = x, 1 = y
                                  int axis_bdr_attr)        // 0 if not checking boundary tag
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

    // 3) Optional: sanity check that boundary tag for axis is located at r = 0
    if (axis_bdr_attr > 0)
    {
        bool found = false;
        for (int be = 0; be < mesh.GetNBE(); ++be)
        {
            if (mesh.GetBdrAttribute(be) == axis_bdr_attr)
            {
                found = true;
                mfem::Array<int> v;
                mesh.GetBdrElementVertices(be, v);
                for (int k = 0; k < v.Size(); ++k)
                {
                    const double* X = mesh.GetVertex(v[k]);
                    if (std::fabs(X[radial_coord_index]) > 1e-12)
                    {
                        std::cerr << "[Axisym] WARNING: Boundary attribute "
                                  << axis_bdr_attr << " is not located exactly at r=0.\n";
                        break;
                    }
                }
            }
        }
        if (!found)
        {
            std::cerr << "[Axisym] WARNING: No boundary elements found with axis_bdr_attr="
                      << axis_bdr_attr << ".\n";
        }
    }

    std::cout << "[Axisym] Mesh OK. r-range: [" << rmin << ", " << rmax << "]\n";
}
