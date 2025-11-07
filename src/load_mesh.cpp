#include "load_mesh.h"
#include <iostream>

std::unique_ptr<mfem::Mesh>
CreateSimulationDomain(const std::string &path, bool use_distributed,
#ifdef MFEM_USE_MPI
                       MPI_Comm comm
#else
                       int /*comm*/
#endif
) {
  // Load Mesh 
  auto serial = std::make_unique<mfem::Mesh>(path.c_str());
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

