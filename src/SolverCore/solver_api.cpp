
#include "mfem.hpp"
#include "load_mesh.h"
#include "boundary_conditions.h"
#include "solver.h"
#include "ComputeElectricField.h"
#include "config/Config.h"
#include "cmdLineParser.h"
#include "solver_api.h"


#ifdef MFEM_USE_OPENMP
  #include <omp.h>
#endif


using namespace mfem;

SimulationResult run_simulation()
{
// -------------------------------- Parallelization ----------------------------------------------
  bool use_distributed = cfg->compute.mpi.enabled;
  bool use_threads = cfg->compute.threads.enabled;  
  bool use_device  = (cfg->compute.device.type != "none" && cfg->compute.device.type != "cpu");
  // MPI Sesion 
  #ifdef MFEM_USE_MPI
    mfem::MPI_Session mpi(argc, argv);
    MPI_Comm comm = MPI_COMM_WORLD;
  #else
    MPI_Comm comm = 0;
  #endif
  // Open MP (multithreadding)
  #ifdef MFEM_USE_OPENMP
  if (use_threads) {
    int threads = cfg->compute.threads.num;
    if (threads <= 0) {
      omp_set_dynamic(1);
      threads = std::thread::hardware_concurrency();
      omp_set_num_threads(threads);
      std::cout << "Using automatic thread assignment" << std::endl;
      #pragma omp parallel
      {
        #pragma omp master
        std::cout << "[Threads] running with " << omp_get_num_threads()
                  << " OpenMP thread cap\n";
      }
    }
    else {
      omp_set_dynamic(0);
      omp_set_num_threads(threads);
      #pragma omp parallel
      {
        #pragma omp master
        std::cout << "[Threads] running with " << omp_get_num_threads()
                  << " OpenMP threads\n";
      }
    }
    // Map affinity: compact | scatter | none
    if      (cfg->compute.threads.affinity == "compact") 
    { setenv("OMP_PROC_BIND", "close", 1); setenv("OMP_PLACES", "cores", 1); }
    else if (cfg->compute.threads.affinity == "scatter") 
    { setenv("OMP_PROC_BIND", "spread", 1); setenv("OMP_PLACES", "cores", 1); }
    else                                                 
    { setenv("OMP_PROC_BIND", "false", 1); }
    }
    else {
      // If OpenMP is compiled in but disabled in config, force 1 thread:
      omp_set_num_threads(1);
      setenv("OMP_PROC_BIND", "false", 1);
    }
  #endif
  // ------------------------------ End Multiprocessing ------------------------------------s

 
  // 1. Create the mesh
  auto mesh = CreateSimulationDomain(model_path, use_distributed, comm); 
  if (cfg->solver.axisymmetric) {
    if (mesh->Dimension() != 2) { std::cerr << "Axisymmetric Simulation Geometry 3D" << std::endl;  }
    // Check r axis starts at null
    CheckAxisymmetricMesh(*mesh,
                          0,
                          cfg->solver.axisymmetric_r0_bd_attribute);
  }

  // 2. Create finite element collection and space
  H1_FECollection fec(cfg->solver.order, mesh->Dimension());
  FiniteElementSpace fespace(mesh.get(), &fec);

  // 3. Get Dirichlet boundary attributes
  Array<int> dirichlet_arr = GetDirichletAttributes(mesh.get(), cfg);
  // FIXME: For Neumann and Robin and axisymmetric we still need to supply boundary markers


  // 4. Solve Poisson
  auto V = SolvePoisson(fespace, dirichlet_arr, cfg);

  // Pick classes based on parallelization
  std::unique_ptr<mfem::FiniteElementSpace> vec_fes, scalar_fes;
  #ifdef MFEM_USE_MPI
  if (use_distributed) {
    auto *pmesh = static_cast<mfem::ParMesh*>(mesh.get());
    vec_fes    = std::make_unique<mfem::ParFiniteElementSpace>(pmesh, &fec, pmesh->Dimension());
    scalar_fes = std::make_unique<mfem::ParFiniteElementSpace>(pmesh, &fec, 1);
  } else
  #endif
  {
    vec_fes    = std::make_unique<mfem::FiniteElementSpace>(mesh.get(), &fec, mesh->Dimension());
    scalar_fes = std::make_unique<mfem::FiniteElementSpace>(mesh.get(), &fec, 1);
  }

  // 1) Initialize the postprocessor
  InitFieldPostprocessor(fespace, /*smooth_output=*/false);
  // 2) Allocate output fields on the right spaces
  auto E    = CreateE();     // vector L2 (vdim=dim)
  auto Emag = CreateEmag();  // scalar L2

  // 3) Compute E and |E|
  ComputeElectricField(*V, *E, /*scale=*/-1.0); // or -0.01 for V/cm
  ComputeFieldMagnitude(*E, *Emag);

  // Make return struct 
  SimulationResult result;
  result.mesh = std::move(mesh);
  result.V    = std::move(V);
  result.E    = std::move(E);
  result.Emag = std::move(Emag);

  result.success = true;

  return result;
}