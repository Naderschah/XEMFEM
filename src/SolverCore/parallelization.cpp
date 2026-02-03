// parallel/runtime.cpp
#include "parallelization.h"

#include <stdexcept>
#include <string>

#include <omp.h>
#include "HYPRE.h"

namespace parallel
{
    // Guard against multiple initializations in one process.
    static bool g_initialized = false;

    static mfem::MPI_Session *g_mpi_session = nullptr;
    static mfem::Device      *g_device      = nullptr;

    void init_environment(Config &cfg,
                          int &argc, char **&argv)
    {
        if (g_initialized)
        {
            return; // already configured
        }
        if (cfg.debug.debug) std::cout << "[DEBUG] Configuring solver device" << std::endl;

        // 1) Initialize MPI via MFEM wrapper
        if (cfg.compute.mpi.enabled)
        {
            g_mpi_session = new mfem::MPI_Session(argc, argv);
            if (cfg.debug.debug) std::cout << "[DEBUG:MPI] MPI Session: " << g_mpi_session << std::endl;
            // TODO MPI 
            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            int world_size = 0;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);
            // Lazy non host output suppression 
            if (rank != 0)
            {
                std::cout.setstate(std::ios::failbit);
            }
            if (rank == 0) { std::cout << "[MPI] number of processes = " << world_size << "\n"; }
        }
        
        // Threading 
        if (false)
        {
            int num_threads = 1;
            if (cfg.debug.debug) std::cout << "[DEBUG:OpenMP] Reports Max Thread count: " << omp_get_max_threads() << std::endl;
            if (cfg.compute.threads.num_auto || cfg.compute.threads.num <= 0)
            {
                num_threads = omp_get_max_threads();
                cfg.compute.threads.num = num_threads;
            }
            else
            {
                num_threads = cfg.compute.threads.num;
            }

            // Set num threads via OpenMP
            omp_set_num_threads(num_threads);
            // Sidenote HYPREs: hypre_SetNumThreads  calls  omp_set_num_threads

            if (cfg.compute.threads.affinity == "compact") {
                setenv("OMP_PROC_BIND", "close", 1);   // GCC/LLVM/OpenMP
            } else if (cfg.compute.threads.affinity == "scatter") {
                setenv("OMP_PROC_BIND", "spread", 1);
            } else {
                throw std::runtime_error("Unknown threads.affinity: " + cfg.compute.threads.affinity);
            }
            std::cout << "[OpenMP] Thread Count: " << num_threads 
                      << "; Affinity: " << getenv("OMP_PROC_BIND") << std::endl;
            
        }

        if (cfg.compute.threads.enabled)
        {
            cfg.compute.threads.enabled = false;
            std::cout << "[WARNING] Disabled OMP Functionality is to be removed" << std::endl;
        }
        
        g_device = new mfem::Device("cpu");

        g_initialized = true;
    }
    

} // namespace parallel
