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

    void init_environment(Config &cfg)
    {
        if (g_initialized)
        {
            return; // already configured
        }
        if (cfg.debug.debug) std::cout << "[DEBUG] Configuring solver device" << std::endl;

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int world_size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        if (cfg.compute.mpi.enabled && (world_size > 1))
        {
            if (rank == 0) { std::cout << "[MPI] number of processes = " << world_size << "\n"; }
            omp_set_num_threads(1);

        }
        
        // Threading 
        if (cfg.compute.threads.enabled)
        {
            // TODO THis might not set threads correctly for all dependencies due to the env variables for MPI, otherwise MPI runs with all threads
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

            g_device = new mfem::Device("omp");
            g_initialized = true;
            
        }
        // Setting single threaded
        if (not g_initialized) { 
            g_device = new mfem::Device("cpu");
            g_initialized = true;
        }
        
    }

    void init_mpi(int &argc, char **&argv)
    {
        g_mpi_session = new mfem::MPI_Session(argc, argv);
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
    }
    

} // namespace parallel
