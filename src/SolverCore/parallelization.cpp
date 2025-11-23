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

    void init_environment(const Config &cfg,
                          int &argc, char **&argv)
    {
        if (g_initialized)
        {
            return; // already configured
        }
        if (cfg.debug.debug) std::cout << "[DEBUG] Configuring solver device" << std::endl;

        // 1) Initialize MPI via MFEM wrapper
        g_mpi_session = new mfem::MPI_Session(argc, argv);
        if (cfg.debug.debug) std::cout << "[DEBUG:MPI] MPI Session: " << g_mpi_session << std::endl;
        // TODO MPI 

        int world_size = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        
        // Threading 
        if (cfg.compute.threads.enabled)
        {
            int num_threads = 1;
            if (cfg.debug.debug) std::cout << "[DEBUG:OpenMP] Reports Max Thread count: " << omp_get_max_threads() << std::endl;
            if (cfg.compute.threads.num_auto || cfg.compute.threads.num <= 0)
            {
                num_threads = omp_get_max_threads();
            }
            else
            {
                num_threads = cfg.compute.threads.num;
            }

            // Set num threads via OpenMP
            omp_set_num_threads(num_threads);
            std::cout << "[OpenMP] Thread Count: " << num_threads << std::endl;
            // Sidenote HYPREs: hypre_SetNumThreads  calls  omp_set_num_threads


            if (cfg.compute.threads.affinity == "compact") {
                setenv("OMP_PROC_BIND", "close", 1);   // GCC/LLVM/OpenMP
            } else if (cfg.compute.threads.affinity == "scatter") {
                setenv("OMP_PROC_BIND", "spread", 1);
            } else {
                throw std::runtime_error("Unknown threads.affinity: " + cfg.compute.threads.affinity);
            }
            std::cout << "[OpenMP] Affinity: " << getenv("OMP_PROC_BIND") << std::endl;
            
        }
        
        /*
        Further device configuration

        MPI: 
        - On by default but only at device level 0 - i am not sure what happens if we use more MPI processes TODO See what I actually do with MPI
        - We can not mix devices between MPI threads, ie all must be GPU or CPU not a mix of them - allthough I do want to try to see if this can work
        - While n MPI produces similar idea as n OpenMP threads it is slower
        - Not really relevant at the moment 
        - On NUMA devices (2 CPU) using 2 MPI processes might be worthwhile? Should test 

        MFEM distinguishes between:
        - "device types": cpu / omp / cuda / hip (+ debug)
        - "implementation layers": native, RAJA, OCCA, CEED
        "cpu"
        - Pure host execution.
        - Always available, no extra build flags.
        - Sequential loops on each MPI rank.
        "omp"
        - Host CPU with OpenMP parallelism.
        - Enabled when MFEM_USE_OPENMP=ON.
        - Uses OpenMP threads on the CPU.
        "cuda"
        - NVIDIA GPU via native CUDA kernels.
        - Enabled when MFEM_USE_CUDA=ON.
        - Requires CUDA toolkit + NVIDIA driver.
        - Used either directly ("cuda") or indirectly (RAJA/OCCA/CEED CUDA backends).
        "hip"
        - AMD GPU via native HIP kernels.
        - Enabled when MFEM_USE_HIP=ON.
        - Requires ROCm/HIP stack on the system.
        "debug"
        - Special backend for debugging device/host memory consistency.
        - Very slow; used for correctness checks, not production.
        - Should not be combined with other device backends.

        RAJA backends : RAJA is a performance-portable loop abstraction.
        "raja-cpu" 
        "raja-omp" 
        "raja-cuda"
        "raja-hip" 

        OCCA backends : JIT-compiles kernels at runtime
        "occa-cpu" 
        "occa-omp" 
        "occa-cuda"
        "occa-hip" 

        CEED backends : CEED backends delegate operator evaluation to libCEED (CEED = Center for Efficient Exascale Discretization)
        Focuses on high-order finite elements for large problems to utilize GPU acceleration CPU vectorization and parallel scaling
        Idk this seems curious but i still dont quite understand 
        Also seems nontrivial to implement in the solver since it apparently applies y = A * x without building A which I dont understand. 
        "ceed-cpu" 
        "ceed-omp" 
        "ceed-cuda"
        "ceed-hip" 


        # Complete List
        CPU / host:
        "cpu", "omp",
        "raja-cpu", "raja-omp",
        "occa-cpu", "occa-omp",
        "ceed-cpu",
        "debug"

        NVIDIA GPUs:
        "cuda",
        "raja-cuda",
        "occa-cuda",
        "ceed-cuda" (libCEED + CUDA, also enables "cuda")

        AMD GPUs:
        "hip",
        "raja-hip",
        "occa-hip" (if built),
        "ceed-hip" (libCEED + HIP, also enables "hip")
        */
        g_device = new mfem::Device("omp");//device_type.c_str());

        g_initialized = true;
    }

} // namespace parallel
