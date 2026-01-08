#include <memory>
#include <filesystem>
#include <string>

#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "solver_api.h"
#include "Config.h"
#include "cmdLineParser.h"
#include "parallelization.h"
#include "ComputeElectricField.h"
#include "config_modification.h"
#include "optimization.h"
#include "sweeps.h"

// -----------------------------------------------------------------------------
// main()
// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    cli::InputParser args(argc, argv);
    if (args.has("-h") || args.has("--help"))
    {
        cli::print_usage(argv[0]);
        return 0;
    }

    // ----------------------   Read options ---------------------------

    auto config_str_opt = args.get("-c");
    if (!config_str_opt) { config_str_opt = args.get("--config"); }
    if (!config_str_opt)
    {
        std::cerr << "Error: missing required argument -c/--config\n";
        cli::print_usage(argv[0]);
        return 1;
    }
    std::filesystem::path config_path = cli::to_absolute(*config_str_opt);
    if (!std::filesystem::exists(config_path))
    {
        std::cerr << "Error: config file not found: " << config_path << "\n";
        return 1;
    }

    std::cout << "[Config] " << config_path << "\n";

    // ----------------------   Load Config ---------------------------

    Config init_cfg = Config::Load(config_path.string());
    if (init_cfg.debug.debug && init_cfg.solver.axisymmetric) {std::cout<<"[DEBUG:MAIN] Running Axisymmetric Simulation" <<std::endl;}
    // Initialize parallel/MPI environment
    parallel::init_environment(init_cfg, argc, argv);

    const auto &sweeps = init_cfg.sweeps;

    // -----------  Single / Optimization / Sweep dispatch -------------

    // Do optimization
    if (init_cfg.optimize.enabled && (!init_cfg.optimize.metrics_only))
    {
        std::cout << "[MAIN] Executing an optimization" << std::endl;
        run_optimization(init_cfg);
        return 0;
    }
    // Do sweep
    else if (!sweeps.empty())
    {
        std::cout << "[MAIN] Executing a sweep simulation" << std::endl;
        std::size_t run_counter = 0;
        std::vector<std::pair<std::string, std::string>> active_params;
        std::vector<RunRecord> records;
        std::vector<Assignment> assignments;

        sweep_recursive_cfg(init_cfg,
                            sweeps,
                            /* idx        */ 0,
                            active_params,
                            assignments,
                            run_counter,
                            records);

        std::filesystem::path save_root(init_cfg.save_path);
        write_sweep_meta(init_cfg.geometry_id, save_root, records);
        return 0;
    }
    // Single Run     
    else
    {
        std::cout << "[MAIN] Executing a single simulation" << std::endl;
        std::size_t run_counter = 0;
        std::vector<std::pair<std::string, std::string>> active_params;
        std::vector<RunRecord> records;

        // Single run uses the same run_one + RunRecord machinery as a 1-point sweep
        std::string run_dir_name = run_one(init_cfg, active_params, run_counter);
        if (!run_dir_name.empty())
        {
            RunRecord rec;
            rec.run_dir_name = run_dir_name;
            rec.params       = active_params; // empty -> (none)
            records.push_back(std::move(rec));
        }
        ++run_counter; // used only for naming run_0001, harmless but consistent

        std::filesystem::path save_root(init_cfg.save_path);
        write_sweep_meta(init_cfg.geometry_id, save_root, records);
        return 0;
    }

    std::cerr << "Did Nothing" <<std::endl;
}
