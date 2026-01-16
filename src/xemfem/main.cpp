#include <iostream>
#include <string>

#include "cmdLineParser.h"
#include "Config.h"
#include "optimization.h"
#include "sweeps.h"
#include "parallelization.h"


static int run_sim(Config init_cfg) {
    const auto &sweeps = init_cfg.sweeps;
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
        std::cout << "[MAIN] Executing a sweep" << std::endl;
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
        std::cout << "[MAIN] Executing a simulation" << std::endl;
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
        ++run_counter; // used only for naming run_0001

        std::filesystem::path save_root(init_cfg.save_path);
        return 0;
    }
    return 0;
}

static int run_metrics(Config init_cfg) {
    run_metrics_only(init_cfg);
    return 0;
}

static int run_plot(Config init_cfg) {
    std::cout << "[XEMFEM] Plotting to be implemented\n";
    // TODO: return xemfem_plot_main(argc, argv);
    return 0;
}

int main(int argc, char** argv)
{
    // -------------------- Check Subcommand ----------------------
    cli::InputParser pre_args(argc, argv);

    if (pre_args.has("-h") || pre_args.has("--help")) {
        cli::print_usage(argv[0]);
        return 0;
    }

    const std::string cmd = pre_args.subcommand().value_or("sim");

    // Strip subcommand 
    cli::InputParser::strip_subcommand(argc, argv);

    // ---------------------- Read options -----------------------------
    cli::InputParser args(argc, argv);

    auto config_str_opt = args.get("-c");
    if (!config_str_opt) { config_str_opt = args.get("--config"); }

    // Plot may not require a config
    std::filesystem::path config_path;
    if (!config_str_opt) {
        std::cout << "Using Default config path ../geometry/config.yaml \n";
        config_path = cli::to_absolute("../geometry/config.yaml");
    }
    else {
        config_path = cli::to_absolute(*config_str_opt);
    }

    if (!std::filesystem::exists(config_path)) {
        std::cout << "Error: config file not found: " << config_path << "\n";
        return 1;
    }
    std::cout << "[Config] " << config_path << "\n";

    // ---------------------- Load Config ------------------------------
    Config init_cfg;
    init_cfg = Config::Load(config_path.string());

    // MPI / OMP?
    parallel::init_environment(init_cfg, argc, argv);

    // ------------------------- Dispatch ------------------------------
    if (cmd == "sim") {
        return run_sim(init_cfg);
    }
    if (cmd == "metrics") {
        return run_metrics(init_cfg);
    }
    if (cmd == "plot") {
        // TODO Needs extra args?
        return run_plot(init_cfg);
    }

    std::cerr << "Error: unknown subcommand '" << cmd << "'\n";
    cli::print_usage(argv[0]);
    return 1;
}