#include <iostream>
#include <string>
#include <omp.h>

#include "cmdLineInteraction.h"
#include "Config.h"
#include "optimization.h"
#include "sweeps.h"
#include "parallelization.h"
#ifdef HAVE_VTK
    #include "plotting_api.h"
#endif
#include "interpolator.h"
#include "path_handler.h"

std::string ReplaceMeshPathInConfig(const std::string& config_str,
                                    const std::string& new_mesh_path)
{
    YAML::Node config = YAML::Load(config_str);
    config["mesh"]["path"] = new_mesh_path;
    YAML::Emitter out;
    out << config;
    return std::string(out.c_str());
}

static int run_sim(Config init_cfg, std::string config_str) {
    const auto &sweeps = init_cfg.sweeps;

    if (init_cfg.mesh.amr.enable)
    {
        std::cout << "[MAIN] Producing AMR Mesh" << std::endl;
        std::string mesh_path = PrecomputeAMRMesh(init_cfg);
        init_cfg.mesh.path = mesh_path;
        config_str = ReplaceMeshPathInConfig(config_str, mesh_path);
    }

    // Do optimization
    if (init_cfg.optimize.enabled && (!init_cfg.optimize.metrics_only))
    {
        std::cout << "[MAIN] Executing an optimization" << std::endl;
        run_optimization(init_cfg, config_str);
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
                            config_str, 
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
        auto root = YAML::Load(config_str);
        std::string run_dir_name = run_one(init_cfg, root, active_params, run_counter, false);
        if (!run_dir_name.empty())  
        {
            RunRecord rec;
            rec.run_dir_name = run_dir_name;
            rec.params       = active_params; // empty -> (none)
            records.push_back(std::move(rec));
        }
        ++run_counter; // used only for naming run_0001

        std::filesystem::path save_root(init_cfg.save_path);
        write_sweep_meta(init_cfg.geometry_id, save_root, records);
        return 0;
    }
    return 0;
}

static int run_metrics(Config init_cfg) {
    run_metrics_only(init_cfg);
    return 0;
}

static int run_plot(Config init_cfg) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //if (rank == 0){ make_plot_api(init_cfg); }
    std::cout << "Plotting Currently Broken and removed from runtime in favor of compiletime" << std::endl;
    return 0;
}

static int run_interpolate(Config init_cfg) {
    do_interpolate(init_cfg);
    return 0;
}

int main(int argc, char** argv)
{    
    // --------------------- MPI needs to be innited early ----------------------------
    parallel::init_mpi(argc, argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_size);
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

    std::string config_path_str;
    int config_ok = 1;

    if (rank == 0) {
        auto config_str_opt = args.get("-c");
        if (!config_str_opt) { config_str_opt = args.get("--config"); }

        std::filesystem::path config_path;
        if (!config_str_opt) {
            std::cout << "Using Default config path ../geometry/config.yaml \n";
            config_path = cli::to_absolute("../geometry/config.yaml");
        } else {
            config_path = cli::to_absolute(*config_str_opt);
        }

        config_path_str = config_path.string();

        if (!std::filesystem::exists(config_path)) {
            std::cout << "Error: config file not found: " << config_path << "\n";
            config_ok = 0;
        } else {
            std::cout << "[Config] " << config_path << "\n";
        }
    }

    // Broadcast config_ok
    MPI_Bcast(&config_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!config_ok) {
        return 1;
    }

    // Broadcast config path string
    std::uint64_t path_n = 0;
    if (rank == 0) path_n = static_cast<std::uint64_t>(config_path_str.size());
    MPI_Bcast(&path_n, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    config_path_str.resize(static_cast<std::size_t>(path_n));
    if (path_n > 0) {
        MPI_Bcast(config_path_str.data(), static_cast<int>(path_n), MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    // ---------------------- Load Config ------------------------------
    std::string config_str = ReadConfigString(config_path_str, MPI_COMM_WORLD);
    Config init_cfg;
    init_cfg = Config::LoadFromString(config_str, cmd); // MPI handled internally
    // MPI Set Up
    parallel::init_environment(init_cfg);

    // Mark the number of ranks used 
    {
        YAML::Node root_tmp = YAML::Load(config_str);
        root_tmp["mpi"]["ranks"] = world_size;
        YAML::Emitter out;
        out << root_tmp;
        config_str = std::string(out.c_str());
    }
    // ------------------------- Dispatch ------------------------------
    if (cmd == "sim") {
        ensure_directory(init_cfg);
        run_sim(init_cfg, config_str);
        std::cout << "Done" <<std::endl;
        return 0;
    }
    if (cmd == "metrics") {
        return run_metrics(init_cfg);
        std::cout << "Done" <<std::endl;
        return 0;
    }
    if (cmd == "plot") {
        // TODO Needs extra args?
        if (rank == 0)
            run_plot(init_cfg);
            std::cout << "Done" <<std::endl;
        return 0;
    }
    if (cmd == "interpolate") {
        // TODO Needs extra args?
        run_interpolate(init_cfg);
        std::cout << "Done" <<std::endl;
        return 0;
    }

    std::cerr << "Error: unknown subcommand '" << cmd << "'\n";
    cli::print_usage(argv[0]);
    return 1;
}