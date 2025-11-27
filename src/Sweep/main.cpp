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


// -----------------------------------------------------------------------------
// One run: given Config and metadata, run simulation and save outputs + meta
// -----------------------------------------------------------------------------
static std::string run_one(const Config &cfg,
                           const std::vector<std::pair<std::string, std::string>> &active_params,
                           std::size_t run_index) // 0-based
{
    std::filesystem::path model_path = cfg.mesh.path;

    // Determine root save directory from Config
    std::filesystem::path save_root(cfg.save_path);
    std::error_code ec;

    std::filesystem::create_directories(save_root, ec);
    if (ec)
    {
        std::cerr << "Warning: could not create save_root directory "
                  << save_root << " : " << ec.message() << "\n";
    }

    // Per-run directory: run_0001, run_0002, ...
    std::size_t display_index = run_index + 1;  // make it 1-based for naming
    std::ostringstream run_name;
    run_name << "run_" << std::setw(4) << std::setfill('0') << display_index;

    std::filesystem::path run_dir = save_root / run_name.str();
    std::filesystem::create_directories(run_dir, ec);
    if (ec)
    {
        std::cerr << "Warning: could not create run directory "
                  << run_dir << " : " << ec.message() << "\n";
    }

    // Copy config and override solver output paths so they land in run_dir
    Config cfg_copy = cfg;

    auto cfg_ptr = std::make_shared<Config>(cfg_copy);
    if (!cfg.debug.dry_run)
    {
        SimulationResult result = run_simulation(cfg_ptr, model_path);
        if (!result.success)
        {
            std::cerr << "Simulation failed for " << run_name.str()
                      << ": " << result.error_message << "\n";
            return {};  // signal failure to caller
        }

        save_results(result, save_root);
    }

    if (cfg.debug.debug){
        std::cout << "[DEBUG] Boundary Condition values" << std::endl;
        for (const auto& [name, b] : cfg_copy.boundaries) {
            std::cout << name << ": " << b.value << '\n';
        }
    }
    return run_name.str();
}

// -----------------------------------------------------------------------------
// Recursive sweep over Config
// -----------------------------------------------------------------------------
static void sweep_recursive_cfg(const Config &base_cfg,
                                const std::vector<SweepEntry> &sweeps,
                                std::size_t idx,
                                std::vector<std::pair<std::string, std::string>> &active_params, // (label, value)
                                std::vector<Assignment> &assignments,                           // (path, value)
                                std::size_t &run_counter,
                                std::vector<RunRecord> &records)
{
    if (idx == sweeps.size())
    {
        // All sweep parameters fixed -> apply to a copy and run
        Config cfg = base_cfg;

        for (const auto &a : assignments)
        {
            set_cfg_value_from_string(cfg, a.path, a.value);
        }
        std::string run_dir_name = run_one(cfg, active_params, run_counter);

        if (!run_dir_name.empty())
        {
            RunRecord rec;
            rec.run_dir_name = run_dir_name;
            rec.params = active_params;
            records.push_back(std::move(rec));
        }

        ++run_counter;
        return;
    }

    const SweepEntry &sw = sweeps[idx];
    const std::string label = sw.name.empty() ? sw.path : sw.name;

    switch (sw.kind)
    {
        case SweepEntry::Kind::Discrete:
        {
            for (const auto &val : sw.values)
            {
                active_params.emplace_back(label, val);
                assignments.push_back({ sw.path, val });

                sweep_recursive_cfg(base_cfg, sweeps, idx + 1,
                                    active_params, assignments,
                                    run_counter, records);

                assignments.pop_back();
                active_params.pop_back();
            }
            break;
        }

        case SweepEntry::Kind::Range:
        {
            if (sw.steps <= 1)
            {
                std::string val_str = std::to_string(sw.start);

                active_params.emplace_back(label, val_str);
                assignments.push_back({ sw.path, val_str });

                sweep_recursive_cfg(base_cfg, sweeps, idx + 1,
                                    active_params, assignments,
                                    run_counter, records);

                assignments.pop_back();
                active_params.pop_back();
            }
            else
            {
                double step = (sw.end - sw.start) / double(sw.steps - 1);
                for (int i = 0; i < sw.steps; ++i)
                {
                    double v = sw.start + i * step;
                    std::string val_str = std::to_string(v);

                    active_params.emplace_back(label, val_str);
                    assignments.push_back({ sw.path, val_str });

                    sweep_recursive_cfg(base_cfg, sweeps, idx + 1,
                                        active_params, assignments,
                                        run_counter, records);

                    assignments.pop_back();
                    active_params.pop_back();
                }
            }
            break;
        }
    }
}

// -----------------------------------------------------------------------------
// Global meta.txt in save_root
// -----------------------------------------------------------------------------
static void write_sweep_meta(const std::string &geometry_id,
                              const std::filesystem::path &save_root,
                              const std::vector<RunRecord> &records)
{
    std::error_code ec;
    std::filesystem::create_directories(save_root, ec);
    if (ec)
    {
        std::cerr << "Warning: could not create save_root directory "
                  << save_root << " : " << ec.message() << "\n";
    }

    std::ofstream meta(save_root / "meta.txt");
    if (!meta)
    {
        std::cerr << "Warning: could not open meta.txt in " << save_root << "\n";
        return;
    }

    // Header: geometry_id as the "name" of this sweep
    meta << geometry_id << "\n\n";

    if (records.empty())
    {
        meta << "(no runs)\n";
        return;
    }

    for (const auto &rec : records)
    {
        meta << rec.run_dir_name << ":\n";
        if (rec.params.empty())
        {
            meta << "  (none)\n\n";
        }
        else
        {
            for (const auto &p : rec.params)
            {
                meta << "  " << p.first << " = " << p.second << "\n";
            }
            meta << "\n";
        }
    }
}

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

    // Initialize parallel/MPI environment
    parallel::init_environment(init_cfg, argc, argv);

    const auto &sweeps = init_cfg.sweeps;

    // -----------  Single / Optimization / Sweep dispatch -------------

    if (init_cfg.optimize.enabled)
    {
        std::cout << "[MAIN] Executing an optimization" << std::endl;
        run_optimization(init_cfg);
        return 0;
    }

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
