
// ===================================== Nelder Mead ==================================
/* Uses Git YibaiMeng/nelder-mead header only library 
Requires a function that takes const std::vector<T>& returning T 
No built in bounds or constraints -> We need to handle these manually
*/
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <stdexcept>
#include <limits>
#include <cmath> 
#include "mfem.hpp"

#include "Config.h"
#include "solver_api.h"
#include "ComputeElectricField.h"
#include "config_modification.h"
#include "optimization.h"
#include "trace_fieldlines.h"

#include "nelder-mead.h" 


// ===================================== Main Loop ==================================== 

void run_optimization(const Config &init_cfg)
{
    const auto &opt = init_cfg.optimize;
    if (!opt.enabled) return;

    if (opt.variables.empty()) {
        throw std::runtime_error("optimize.enabled=true but no variables specified");
    }

    // Base config from which all evaluations start
    Config base_cfg = init_cfg;

    // Build initial point x0 from YAML
    std::vector<double> x0;
    x0.reserve(opt.variables.size());
    for (const auto &v : opt.variables) {
        double x_init = v.initial;
        if (x_init < v.lower || x_init > v.upper) {
            x_init = 0.5 * (v.lower + v.upper);
        }
        x0.push_back(x_init);
    }

    // Build objective function (metrics -> scalar)
    ObjectiveFn objective_fn = make_objective_function(opt);

    // Storage for all runs (for meta.txt)
    std::vector<OptRunRecord> records;
    records.reserve(opt.max_fun_evals > 0 ? opt.max_fun_evals : 100);

    // Counter to label runs as run_0001, run_0002, ...
    std::size_t eval_counter = 0;

    // Wrapper that Nelder–Mead will call for each x
    auto nm_func = [&](const std::vector<double> &x) -> double {
        double f = evaluate_one_optimization_point(
            base_cfg,
            opt,
            x,
            eval_counter,
            objective_fn,
            records
        );
        ++eval_counter;
        return f;
    };

    // Run Nelder–Mead
    std::vector<double> x_opt = nelder_mead::find_min(
        nm_func,
        x0,
        opt.adaptive,
        /* initial_simplex */ {},
        opt.tol_fun,
        opt.tol_x,
        static_cast<unsigned int>(opt.max_iters),
        static_cast<unsigned int>(opt.max_fun_evals)
    );

    // find the best recorded run
    double best_f = std::numeric_limits<double>::infinity();
    const OptRunRecord *best_rec = nullptr;
    for (const auto &rec : records) {
        if (rec.objective_value < best_f) {
            best_f = rec.objective_value;
            best_rec = &rec;
        }
    }

    if (best_rec) {
        std::cout << "[OPT] Best objective: " << best_f
                  << " in " << best_rec->run_dir_name << "\n";
    } else {
        std::cout << "[OPT] No successful runs recorded.\n";
    }

    // Write meta.txt summarizing all optimization runs
    std::filesystem::path save_root(init_cfg.save_path);
    write_optimization_meta(init_cfg.geometry_id, save_root, records);
}


RunAndMetricsResult run_simulation_and_save(const Config &cfg, std::size_t eval_index)
{
    RunAndMetricsResult out;

    std::filesystem::path model_path = cfg.mesh.path;
    std::filesystem::path save_root(cfg.save_path);

    std::error_code ec;
    std::filesystem::create_directories(save_root, ec);
    if (ec) {
        std::cerr << "Warning: could not create save_root directory "
                  << save_root << " : " << ec.message() << "\n";
    }

    // Make run name
    std::size_t display_index = eval_index + 1;
    std::ostringstream run_name;
    run_name << "run_" << std::setw(4) << std::setfill('0') << display_index;
    out.run_dir_name = run_name.str();

    std::filesystem::path run_dir = save_root / out.run_dir_name;
    std::filesystem::create_directories(run_dir, ec);
    if (ec) {
        std::cerr << "Warning: could not create run directory "
                  << run_dir << " : " << ec.message() << "\n";
    }

    // Copy config and override solver output paths
    Config cfg_copy = cfg;
    cfg_copy.save_path = run_dir.string();

    auto cfg_ptr = std::make_shared<Config>(cfg_copy);

    if (cfg.debug.dry_run) {
        // no actual simulation; out.success is false, but we still return run_dir_name
        out.success = false;
        return out;
    }

    SimulationResult result = run_simulation(cfg_ptr, model_path);
    if (!result.success) {
        std::cerr << "Simulation failed for " << out.run_dir_name
                  << ": " << result.error_message << "\n";
        out.success = false;
        return out;
    }

    // Save Results 
    save_results(result, run_dir);

    out.success = true;
    out.sim     = std::move(result);
    return out;
}


// ================================ Optimization Targets ==============================

double evaluate_one_optimization_point(const Config &base_cfg,
                                       const OptimizationSettings &opt,
                                       const std::vector<double> &x,
                                       std::size_t eval_index,
                                       const ObjectiveFn &objective_fn,
                                       std::vector<OptRunRecord> &records)
{
    // 1. Build config for this eval
    Config cfg = base_cfg;

    // apply optimization variables
    apply_opt_vars(cfg, opt, x);

    // apply field cage network (if enabled)
    apply_fieldcage_network(cfg);

    // 2. Run simulation and save outputs
    RunAndMetricsResult r = run_simulation_and_save(cfg, eval_index);

    if (!r.success) {
        // Penalize failures
        double penalty = 1e30;

        OptRunRecord rec;
        rec.run_dir_name    = r.run_dir_name; // may be empty if something went really wrong
        rec.objective_value = penalty;
        rec.vars            = make_var_list(opt, x);
        records.push_back(std::move(rec));

        return penalty;
    }

    // 3. Compute metrics
    OptimizationMetrics metrics = compute_metrics(cfg, r.sim);

    // 4. Evaluate objective
    double f = objective_fn(metrics);

    // 5. Record for meta.txt
    OptRunRecord rec;
    rec.run_dir_name    = r.run_dir_name;
    rec.objective_value = f;
    rec.vars            = make_var_list(opt, x);
    records.push_back(std::move(rec));

    return f;
}

ObjectiveFn make_objective_function(const OptimizationSettings &opt)
{
    const std::string &name = opt.objective;

    if (name == "CIV") {
        // Minimize CIV
        return [](const OptimizationMetrics &m) {
            return m.CIV;
        };
    }

    if (name == "weighted") {
        const double wC = opt.w_CIV;

        return [wC](const OptimizationMetrics &m) {
            return wC * m.CIV;
        };
    }

    // Default / unknown objective: error
    throw std::runtime_error("Unknown optimization objective: '" + name + "'");
}

OptimizationMetrics compute_metrics(const Config &cfg, const SimulationResult &result)
{
    OptimizationMetrics m;

    m.CIV = compute_civ(cfg, result);

    if (cfg.optimize.print_results)
    {
        std::cout << "[OPT] CIV: " << m.CIV << std::endl;
    }

    return m;
}

double compute_civ(const Config &cfg, const SimulationResult &result)
/*
Compute Charge Insensitive Volume on the mesh
Ie 1-(Volume where electrons reach liquid gas interface) / (Total Volume)
*/
{
    using namespace mfem;

    MFEM_VERIFY(result.mesh, "compute_civ: mesh is null");
    MFEM_VERIFY(result.E,    "compute_civ: E field is null");

    ParMesh &pmesh = *result.mesh;
    const int dim  = pmesh.Dimension();
    if (dim != 2) {
        throw std::runtime_error("compute_civ currently implemented for 2D axisymmetric (r,z) only.");
    }

    CivSeeds seeds = ExtractCivSeeds(cfg, result);

    // If no seeds throw error
    if (seeds.positions.empty()) {
        throw std::runtime_error("No seeds to compute CIV for");
    }
    
    // TODO Add config entry for this
    ElectronTraceParams trace_params;  
    std::vector<ElectronTraceResult> trace_results;

    TraceElectronFieldLines(result, cfg, seeds, trace_results);

    MFEM_VERIFY(trace_results.size() == seeds.positions.size(), "TraceElectronFieldLines: result size mismatch");

    double V_total = 0.0; // active TPC volume
    double V_civ   = 0.0; // Volume where electrons do not reach liquid gas interface

    for (std::size_t i = 0; i < seeds.positions.size(); ++i)
    {
        const double dV = seeds.volumes[i];
        const ElectronTraceResult &res = trace_results[i];

        V_total += dV;

        // Check if hit liquid gas or not 
        if (res.exit_code != ElectronExitCode::HitLiquidGas)
        {
            V_civ += dV;
        }
    }
    // Return volume fraction
    return V_civ / V_total;
}


// ================================ Helper function =================================
void apply_opt_vars(Config &cfg, const OptimizationSettings &opt, const std::vector<double> &x)
{
    for (std::size_t i = 0; i < opt.variables.size(); ++i) {
        const auto &var = opt.variables[i];
        double v = x[i];

        // enforce bounds
        if (v < var.lower) v = var.lower;
        if (v > var.upper) v = var.upper;

        // write into Config via existing path mechanism
        set_cfg_value_from_string(cfg, var.path, std::to_string(v));
    }
}

std::vector<std::pair<std::string, std::string>> make_var_list(const OptimizationSettings &opt, const std::vector<double> &x)
{
    std::vector<std::pair<std::string, std::string>> vars;
    vars.reserve(opt.variables.size());

    for (std::size_t i = 0; i < opt.variables.size(); ++i) {
        const auto &v = opt.variables[i];
        double val = x[i];
        if (val < v.lower) val = v.lower;
        if (val > v.upper) val = v.upper;
        vars.emplace_back(v.name.empty() ? v.path : v.name,
                          std::to_string(val));
    }
    return vars;
}

// ============================= Write Meta function ===================================
void write_optimization_meta(const std::string &geometry_id,
                             const std::filesystem::path &save_root,
                             const std::vector<OptRunRecord> &records)
{
    std::error_code ec;
    std::filesystem::create_directories(save_root, ec);
    if (ec) {
        std::cerr << "Warning: could not create save_root directory "
                  << save_root << " : " << ec.message() << "\n";
    }

    std::ofstream meta(save_root / "meta.txt");
    if (!meta) {
        std::cerr << "Warning: could not open meta.txt in " << save_root << "\n";
        return;
    }

    meta << geometry_id << "\n\n";

    if (records.empty()) {
        meta << "(no runs)\n";
        return;
    }

    for (const auto &rec : records) {
        std::filesystem::path run_dir = save_root / rec.run_dir_name;

        meta << rec.run_dir_name << ":\n";
        meta << "  objective: " << rec.objective_value << "\n";

        meta << "  params:\n";
        if (rec.vars.empty()) {
            meta << "    (none)\n";
        } else {
            for (const auto &p : rec.vars) {
                meta << "    " << p.first << " = " << p.second << "\n";
            }
        }

        meta << "  outputs:\n";
        meta << "    field_ex: " << (run_dir / "field_ex.gf").string() << "\n";
        meta << "    field_ey: " << (run_dir / "field_ey.gf").string() << "\n";
        meta << "    field_mag: " << (run_dir / "field_mag.gf").string() << "\n";
        meta << "    mesh: " << (run_dir / "simulation_mesh.msh.000000").string() << "\n";
        meta << "    V_solution: " << (run_dir / "solution_V.gf.000000").string() << "\n";
        meta << "\n";
    }
}
