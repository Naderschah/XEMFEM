#pragma once

#include <functional>
#include <filesystem>
#include <vector>
#include <string>
#include <utility>

// Holds the optimization metrics
struct OptimizationMetrics {
    double CIV = 0.0;          // Charge Insensitive Volume
};
// For writing the meta file 
struct OptRunRecord
{
    std::string run_dir_name;    
    double objective_value = 0.;
    std::vector<std::pair<std::string, std::string>> vars;
};
// Internal result object
struct RunAndMetricsResult {
    std::string  run_dir_name;
    bool         success = false;
    SimulationResult sim; 
};

enum class FieldLineDest {
    Gate,       // reaches liquid-gas interface in TPC aperture
    PTFE,       // ends on PTFE
    Cathode,
    Other       // anything else (lost, ambiguous, etc.)
};


using ObjectiveFn = std::function<double(const OptimizationMetrics&)>;

void run_optimization(const Config &init_cfg);

RunAndMetricsResult run_simulation_and_save(const Config &cfg, std::size_t eval_index);

double evaluate_one_optimization_point(const Config &base_cfg, const OptimizationSettings &opt, const std::vector<double> &x, std::size_t eval_index, const ObjectiveFn &objective_fn, std::vector<OptRunRecord> &records);

ObjectiveFn make_objective_function(const OptimizationSettings &opt);

OptimizationMetrics compute_metrics(const Config &cfg, const SimulationResult &result);

double compute_civ(const Config &cfg, const SimulationResult &result);

void apply_opt_vars(Config &cfg, const OptimizationSettings &opt, const std::vector<double> &x);

std::vector<std::pair<std::string, std::string>> make_var_list(const OptimizationSettings &opt, const std::vector<double> &x);

void write_optimization_meta(const std::string &geometry_id, const std::filesystem::path &save_root, const std::vector<OptRunRecord> &records);