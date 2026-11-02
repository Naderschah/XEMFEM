#include "optimization.h"

// ===================================== Main Loop ==================================== 

void run_optimization(const Config &init_cfg, const ConfigDocument &base_doc)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, world_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    const auto &opt = init_cfg.optimize;

    // Single exit flag; function returns once at the end.
    bool do_run = opt.enabled;

    // Command protocol:
    const int CMD_EVAL = 1;
    const int CMD_STOP = 0;

    // Variables used across both roles (root/worker)
    ObjectiveFn objective_fn;                 // constructed only if do_run
    std::vector<OptRunRecord> records;
    OptimizationLogger *logger_ptr = nullptr;
    std::unique_ptr<OptimizationLogger> logger_owner;
    std::filesystem::path save_root(init_cfg.save_path);

    std::vector<double> x0;
    int nvars = 0;

    // -------- Setup (only if enabled) --------
    if (do_run)
    {
        if (opt.variables.empty()) {
            throw std::runtime_error("optimize.enabled=true but no variables specified");
        }

        objective_fn = make_objective_function(opt);

        if (rank == 0) {
            records.reserve(opt.max_fun_evals > 0 ? opt.max_fun_evals : 100);
            logger_owner = std::make_unique<OptimizationLogger>(init_cfg.geometry_id, save_root);
            logger_owner->initialize();
            logger_ptr = logger_owner.get();
        } else {
            records.clear();
            records.shrink_to_fit();
        }

        nvars = static_cast<int>(opt.variables.size());

        // Build initial point x0 (rank 0), broadcast to all ranks
        if (rank == 0) {
            x0.reserve(opt.variables.size());
            for (const auto &v : opt.variables) {
                double x_init = v.initial;
                if (x_init < v.lower || x_init > v.upper) {
                    x_init = 0.5 * (v.lower + v.upper);
                }
                x0.push_back(x_init);
            }
        } else {
            x0.resize(opt.variables.size());
        }

        MPI_Bcast(x0.data(), nvars, MPI_DOUBLE, 0, comm);
    }

    // -------- Worker role --------
    if (do_run && rank != 0)
    {
        bool running = true;

        while (running)
        {
            int cmd = CMD_STOP;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);

            if (cmd == CMD_STOP) {
                running = false;
                continue;
            }

            std::uint64_t eval_index_u64 = 0;
            MPI_Bcast(&eval_index_u64, 1, MPI_UINT64_T, 0, comm);

            std::vector<double> x(nvars, 0.0);
            MPI_Bcast(x.data(), nvars, MPI_DOUBLE, 0, comm);

            (void)evaluate_one_optimization_point(
                base_doc,
                opt,
                x,
                static_cast<std::size_t>(eval_index_u64),
                objective_fn,
                records,
                /*logger*/ nullptr
            );
        }


        // Optional receive final x_opt
        int xopt_n = 0;
        MPI_Bcast(&xopt_n, 1, MPI_INT, 0, comm);
        std::vector<double> x_opt(xopt_n, 0.0);
        if (xopt_n > 0) {
            MPI_Bcast(x_opt.data(), xopt_n, MPI_DOUBLE, 0, comm);
        }
    }

    // -------- Root role --------
    if (do_run && rank == 0)
    {
        std::size_t eval_counter = 0;

        auto nm_func = [&](const std::vector<double> &x) -> double
        {
            int cmd = CMD_EVAL;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);

            std::uint64_t eval_index_u64 = static_cast<std::uint64_t>(eval_counter);
            MPI_Bcast(&eval_index_u64, 1, MPI_UINT64_T, 0, comm);

            MFEM_VERIFY(static_cast<int>(x.size()) == nvars, "nm_func: x size mismatch");
            MPI_Bcast(const_cast<double*>(x.data()), nvars, MPI_DOUBLE, 0, comm);

            double f = evaluate_one_optimization_point(
                base_doc,
                opt,
                x,
                eval_counter,
                objective_fn,
                records,
                logger_ptr
            );

            ++eval_counter;
            return f;
        };

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

        // Tell workers to stop
        {
            int cmd = CMD_STOP;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);
        }

        // Broadcast final x_opt
        {
            int xopt_n = static_cast<int>(x_opt.size());
            MPI_Bcast(&xopt_n, 1, MPI_INT, 0, comm);
            if (xopt_n > 0) {
                MPI_Bcast(x_opt.data(), xopt_n, MPI_DOUBLE, 0, comm);
            }
        }

        // Best run summary (rank 0 only)
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


        if (rank == 0){ write_optimization_meta(init_cfg.geometry_id, save_root, records); }
    }

    // -------- Single return --------
    return;
}



RunAndMetricsResult run_simulation_and_save(YAML::Node yaml_config, std::size_t eval_index)
{
    RunAndMetricsResult out;

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, world_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    if (!yaml_config["save_path"])
        throw std::runtime_error("run_simulation_and_save: YAML missing save_path");
    const std::filesystem::path run_dir(yaml_config["save_path"].as<std::string>());
    out.run_dir_name = yaml_config["save_path"].as<std::string>();
    // -------------------------------------------------------------------------
    // 3) Build Config from YAML (no struct-level mutation for these fields)
    // -------------------------------------------------------------------------
    Config cfg_copy = Config::LoadFromNode(yaml_config);

    // Derive paths from cfg_copy (now reflects YAML mutations)
    std::filesystem::path model_path(cfg_copy.mesh.path);

    // -------------------------------------------------------------------------
    // 4) Run simulation
    // -------------------------------------------------------------------------
    auto cfg_ptr = std::make_shared<Config>(cfg_copy);

    if (cfg_copy.debug.dry_run) {
        out.success = false;
        return out;
    }

    SimulationResult result = run_simulation(cfg_ptr, model_path, true);
    if (!result.success) {
        std::cerr << "Simulation failed for " << out.run_dir_name
                  << ": " << result.error_message << "\n";
        out.success = false;
        return out;
    }

    // -------------------------------------------------------------------------
    // 5) Save results (includes YAML dump; save_results is already MPI-aware)
    // -------------------------------------------------------------------------
    save_results(result, run_dir, yaml_config);

    out.success = true;
    out.sim     = std::move(result);
    return out;
}


// ================================ Optimization Targets ==============================

double evaluate_one_optimization_point(const ConfigDocument &base_doc,
                                       const OptimizationSettings &opt,
                                       const std::vector<double> &x,
                                       std::size_t eval_index,
                                       const ObjectiveFn &objective_fn,
                                       std::vector<OptRunRecord> &records,
                                       OptimizationLogger *logger)
{
    std::cout << "[OPTIMIZATION] Running iteration " << eval_index << std::endl;
    // 1. Build config for this eval
    YAML::Node root = base_doc.Root();

    // apply optimization variables
    apply_opt_vars(root, opt, x);

    // apply field cage network (if enabled)
    apply_fieldcage_network(root);

    // 2. Run simulation and save outputs
    root["save_path"] = make_run_folder(root, eval_index).string();
    RunAndMetricsResult r = run_simulation_and_save(root, eval_index);

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
    Config tmp_cfg = Config::LoadFromNode(root);
    OptimizationMetrics metrics = compute_metrics(tmp_cfg, r.sim);

    // 4. Evaluate objective
    double f = objective_fn(metrics);

    // 5. Record for meta.txt
    OptRunRecord rec;
    rec.run_dir_name    = r.run_dir_name;
    rec.objective_value = f;
    rec.vars            = make_var_list(opt, x);
    records.push_back(std::move(rec));

    if (logger) {
      logger->append_record(records.back());
    }

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
    if (name == "FieldSpread") {
      return [](const OptimizationMetrics &m) {
            return m.FieldSpread;
        };
    }


    if (name == "self_weighting") {
        return [CIV0 = 0.0, FS0 = 0.0]
              (const OptimizationMetrics &m) mutable
        {
            if (CIV0 == 0.0) CIV0 = m.CIV;
            if (FS0  == 0.0) FS0  = m.FieldSpread;

            const double CIV_hat = m.CIV / CIV0;
            const double FS_hat  = m.FieldSpread / FS0;

            return CIV_hat + FS_hat;
        };
    }


    // Default / unknown objective: error
    throw std::runtime_error("Unknown optimization objective: '" + name + "'");
}

double compute_civ(const Config &cfg,
                   const SimulationResult &result,
                   ElectronFieldLineTracer::Statistics *trace_stats)
/*
Compute Charge Insensitive Volume on the mesh
Ie 1-(Volume where electrons reach liquid gas interface) / (Total Volume)
*/
{
    std::chrono::steady_clock::time_point t_start;
    if (cfg.debug.timing)
    {
        t_start = std::chrono::steady_clock::now();
    }

    const std::string &method = cfg.civ_params.method; // "InformedSweep" or "RandomSample"
    double civ = 1.0;

    if (method == "RandomSample")
    {
        civ = ComputeCIV_RandomSample(cfg, result, trace_stats);
    }
    else if (method == "Grid")
    {
        civ = ComputeCIV_FixedGrid(cfg, result, trace_stats);
    }
    else if (method == "AdaptiveGrid")
    {
        civ = ComputeCIV_AdaptiveGrid(cfg, result, trace_stats);
    }
    else if (method == "RowSweep")
    {
        civ = ComputeCIV_RowSweep(cfg, result, trace_stats);
    }
    else
    {
        std::cerr << "[WARN:OPTIMIZATION] Unknown CIV method '" << method << "',\n";
        civ = 1.0;
    }

    if (cfg.debug.timing)
    {
        auto t_end = std::chrono::steady_clock::now();
        auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
        std::cout << "[Timing]: end CIV (" << dt << " ms)" << std::endl;
    }

    return civ;
}



OptimizationMetrics compute_metrics(const Config &cfg, const SimulationResult &result)
{
    OptimizationMetrics m;
    if ((cfg.optimize.objective == "CIV") || (cfg.optimize.objective == "self_weighting"))
    {
        ElectronFieldLineTracer::Statistics trace_stats;
        m.CIV = compute_civ(cfg, result, &trace_stats);
        m.traced_electrons = trace_stats.electron_count;
        m.tracing_seconds = trace_stats.elapsed_seconds;
        if (cfg.optimize.print_results)
        {
            std::cout << "[OPT] CIV: " << m.CIV << std::endl;
        }
    }
    

    if ((cfg.optimize.objective == "FieldSpread") || (cfg.optimize.objective == "self_weighting"))
    {
        m.FieldSpread = computeFieldSpreadMetric(
            cfg, result, &m.field_homogeneity);
        if (cfg.optimize.print_results)
        {
            std::cout << "[OPT] FieldSpread: " << m.FieldSpread << std::endl;
        }
    }
    
    return m;
}

// ================================ Helper function =================================
void apply_opt_vars(YAML::Node &root,
                    const OptimizationSettings &opt,
                    const std::vector<double> &x)
{
    MFEM_VERIFY(opt.variables.size() == x.size(),
                "apply_opt_vars: size mismatch between variables and x");

    for (std::size_t i = 0; i < opt.variables.size(); ++i)
    {
        const auto &var = opt.variables[i];
        double v = x[i];

        if (v < var.lower) v = var.lower;
        if (v > var.upper) v = var.upper;

        if (!root) std::cout << "<root is undefined/null>\n";
        else if (!root.IsMap()) std::cout << "<root is not a map>\n";

        // Split path
        std::stringstream ss(var.path);
        std::string key;
        std::vector<std::string> keys;
        while (std::getline(ss, key, '.')) keys.push_back(key);

        MFEM_VERIFY(!keys.empty(),
                    "apply_opt_vars: empty path for variable " + var.name);

        // Rebind without operator= using optional.emplace()
        std::optional<YAML::Node> cur;
        cur.emplace(root);  // copy-construct handle

        for (std::size_t k = 0; k + 1 < keys.size(); ++k)
        {
            // Existence check WITHOUT creating nodes: use const operator[]
            const YAML::Node cur_const = *cur;
            YAML::Node child = cur_const[keys[k]];
            if (!child)
            {
                std::ostringstream msg;
                msg << "apply_opt_vars: invalid path '" << var.path
                    << "' (missing key '" << keys[k] << "')\n"
                    << "Available keys at this level:\n";

                if (cur_const.IsMap()) {
                    for (auto it = cur_const.begin(); it != cur_const.end(); ++it) {
                        try { msg << it->first.as<std::string>() << ", "; }
                        catch (...) { msg << "<non-string-key>, "; }
                    }
                } else {
                    msg << "<node is not a map>";
                }

                MFEM_VERIFY(false, msg.str());
            }

            // Now rebind to mutable child WITHOUT operator=
            cur.emplace((*cur)[keys[k]]);
        }

        // Set final value
        (*cur)[keys.back()] = v;
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

// ============================= Write Meta functions ===================================
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
    }

    meta << geometry_id << "\n\n";

    if (records.empty()) {
        meta << "(no runs)\n";
        meta.flush();
    }
    else
    {
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
    }}

    meta.flush();
}


void write_single_opt_record_block(std::ostream &meta,
                                   const std::filesystem::path &save_root,
                                   const OptRunRecord &rec)
{
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
    meta << "    field_ex: "   << (run_dir / "field_ex.gf").string() << "\n";
    meta << "    field_ey: "   << (run_dir / "field_ey.gf").string() << "\n";
    meta << "    field_mag: "  << (run_dir / "field_mag.gf").string() << "\n";
    meta << "    mesh: "       << (run_dir / "simulation_mesh.msh.000000").string() << "\n";
    meta << "    V_solution: " << (run_dir / "solution_V.gf.000000").string() << "\n";
    meta << "\n";
}

// ------------------------------
// OptimizationLogger
// ------------------------------

OptimizationLogger::OptimizationLogger(std::string geometry_id,
                                       std::filesystem::path save_root)
    : geometry_id_(std::move(geometry_id))
    , save_root_(std::move(save_root))
    , initialized_(false)
{
}

void OptimizationLogger::ensure_initialized()
{
    if (!initialized_) {
        initialize();
    }
}

void OptimizationLogger::initialize()
{
    if (initialized_) {
        return;
    }

    std::error_code ec;
    std::filesystem::create_directories(save_root_, ec);
    if (ec) {
        std::cerr << "Warning: could not create save_root directory "
                  << save_root_ << " : " << ec.message() << "\n";
        // continue; opening meta.txt may still work if directory exists already
    }

    std::ofstream meta(save_root_ / "meta.txt");
    if (!meta) {
        std::cerr << "Warning: could not open meta.txt for header in "
                  << save_root_ << "\n";
        return; // leave initialized_ = false
    }

    // Header: geometry_id + blank line, no "(no runs)" here
    meta << geometry_id_ << "\n\n";

    initialized_ = true;
}

void OptimizationLogger::append_record(const OptRunRecord &rec)
{
    ensure_initialized();
    if (!initialized_) {
        return;
    }

    std::ofstream meta(save_root_ / "meta.txt", std::ios::app);
    if (!meta) {
        std::cerr << "Warning: could not open meta.txt for appending in "
                  << save_root_ << "\n";
        return;
    }

    write_single_opt_record_block(meta, save_root_, rec);
    meta.flush();
}


void run_metrics_only(const Config &cfg)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    const bool requested_civ =
        cfg.optimize.objective == "CIV" ||
        cfg.optimize.objective == "self_weighting";
    const bool requested_field_metrics =
        cfg.optimize.objective == "FieldSpread" ||
        cfg.optimize.objective == "self_weighting";
    const bool field_metrics_only =
        requested_field_metrics && !requested_civ;

    const auto targets = targets_from_save_root(cfg);
    struct MetricRecord
    {
        std::filesystem::path run_dir;
        OptimizationMetrics metrics;
        double top_ring_voltage = std::numeric_limits<double>::quiet_NaN();
        std::string wall_charge_state;
    };
    std::vector<MetricRecord> records;
    records.reserve(targets.size());

    auto load_preserved_metrics = [&](const std::filesystem::path &run_dir)
    {
        OptimizationMetrics preserved;
        const std::filesystem::path combined_path = run_dir / "metrics.yaml";
        if (!std::filesystem::exists(combined_path))
        {
            throw std::runtime_error(
                "Selected metrics reuse requires existing file: " +
                combined_path.string());
        }
        const YAML::Node combined = YAML::LoadFile(combined_path.string());
        if (requested_civ)
        {
            if (!combined["CIV"])
            {
                throw std::runtime_error(
                    "Existing metrics file has no CIV: " + combined_path.string());
            }
            preserved.CIV = combined["CIV"].as<double>();
            if (combined["tracing"])
            {
                preserved.traced_electrons =
                    combined["tracing"]["electron_count"].as<std::uint64_t>(0);
                preserved.tracing_seconds =
                    combined["tracing"]["elapsed_seconds"].as<double>(0.0);
            }
        }

        if (requested_field_metrics)
        {
            YAML::Node field = combined;
            const std::filesystem::path field_path =
                run_dir / "field_metrics.yaml";
            if (std::filesystem::exists(field_path))
            {
                field = YAML::LoadFile(field_path.string());
            }
            if (!field["FieldSpread"] || !field["field_homogeneity"])
            {
                throw std::runtime_error(
                    "Existing metrics have no field homogeneity data for " +
                    run_dir.string());
            }
            preserved.FieldSpread = field["FieldSpread"].as<double>();
            const YAML::Node h = field["field_homogeneity"];
            preserved.field_homogeneity.integration_volume =
                h["integration_volume"].as<double>();
            preserved.field_homogeneity.mean_magnitude =
                h["mean_magnitude_V_per_m"].as<double>();
            preserved.field_homogeneity.lower_quantile =
                h["lower_quantile_V_per_m"].as<double>();
            preserved.field_homogeneity.upper_quantile =
                h["upper_quantile_V_per_m"].as<double>();
            preserved.field_homogeneity.relative_rms_magnitude =
                h["relative_rms_magnitude"].as<double>();
            preserved.field_homogeneity.relative_rms_radial =
                h["relative_rms_radial_to_mean_abs_axial"].as<double>();
            preserved.field_homogeneity.relative_rms_axial =
                h["relative_rms_axial_to_mean_abs_axial"].as<double>();
            preserved.field_homogeneity.mean_axial =
                h["mean_axial_V_per_m"].as<double>();
            preserved.field_homogeneity.mean_abs_axial =
                h["mean_abs_axial_V_per_m"].as<double>();
        }
        return preserved;
    };

    for (const auto& run_dir : targets)
    {
        std::cout << "[METRICS] Loading results from " << run_dir << "\n";

        // Prefer the resolved leaf configuration saved with a sweep point. It
        // contains the actual voltage/wall-charge assignments and is the
        // authoritative record if the source sweep config changes later.
        Config local_cfg = cfg;
        const std::filesystem::path leaf_config = run_dir / "config.yaml";
        if (std::filesystem::exists(leaf_config))
        {
            // Replace a historical analysis method before validation. In
            // particular, completed simulations may predate the AdaptiveGrid
            // ban; their simulated state remains valid and reusable.
            YAML::Node leaf_node = YAML::LoadFile(leaf_config.string());
            leaf_node["optimize"]["civ_params"]["method"] =
                cfg.civ_params.method;
            local_cfg = Config::LoadFromNode(leaf_node, "metrics");
        }
        // The leaf config remains authoritative for the simulated state
        // (voltages, wall charge, mesh, and solver settings). Analysis controls
        // come from the config supplied to the metrics command so CIV methods
        // and sampling can be changed without rerunning the field solves.
        local_cfg.tracing_params = cfg.tracing_params;
        local_cfg.civ_params = cfg.civ_params;
        local_cfg.field_spread = cfg.field_spread;
        local_cfg.optimize.objective = cfg.optimize.objective;
        local_cfg.optimize.print_results = cfg.optimize.print_results;
        local_cfg.optimize.r_min = cfg.optimize.r_min;
        local_cfg.optimize.r_max = cfg.optimize.r_max;
        local_cfg.optimize.z_min = cfg.optimize.z_min;
        local_cfg.optimize.z_max = cfg.optimize.z_max;
        local_cfg.debug = cfg.debug;
        local_cfg.run_mode = "metrics";
        local_cfg.save_path = run_dir.string();

        double top_ring_voltage = std::numeric_limits<double>::quiet_NaN();
        const auto top_ring = local_cfg.boundaries.find("BC_FieldShapingRings_0");
        if (top_ring != local_cfg.boundaries.end()) {
            top_ring_voltage = top_ring->second.value;
        }
        const auto wall = local_cfg.boundaries.find("BC_WallCharge");
        const bool wall_is_charged =
            wall != local_cfg.boundaries.end() && wall->second.depth_dependent;
        const std::string wall_charge_state =
            wall_is_charged ? "charged" : "uncharged";

        if (cfg.run_selection.wall_charge != "all" &&
            wall_charge_state != cfg.run_selection.wall_charge)
        {
            if (rank == 0)
            {
                std::cout << "[METRICS] Reusing existing "
                          << wall_charge_state << " metrics for "
                          << run_dir.filename() << "\n";
                records.push_back(MetricRecord{
                    run_dir,
                    load_preserved_metrics(run_dir),
                    top_ring_voltage,
                    wall_charge_state});
            }
            continue;
        }

        SimulationResult result = load_results(local_cfg, run_dir);
        if (!result.success) {
            std::cerr << "[METRICS] Failed to load results: "
                      << result.error_message << "\n";
            std::exit(1);
        }

        OptimizationMetrics m = compute_metrics(local_cfg, result);
        if (rank == 0)
        {
            if (requested_civ)
            {
                std::cout << "[METRICS] CIV = " << m.CIV << "\n";
                std::cout << "[METRICS] Electron tracing: electrons="
                          << m.traced_electrons
                          << " elapsed_s=" << m.tracing_seconds << "\n";
            }
            if (requested_field_metrics)
            {
                std::cout << "[METRICS] FieldSpread = "
                          << m.FieldSpread << "\n";
                std::cout << "[METRICS] Field homogeneity: mean|E|="
                          << m.field_homogeneity.mean_magnitude
                          << " relative_rms|E|="
                          << m.field_homogeneity.relative_rms_magnitude
                          << " relative_rms_Er="
                          << m.field_homogeneity.relative_rms_radial
                          << " relative_rms_Ez="
                          << m.field_homogeneity.relative_rms_axial << "\n";
            }
        }

        if (rank == 0)
        {
            YAML::Node metrics;
            metrics["top_field_shaping_V"] = top_ring_voltage;
            metrics["wall_charge"] = wall_is_charged ? "charged" : "uncharged";
            if (wall_is_charged) {
                metrics["wall_charge_profile"]["z_bot"] = wall->second.z_bot;
                metrics["wall_charge_profile"]["z_top"] = wall->second.z_top;
                metrics["wall_charge_profile"]["value_bot_C_per_m2"] = wall->second.value_bot;
                metrics["wall_charge_profile"]["value_top_C_per_m2"] = wall->second.value_top;
            }
            if (requested_civ)
            {
                metrics["CIV"] = m.CIV;
                metrics["civ_method"] = local_cfg.civ_params.method;
                metrics["tracing"]["provider"] =
                    local_cfg.tracing_params.provider;
                metrics["tracing"]["ram_expensive"] =
                    local_cfg.tracing_params.ram_expensive;
                metrics["tracing"]["local_element_step"] =
                    local_cfg.tracing_params.local_element_step;
                metrics["tracing"]["c_step"] =
                    local_cfg.tracing_params.c_step;
                metrics["tracing"]["electron_count"] =
                    static_cast<unsigned long long>(m.traced_electrons);
                metrics["tracing"]["elapsed_seconds"] =
                    m.tracing_seconds;
                if (local_cfg.civ_params.method == "RandomSample") {
                    metrics["civ_num_seed_elements"] =
                        local_cfg.civ_params.num_seed_elements;
                } else {
                    metrics["civ_nr"] = local_cfg.civ_params.nr;
                    metrics["civ_nz"] = local_cfg.civ_params.nz;
                    if (std::isfinite(local_cfg.civ_params.sample_r_min)) {
                        metrics["civ_sample_r_min"] =
                            local_cfg.civ_params.sample_r_min;
                    }
                }
                metrics["civ_roi"]["r_min"] = local_cfg.tracing_params.r_min;
                metrics["civ_roi"]["r_max"] = local_cfg.tracing_params.r_max;
                metrics["civ_roi"]["z_min"] = local_cfg.tracing_params.z_min;
                metrics["civ_roi"]["z_max"] = local_cfg.tracing_params.z_max;
            }
            if (requested_field_metrics)
            {
                const auto &h = m.field_homogeneity;
                metrics["FieldSpread"] = m.FieldSpread;
                metrics["field_homogeneity"]["weighting"] =
                    local_cfg.solver.axisymmetric
                        ? "axisymmetric_volume_2pi_r_dr_dz"
                        : "cartesian_volume";
                metrics["field_homogeneity"]["integration_volume"] =
                    h.integration_volume;
                metrics["field_homogeneity"]["mean_magnitude_V_per_m"] =
                    h.mean_magnitude;
                metrics["field_homogeneity"]["lower_quantile_V_per_m"] =
                    h.lower_quantile;
                metrics["field_homogeneity"]["upper_quantile_V_per_m"] =
                    h.upper_quantile;
                metrics["field_homogeneity"]["lower_quantile_fraction"] =
                    local_cfg.field_spread.lower;
                metrics["field_homogeneity"]["upper_quantile_fraction"] =
                    local_cfg.field_spread.upper;
                metrics["field_homogeneity"]["relative_rms_magnitude"] =
                    h.relative_rms_magnitude;
                metrics["field_homogeneity"]["relative_rms_radial_to_mean_abs_axial"] =
                    h.relative_rms_radial;
                metrics["field_homogeneity"]["relative_rms_axial_to_mean_abs_axial"] =
                    h.relative_rms_axial;
                metrics["field_homogeneity"]["mean_axial_V_per_m"] =
                    h.mean_axial;
                metrics["field_homogeneity"]["mean_abs_axial_V_per_m"] =
                    h.mean_abs_axial;
                metrics["field_homogeneity"]["roi"]["r_min"] =
                    local_cfg.optimize.r_min;
                metrics["field_homogeneity"]["roi"]["r_max"] =
                    local_cfg.optimize.r_max;
                metrics["field_homogeneity"]["roi"]["z_min"] =
                    local_cfg.optimize.z_min;
                metrics["field_homogeneity"]["roi"]["z_max"] =
                    local_cfg.optimize.z_max;
            }

            const std::filesystem::path metrics_path =
                run_dir / (field_metrics_only
                    ? "field_metrics.yaml"
                    : "metrics.yaml");
            std::ofstream out(metrics_path);
            if (!out)
            {
                throw std::runtime_error(
                    "Unable to write metrics file " + metrics_path.string());
            }
            out << metrics << '\n';
            records.push_back(MetricRecord{
                run_dir,
                m,
                top_ring_voltage,
                wall_charge_state});
        }
    }

    if (rank == 0)
    {
        const std::filesystem::path summary_path =
            std::filesystem::path(cfg.save_path) /
            (field_metrics_only ? "field_metrics.csv" : "metrics.csv");
        std::ofstream summary(summary_path);
        if (!summary)
        {
            throw std::runtime_error(
                "Unable to write aggregate metrics file: " + summary_path.string());
        }
        summary << "run,top_field_shaping_V,wall_charge,CIV,FieldSpread,"
                   "FieldMeanMagnitude,FieldRMSRelative,RadialRMSRelative,"
                   "AxialRMSRelative\n";
        summary << std::setprecision(17);
        for (const auto &record : records)
        {
            summary << record.run_dir.filename().string() << ','
                    << record.top_ring_voltage << ','
                    << record.wall_charge_state << ',';
            if (requested_civ) { summary << record.metrics.CIV; }
            summary << ',';
            if (requested_field_metrics)
            {
                summary << record.metrics.FieldSpread << ','
                        << record.metrics.field_homogeneity.mean_magnitude << ','
                        << record.metrics.field_homogeneity.relative_rms_magnitude << ','
                        << record.metrics.field_homogeneity.relative_rms_radial << ','
                        << record.metrics.field_homogeneity.relative_rms_axial;
            }
            else
            {
                summary << ",,,,";
            }
            summary << '\n';
        }
        std::cout << "[METRICS] Wrote aggregate summary to " << summary_path << "\n";

        if (requested_civ)
        {
            const std::filesystem::path civ_path =
                std::filesystem::path(cfg.save_path) / "civ.csv";
            std::ofstream civ(civ_path);
            if (!civ) {
                throw std::runtime_error(
                    "Unable to write aggregate CIV file: " + civ_path.string());
            }
            civ << "top_field_shaping_V,wall_charge,CIV\n"
                << std::setprecision(17);
            for (const auto &record : records) {
                civ << record.top_ring_voltage << ','
                    << record.wall_charge_state << ','
                    << record.metrics.CIV << '\n';
            }
        }

        if (requested_field_metrics && !field_metrics_only)
        {
            const std::filesystem::path field_path =
                std::filesystem::path(cfg.save_path) / "field_metrics.csv";
            std::ofstream field(field_path);
            if (!field) {
                throw std::runtime_error(
                    "Unable to write aggregate field metrics file: " +
                    field_path.string());
            }
            field << "run,top_field_shaping_V,wall_charge,FieldSpread,"
                     "FieldMeanMagnitude,FieldRMSRelative,RadialRMSRelative,"
                     "AxialRMSRelative\n"
                  << std::setprecision(17);
            for (const auto &record : records) {
                field << record.run_dir.filename().string() << ','
                      << record.top_ring_voltage << ','
                      << record.wall_charge_state << ','
                      << record.metrics.FieldSpread << ','
                      << record.metrics.field_homogeneity.mean_magnitude << ','
                      << record.metrics.field_homogeneity.relative_rms_magnitude << ','
                      << record.metrics.field_homogeneity.relative_rms_radial << ','
                      << record.metrics.field_homogeneity.relative_rms_axial << '\n';
            }
        }
    }
    MPI_Barrier(comm);
}
