#include "optimization.h"

// ===================================== Main Loop ==================================== 

void run_optimization(const Config &init_cfg, std::string config_str)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, world_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    const auto &opt = init_cfg.optimize;
    if (!opt.enabled) { return; }

    if (opt.variables.empty()) {
        throw std::runtime_error("optimize.enabled=true but no variables specified");
    }

    // Build objective function (metrics -> scalar) on all ranks (must be consistent)
    ObjectiveFn objective_fn = make_objective_function(opt);

    // Storage and logger: only rank 0 owns these
    std::vector<OptRunRecord> records;
    OptimizationLogger *logger_ptr = nullptr;
    std::unique_ptr<OptimizationLogger> logger_owner;

    std::filesystem::path save_root(init_cfg.save_path);

    if (rank == 0) {
        records.reserve(opt.max_fun_evals > 0 ? opt.max_fun_evals : 100);
        logger_owner = std::make_unique<OptimizationLogger>(init_cfg.geometry_id, save_root);
        logger_owner->initialize();
        logger_ptr = logger_owner.get();
    } else {
        // dummy storage on non-root (evaluate may accept refs)
        records.clear();
        records.shrink_to_fit();
    }

    // Build initial point x0 (deterministic; do it on rank 0, broadcast)
    std::vector<double> x0;
    const int nvars = static_cast<int>(opt.variables.size());

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

    // Broadcast x0 so all ranks agree (mainly for sanity / potential debug)
    MPI_Bcast(x0.data(), nvars, MPI_DOUBLE, 0, comm);

    // Command protocol:
    //   cmd = 1 -> evaluate one point (followed by eval_index + x vector)
    //   cmd = 0 -> stop service loop
    const int CMD_EVAL = 1;
    const int CMD_STOP = 0;

    // Non-root ranks enter a service loop:
    // They do NOT run Nelder–Mead. They just receive x from rank 0 and participate
    // in the collective simulation inside evaluate_one_optimization_point.
    if (rank != 0)
    {
        while (true)
        {
            int cmd = CMD_STOP;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);

            if (cmd == CMD_STOP) {
                break;
            }

            std::uint64_t eval_index_u64 = 0;
            MPI_Bcast(&eval_index_u64, 1, MPI_UINT64_T, 0, comm);

            std::vector<double> x(nvars, 0.0);
            MPI_Bcast(x.data(), nvars, MPI_DOUBLE, 0, comm);

            // Participate in evaluation (collective MPI work happens inside)
            // Pass dummy logger; records should ideally only be appended on rank 0
            // (evaluate_one_optimization_point should guard on logger!=nullptr or rank==0).
            (void)evaluate_one_optimization_point(
                config_str,
                opt,
                x,
                static_cast<std::size_t>(eval_index_u64),
                objective_fn,
                records,
                /*logger*/ nullptr
            );
        }

        // Optionally receive final x_opt (useful if downstream needs it)
        // Rank 0 will broadcast it.
        int xopt_n = 0;
        MPI_Bcast(&xopt_n, 1, MPI_INT, 0, comm);
        std::vector<double> x_opt(xopt_n, 0.0);
        if (xopt_n > 0) {
            MPI_Bcast(x_opt.data(), xopt_n, MPI_DOUBLE, 0, comm);
        }

        return;
    }

    // ---------------- Rank 0 runs Nelder–Mead ----------------

    std::size_t eval_counter = 0;

    auto nm_func = [&](const std::vector<double> &x) -> double
    {
        // Broadcast "evaluate" command
        int cmd = CMD_EVAL;
        MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);

        // Broadcast eval index to keep run naming consistent across ranks
        std::uint64_t eval_index_u64 = static_cast<std::uint64_t>(eval_counter);
        MPI_Bcast(&eval_index_u64, 1, MPI_UINT64_T, 0, comm);

        // Broadcast x
        MFEM_VERIFY(static_cast<int>(x.size()) == nvars, "nm_func: x size mismatch");
        MPI_Bcast(const_cast<double*>(x.data()), nvars, MPI_DOUBLE, 0, comm);

        // Evaluate (all ranks will participate; rank 0 gets the returned objective)
        double f = evaluate_one_optimization_point(
            config_str,
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

    // Broadcast final x_opt to all ranks (optional but usually helpful)
    {
        int xopt_n = static_cast<int>(x_opt.size());
        MPI_Bcast(&xopt_n, 1, MPI_INT, 0, comm);
        if (xopt_n > 0) {
            MPI_Bcast(x_opt.data(), xopt_n, MPI_DOUBLE, 0, comm);
        }
    }

    // Find the best recorded run (rank 0 only)
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

    // Write meta.txt summarizing all optimization runs (rank 0 only)
    write_optimization_meta(init_cfg.geometry_id, save_root, records);
}


RunAndMetricsResult run_simulation_and_save(YAML::Node yaml_config, std::size_t eval_index)
{
    RunAndMetricsResult out;

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, world_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    // -------------------------------------------------------------------------
    // 0) Determine save_root from YAML (supports root["save_path"] or root["config"]["save_path"])
    // -------------------------------------------------------------------------
    std::string save_root_str;
    if (yaml_config["config"] && yaml_config["config"].IsMap() && yaml_config["config"]["save_path"])
        save_root_str = yaml_config["config"]["save_path"].as<std::string>();
    else if (yaml_config["save_path"])
        save_root_str = yaml_config["save_path"].as<std::string>();
    else
        throw std::runtime_error("run_simulation_and_save: YAML missing save_path (root.save_path or root.config.save_path)");

    std::filesystem::path save_root(save_root_str);

    // -------------------------------------------------------------------------
    // 1) Create output directories (rank 0), then barrier
    // -------------------------------------------------------------------------
    std::error_code ec;

    if (rank == 0) {
        std::filesystem::create_directories(save_root, ec);
        if (ec) {
            std::cerr << "Warning: could not create save_root directory "
                      << save_root << " : " << ec.message() << "\n";
            ec.clear();
        }
    }
    MPI_Barrier(comm);

    // Make run name
    std::size_t display_index = eval_index + 1;
    std::ostringstream run_name;
    run_name << "run_" << std::setw(4) << std::setfill('0') << display_index;
    out.run_dir_name = run_name.str();

    std::filesystem::path run_dir = save_root / out.run_dir_name;

    if (rank == 0) {
        std::filesystem::create_directories(run_dir, ec);
        if (ec) {
            std::cerr << "Warning: could not create run directory "
                      << run_dir << " : " << ec.message() << "\n";
            ec.clear();
        }
    }
    MPI_Barrier(comm);

    // -------------------------------------------------------------------------
    // 2) Mutate YAML for this run (all ranks do the same deterministic changes)
    // -------------------------------------------------------------------------
    {
        const std::string run_dir_str = run_dir.string();
        if (yaml_config["config"] && yaml_config["config"].IsMap())
            yaml_config["config"]["save_path"] = run_dir_str;
        else
            yaml_config["save_path"] = run_dir_str;
    }

    // Record MPI size in YAML (ensure 'mpi' is a map)
    {
        yaml_config["mpi"]["ranks"] = world_size;
    }

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

    SimulationResult result = run_simulation(cfg_ptr, model_path);
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

double evaluate_one_optimization_point(const std::string &config_str,
                                       const OptimizationSettings &opt,
                                       const std::vector<double> &x,
                                       std::size_t eval_index,
                                       const ObjectiveFn &objective_fn,
                                       std::vector<OptRunRecord> &records,
                                       OptimizationLogger *logger)
{
    std::cout << "[OPTIMIZATION] Running iteration " << eval_index << std::endl;
    // 1. Build config for this eval
    auto root = YAML::Load(config_str);

    // apply optimization variables
    apply_opt_vars(root, opt, x);

    // apply field cage network (if enabled)
    apply_fieldcage_network(root);

    // 2. Run simulation and save outputs
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

double compute_civ(const Config &cfg, const SimulationResult &result)
/*
Compute Charge Insensitive Volume on the mesh
Ie 1-(Volume where electrons reach liquid gas interface) / (Total Volume)
*/
{
    std::chrono::steady_clock::time_point t_start;
    if (cfg.debug.debug)
    {
        t_start = std::chrono::steady_clock::now();
        std::cout << "[DEBUG:CIV] Timing: start CIV" << std::endl;
    }

    const std::string &method = cfg.civ_params.method; // "InformedSweep" or "RandomSample"
    double civ = 1.0;

    if (method == "RandomSample")
    {
        civ = ComputeCIV_RandomSample(cfg, result);
    }
    else if (method == "Grid")
    {
        civ = ComputeCIV_FixedGrid(cfg, result);
    }
    else if (method == "AdaptiveGrid")
    {
        civ = ComputeCIV_AdaptiveGrid(cfg, result);
    }
    else if (method == "RowSweep")
    {
        civ = ComputeCIV_RowSweep(cfg, result);
    }
    else
    {
        std::cerr << "[WARN:OPTIMIZATION] Unknown CIV method '" << method << "',\n";
        civ = 1.0;
    }

    if (cfg.debug.debug)
    {
        auto t_end = std::chrono::steady_clock::now();
        auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
        std::cout << "[DEBUG:CIV] Timing: end CIV (" << dt << " ms)" << std::endl;
    }

    return civ;
}



OptimizationMetrics compute_metrics(const Config &cfg, const SimulationResult &result)
{
    OptimizationMetrics m;
    if ((cfg.optimize.objective == "CIV") || (cfg.optimize.objective == "self_weighting"))
    {
        m.CIV = compute_civ(cfg, result);
        if (cfg.optimize.print_results)
        {
            std::cout << "[OPT] CIV: " << m.CIV << std::endl;
        }
    }
    

    if ((cfg.optimize.objective == "FieldSpread") || (cfg.optimize.objective == "self_weighting"))
    {
        m.FieldSpread = computeFieldSpreadMetric(cfg, result);
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
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (rank != 0) {
        return;
    }

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
        meta.flush();
        MPI_Barrier(comm);
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

    meta.flush();
    MPI_Barrier(comm);
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


// Helper code for only running optimization no mesh creation

static inline std::string rtrim(const std::string &s)
{
    std::string::size_type end = s.find_last_not_of(" \t\r\n");
    if (end == std::string::npos) return "";
    return s.substr(0, end + 1);
}

static std::vector<std::string> read_root_level_runs_from_meta(const std::filesystem::path &save_root)
{
    std::vector<std::string> runs;

    std::filesystem::path meta_path = save_root / "meta.txt";
    if (!std::filesystem::exists(meta_path)) {
        return runs; // empty → no meta.txt
    }

    std::ifstream meta(meta_path);
    if (!meta) {
        std::cerr << "Warning: could not open meta.txt in " << save_root << "\n";
        return runs;
    }

    std::string line;
    while (std::getline(meta, line)) {
        // Ignore empty lines
        if (line.empty()) continue;

        // Ignore lines starting with whitespace (indented blocks)
        if (!line.empty() && (line[0] == ' ' || line[0] == '\t')) {
            continue;
        }

        // We only care about root-level lines like "run_0001:"
        std::string trimmed = rtrim(line);
        if (trimmed.size() < 2) continue;

        if (trimmed.back() == ':') {
            std::string name = trimmed.substr(0, trimmed.size() - 1);
            if (!name.empty()) {
                runs.push_back(name);
            }
        }
    }

    return runs;
}

// TODO add option for multiple runs 
void run_metrics_only(const Config &cfg)
{
    std::filesystem::path save_root(cfg.save_path);

    std::filesystem::path run_dir;

    auto runs = read_root_level_runs_from_meta(save_root);
    if (runs.empty()) {
        run_dir = save_root;
        std::cout << "[METRICS] meta.txt not found or empty; "
                  << "using save_root directly: " << run_dir << "\n";
    } else if (runs.size() == 1) {
        run_dir = save_root / runs.front();
        std::cout << "[METRICS] Found single run in meta.txt: "
                  << runs.front() << " → " << run_dir << "\n";
    } else {
        std::cerr << "[METRICS] meta.txt contains multiple runs under "
                  << save_root << ":\n";
        for (const auto &r : runs) {
            std::cerr << "  - " << r << "\n";
        }
        std::cerr << "Please specify which run to use (CLI/config).\n";
        std::exit(1);
    }

    std::cout << "[METRICS] Loading results from " << run_dir << "\n";

    SimulationResult result = load_results(cfg, run_dir);
    if (!result.success) {
        std::cerr << "[METRICS] Failed to load results: "
                  << result.error_message << "\n";
        std::exit(1);
    }

    OptimizationMetrics m = compute_metrics(cfg, result);
    std::cout << "[METRICS] CIV = " << m.CIV << "\n";
    std::cout << "[METRICS] FieldSpread = " << m.FieldSpread << "\n";
}