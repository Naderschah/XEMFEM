#include "sweeps.h"

// -----------------------------------------------------------------------------
// Recursive sweep over Config
// -----------------------------------------------------------------------------
void sweep_recursive_cfg(const Config &base_cfg,
                         const std::string config_str,
                         const std::vector<SweepEntry> &sweeps,
                         std::size_t idx,
                         std::vector<std::pair<std::string, std::string>> &active_params, // (label, value)
                         std::vector<Assignment> &assignments,                           // (path, value)
                         std::size_t &run_counter,
                         std::vector<RunRecord> &records)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    const int CMD_EVAL = 1;
    const int CMD_STOP = 0;
    const std::uint32_t SWEEP_MAGIC = 0x53E3E9A1u;

    // Parse string -> YAML scalar with best-effort typing (bool/int/double/string).
    auto parse_scalar = [](const std::string &s) -> YAML::Node {
        if (s == "true" || s == "True" || s == "TRUE")    return YAML::Node(true);
        if (s == "false" || s == "False" || s == "FALSE") return YAML::Node(false);

        { std::size_t pos = 0; try { long long v = std::stoll(s, &pos); if (pos == s.size()) return YAML::Node(v); } catch (...) {} }
        { std::size_t pos = 0; try { double v = std::stod(s, &pos);    if (pos == s.size()) return YAML::Node(v); } catch (...) {} }

        return YAML::Node(s);
    };

    // Set root[path] = value_node, where path is dot-separated.
    auto set_by_dot_path = [&](YAML::Node &root,
                               const std::string &path,
                               const YAML::Node &value_node)
    {
        std::stringstream ss(path);
        std::string key;
        std::vector<std::string> keys;
        while (std::getline(ss, key, '.')) if (!key.empty()) keys.push_back(key);

        if (keys.empty())
            throw std::runtime_error("sweep_recursive_cfg: empty assignment path");
        if (!root || !root.IsMap())
            throw std::runtime_error("sweep_recursive_cfg: root is not a map for path '" + path + "'");

        YAML::Node cur = root;
        for (std::size_t i = 0; i + 1 < keys.size(); ++i) {
            const std::string &kk = keys[i];
            if (!cur[kk])
                throw std::runtime_error("sweep_recursive_cfg: invalid path '" + path + "' (missing key '" + kk + "')");
            if (!cur[kk].IsMap())
                throw std::runtime_error("sweep_recursive_cfg: path '" + path + "' traverses non-map key '" + kk + "'");
            cur.reset(cur[kk]);
        }
        cur[keys.back()] = value_node;
    };

    // ---------------------------------------------------------------------
    // Workers: enter the service loop exactly once (top-level call only).
    // Root broadcasts: MAGIC, CMD, then payload for each leaf.
    // ---------------------------------------------------------------------
    if (rank != 0)
    {
        if (idx != 0) return;

        while (true)
        {
            std::uint32_t magic = 0;
            MPI_Bcast(&magic, 1, MPI_UINT32_T, 0, comm);
            if (magic != SWEEP_MAGIC) {
                MPI_Abort(comm, 1);
            }

            int cmd = CMD_STOP;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);
            if (cmd == CMD_STOP) break;

            std::uint64_t run_counter_u64 = 0;
            MPI_Bcast(&run_counter_u64, 1, MPI_UINT64_T, 0, comm);

            std::uint64_t yaml_len_u64 = 0;
            MPI_Bcast(&yaml_len_u64, 1, MPI_UINT64_T, 0, comm);

            if (yaml_len_u64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                MPI_Abort(comm, 2);
            }

            std::string run_config_str;
            run_config_str.resize(static_cast<std::size_t>(yaml_len_u64));
            if (yaml_len_u64 > 0) {
                MPI_Bcast(run_config_str.data(), static_cast<int>(yaml_len_u64), MPI_CHAR, 0, comm);
            }

            std::uint64_t nparams_u64 = 0;
            MPI_Bcast(&nparams_u64, 1, MPI_UINT64_T, 0, comm);
            const std::size_t nparams = static_cast<std::size_t>(nparams_u64);

            auto bcast_string = [&](std::string &s) {
                std::uint64_t len_u64 = 0;
                MPI_Bcast(&len_u64, 1, MPI_UINT64_T, 0, comm);

                if (len_u64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                    MPI_Abort(comm, 3);
                }

                s.resize(static_cast<std::size_t>(len_u64));
                if (len_u64 > 0) {
                    MPI_Bcast(s.data(), static_cast<int>(len_u64), MPI_CHAR, 0, comm);
                }
            };

            std::vector<std::pair<std::string, std::string>> active_params_recv;
            active_params_recv.reserve(nparams);

            for (std::size_t i = 0; i < nparams; ++i) {
                std::string k, v;
                bcast_string(k);
                bcast_string(v);
                active_params_recv.emplace_back(std::move(k), std::move(v));
            }

            // make_run_folder() uses MPI_COMM_WORLD internally, so all ranks must call it here.
            YAML::Node yaml_root = YAML::Load(run_config_str);
            fs::path run_path = make_run_folder(yaml_root, static_cast<int>(run_counter_u64));
            yaml_root["save_path"] = run_path.string();

            (void)run_one(base_cfg,
                          yaml_root,
                          active_params_recv,
                          static_cast<std::size_t>(run_counter_u64),
                          true);

            // Keep protocol lockstep between leaf evaluations.
            MPI_Barrier(comm);
        }

        return;
    }

    // ---------------------------------------------------------------------
    // Root: recurse to enumerate sweep points; at each leaf broadcast work.
    // ---------------------------------------------------------------------
    if (idx == sweeps.size())
    {
        YAML::Node root = YAML::Load(config_str);
        for (const auto &a : assignments) {
            set_by_dot_path(root, a.path, parse_scalar(a.value));
        }

        // Do not call make_run_folder() before broadcasting: it contains WORLD collectives.
        // All ranks call make_run_folder() after receiving the YAML payload.
        YAML::Emitter out;
        out << root;
        std::string run_config_str = out.c_str();

        // Header: MAGIC then CMD (workers validate stream / allow STOP).
        {
            std::uint32_t magic = SWEEP_MAGIC;
            MPI_Bcast(&magic, 1, MPI_UINT32_T, 0, comm);

            int cmd = CMD_EVAL;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);
        }

        std::uint64_t run_counter_u64 = static_cast<std::uint64_t>(run_counter);
        MPI_Bcast(&run_counter_u64, 1, MPI_UINT64_T, 0, comm);

        std::uint64_t yaml_len_u64 = static_cast<std::uint64_t>(run_config_str.size());
        MPI_Bcast(&yaml_len_u64, 1, MPI_UINT64_T, 0, comm);

        if (yaml_len_u64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            MPI_Abort(comm, 2);
        }

        if (yaml_len_u64 > 0) {
            MPI_Bcast(run_config_str.data(), static_cast<int>(yaml_len_u64), MPI_CHAR, 0, comm);
        }

        auto bcast_string_root = [&](const std::string &s) {
            std::uint64_t len_u64 = static_cast<std::uint64_t>(s.size());
            MPI_Bcast(&len_u64, 1, MPI_UINT64_T, 0, comm);

            if (len_u64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                MPI_Abort(comm, 3);
            }

            if (len_u64 > 0) {
                MPI_Bcast(const_cast<char*>(s.data()), static_cast<int>(len_u64), MPI_CHAR, 0, comm);
            }
        };

        std::uint64_t nparams_u64 = static_cast<std::uint64_t>(active_params.size());
        MPI_Bcast(&nparams_u64, 1, MPI_UINT64_T, 0, comm);

        for (const auto &kv : active_params) {
            bcast_string_root(kv.first);
            bcast_string_root(kv.second);
        }

        // All ranks must call make_run_folder() (WORLD collectives inside).
        YAML::Node yaml_root = YAML::Load(run_config_str);
        fs::path run_path = make_run_folder(yaml_root, static_cast<int>(run_counter_u64));
        yaml_root["save_path"] = run_path.string();

        std::string run_dir_name = run_one(base_cfg, yaml_root, active_params, run_counter, true);

        MPI_Barrier(comm);

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

                sweep_recursive_cfg(base_cfg, config_str, sweeps, idx + 1,
                                    active_params, assignments, run_counter, records);

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

                sweep_recursive_cfg(base_cfg, config_str, sweeps, idx + 1,
                                    active_params, assignments, run_counter, records);

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

                    sweep_recursive_cfg(base_cfg, config_str, sweeps, idx + 1,
                                        active_params, assignments, run_counter, records);

                    assignments.pop_back();
                    active_params.pop_back();
                }
            }
            break;
        }

        case SweepEntry::Kind::Fixed:
        {
            for (const auto &cfg : sw.configs)
            {
                active_params.emplace_back(label, cfg.label);

                const std::size_t n_added = cfg.assigns.size();
                for (const auto &a : cfg.assigns) assignments.push_back(a);

                sweep_recursive_cfg(base_cfg, config_str, sweeps, idx + 1,
                                    active_params, assignments, run_counter, records);

                for (std::size_t k = 0; k < n_added; ++k) assignments.pop_back();
                active_params.pop_back();
            }
            break;
        }
    }

    // Top-level only: send STOP. Do not add a barrier here; workers are waiting in Bcast(magic).
    if (idx == 0)
    {
        std::uint32_t magic = SWEEP_MAGIC;
        MPI_Bcast(&magic, 1, MPI_UINT32_T, 0, comm);

        int cmd = CMD_STOP;
        MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);
    }
}




// -----------------------------------------------------------------------------
// Global meta.txt in save_root
// -----------------------------------------------------------------------------
void write_sweep_meta(const std::string &geometry_id,
                      const std::filesystem::path &save_root,
                      const std::vector<RunRecord> &records)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0)
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
        }
        else
        {
            meta << geometry_id << "\n\n";

            if (records.empty())
            {
                meta << "(no runs)\n";
            }
            else
            {
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

            meta.flush();
        }
    }

    MPI_Barrier(comm);
}

