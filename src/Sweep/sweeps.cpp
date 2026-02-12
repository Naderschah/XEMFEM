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
    int rank = 0, world_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    // Helper: parse a string into a YAML scalar with best-effort typing
    auto parse_scalar = [](const std::string &s) -> YAML::Node {
        if (s == "true" || s == "True" || s == "TRUE")   return YAML::Node(true);
        if (s == "false" || s == "False" || s == "FALSE") return YAML::Node(false);

        {
            std::size_t pos = 0;
            try {
                long long v = std::stoll(s, &pos);
                if (pos == s.size()) return YAML::Node(v);
            } catch (...) {}
        }
        {
            std::size_t pos = 0;
            try {
                double v = std::stod(s, &pos);
                if (pos == s.size()) return YAML::Node(v);
            } catch (...) {}
        }
        return YAML::Node(s);
    };

    // Helper: set root[path] = value_node, where path is dot-separated
    auto set_by_dot_path = [&](YAML::Node &root, const std::string &path, const YAML::Node &value_node) {
        std::stringstream ss(path);
        std::string key;
        std::vector<std::string> keys;
        while (std::getline(ss, key, '.')) {
            if (!key.empty()) keys.push_back(key);
        }
        if (keys.empty()) {
            throw std::runtime_error("sweep_recursive_cfg: empty assignment path");
        }

        YAML::Node node = root;
        for (std::size_t k = 0; k + 1 < keys.size(); ++k) {
            const std::string &kk = keys[k];
            if (!node[kk]) {
                throw std::runtime_error(
                    "sweep_recursive_cfg: invalid path '" + path + "' (missing key '" + kk + "')");
            }
            node = node[kk];
        }
        node[keys.back()] = value_node;
    };

    // Command protocol:
    //   Root rank enumerates all sweep points and broadcasts the leaf work.
    //   All ranks call run_one for each leaf (collective heavy compute).
    const int CMD_EVAL = 1;
    const int CMD_STOP = 0;

    // ---------------- Non-root ranks: service loop ----------------
    if (rank != 0)
    {
        while (true)
        {
            int cmd = CMD_STOP;
            MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);

            if (cmd == CMD_STOP) {
                break;
            }

            // Receive run_counter (for consistent run naming)
            std::uint64_t run_counter_u64 = 0;
            MPI_Bcast(&run_counter_u64, 1, MPI_UINT64_T, 0, comm);

            // Receive YAML string length + bytes
            std::uint64_t yaml_len_u64 = 0;
            MPI_Bcast(&yaml_len_u64, 1, MPI_UINT64_T, 0, comm);

            std::string run_config_str;
            run_config_str.resize(static_cast<std::size_t>(yaml_len_u64));
            if (yaml_len_u64 > 0) {
                MPI_Bcast(run_config_str.data(), static_cast<int>(yaml_len_u64), MPI_CHAR, 0, comm);
            }

            // Receive active_params for metadata (label/value pairs)
            std::uint64_t nparams_u64 = 0;
            MPI_Bcast(&nparams_u64, 1, MPI_UINT64_T, 0, comm);
            const std::size_t nparams = static_cast<std::size_t>(nparams_u64);

            std::vector<std::pair<std::string, std::string>> active_params_recv;
            active_params_recv.reserve(nparams);

            auto bcast_string = [&](std::string &s) {
                std::uint64_t len_u64 = 0;
                MPI_Bcast(&len_u64, 1, MPI_UINT64_T, 0, comm);
                s.resize(static_cast<std::size_t>(len_u64));
                if (len_u64 > 0) {
                    MPI_Bcast(s.data(), static_cast<int>(len_u64), MPI_CHAR, 0, comm);
                }
            };

            for (std::size_t i = 0; i < nparams; ++i) {
                std::string k, v;
                bcast_string(k);
                bcast_string(v);
                active_params_recv.emplace_back(std::move(k), std::move(v));
            }

            // Collective heavy compute: all ranks call run_one
            YAML::Node yaml_root = YAML::Load(run_config_str);
            (void)run_one(base_cfg,
                          yaml_root,
                          active_params_recv,
                          static_cast<std::size_t>(run_counter_u64), 
                          true);
        }

        return;
    }

    // ---------------- Rank 0: enumerates sweep points ----------------

    if (idx == sweeps.size())
    {
        // Leaf: build per-run YAML string
        YAML::Node root = YAML::Load(config_str);
        //for (const auto &a : assignments) {
        //    set_by_dot_path(root, a.path, parse_scalar(a.value));
        //}

        YAML::Emitter out;
        out << root;
        std::string run_config_str = out.c_str();

        // Broadcast "evaluate" command and payload to workers
        int cmd = CMD_EVAL;
        MPI_Bcast(&cmd, 1, MPI_INT, 0, comm);

        std::uint64_t run_counter_u64 = static_cast<std::uint64_t>(run_counter);
        MPI_Bcast(&run_counter_u64, 1, MPI_UINT64_T, 0, comm);

        std::uint64_t yaml_len_u64 = static_cast<std::uint64_t>(run_config_str.size());
        MPI_Bcast(&yaml_len_u64, 1, MPI_UINT64_T, 0, comm);
        if (yaml_len_u64 > 0) {
            MPI_Bcast(run_config_str.data(), static_cast<int>(yaml_len_u64), MPI_CHAR, 0, comm);
        }

        // Broadcast active_params so all ranks call run_one with identical metadata
        auto bcast_string_root = [&](const std::string &s) {
            std::uint64_t len_u64 = static_cast<std::uint64_t>(s.size());
            MPI_Bcast(&len_u64, 1, MPI_UINT64_T, 0, comm);
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

        // Collective compute
        YAML::Node yaml_root = YAML::Load(run_config_str);
        std::string run_dir_name = run_one(base_cfg, yaml_root, active_params, run_counter, true);

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

                sweep_recursive_cfg(base_cfg, config_str, sweeps, idx + 1,
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

                    sweep_recursive_cfg(base_cfg, config_str, sweeps, idx + 1,
                                        active_params, assignments,
                                        run_counter, records);

                    assignments.pop_back();
                    active_params.pop_back();
                }
            }
            break;
        }
    }

    // Only rank 0 reaches here. Caller (rank 0) should broadcast CMD_STOP after sweep completes.
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

