#include "sweeps.h"

// -----------------------------------------------------------------------------
// Recursive sweep over Config
// -----------------------------------------------------------------------------
void sweep_recursive_cfg(const Config &base_cfg,
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
void write_sweep_meta(const std::string &geometry_id,
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