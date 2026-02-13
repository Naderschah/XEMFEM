#include "path_handler.h"

void ensure_directory(const Config cfg)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fs::path dir_path(cfg.save_path);
    bool delete_contents = cfg.delete_files_present;

    int error_flag = 0;

    if (rank == 0)
    {
        try
        {
            if (fs::exists(dir_path))
            {
                if (!fs::is_directory(dir_path))
                {
                    throw std::runtime_error("Path exists but is not a directory: " + dir_path.string());
                }

                if (delete_contents)
                {
                    for (const auto& entry : fs::directory_iterator(dir_path))
                    {
                        fs::remove_all(entry.path());
                    }
                }
            }
            else
            {
                fs::create_directories(dir_path);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Rank 0 filesystem error: " << e.what() << std::endl;
            error_flag = 1;
        }
    }

    // Broadcast error status to all ranks
    MPI_Bcast(&error_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Ensure all ranks wait until directory handling is complete
    MPI_Barrier(MPI_COMM_WORLD);

    if (error_flag)
    {
        throw std::runtime_error("Directory setup failed on rank 0.");
    }
}


fs::path make_run_folder(const Config& cfg, int run_id, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::ostringstream ss;
    ss << "run_" << std::setw(4) << std::setfill('0') << run_id;

    fs::path run_path = fs::path(cfg.save_path) / ss.str();

    int error_flag = 0;

    if (rank == 0)
    {
        try {
            if (!fs::exists(run_path)) {
                fs::create_directories(run_path);
            }
        } catch (...) {
            error_flag = 1;
        }
    }

    MPI_Bcast(&error_flag, 1, MPI_INT, 0, comm);
    MPI_Barrier(comm);

    if (error_flag)
        throw std::runtime_error("Failed to create run folder.");

    return run_path;
}

fs::path make_run_folder(const YAML::Node yaml_config, int run_id, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::ostringstream ss;
    ss << "run_" << std::setw(4) << std::setfill('0') << run_id;

    fs::path run_path = fs::path(yaml_config["save_path"].as<std::string>()) / ss.str();

    int error_flag = 0;

    if (rank == 0)
    {
        try {
            if (!fs::exists(run_path)) {
                fs::create_directories(run_path);
            }
        } catch (...) {
            error_flag = 1;
        }
    }

    MPI_Bcast(&error_flag, 1, MPI_INT, 0, comm);
    MPI_Barrier(comm);

    if (error_flag)
        throw std::runtime_error("Failed to create run folder.");

    return run_path;
}



std::vector<fs::path> targets_from_save_root(const Config& cfg)
{
    const fs::path save_root(cfg.save_path);

    if (!fs::exists(save_root) || !fs::is_directory(save_root)) {
        throw std::runtime_error("metrics_targets_from_save_root: save_root does not exist or is not a directory: "
                                 + save_root.string());
    }

    std::vector<fs::path> run_dirs;

    for (const auto& entry : fs::directory_iterator(save_root))
    {
        if (!entry.is_directory())
            continue;

        const std::string name = entry.path().filename().string();

        // Match run_XXXX
        if (name.rfind("run_", 0) == 0) {
            run_dirs.push_back(entry.path());
        }
    }

    if (run_dirs.empty()) {
        return { save_root };
    }

    std::sort(run_dirs.begin(), run_dirs.end());
    return run_dirs;
}