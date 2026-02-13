#pragma once

#include <mpi.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include "Config.h"
#include <yaml-cpp/yaml.h>

namespace fs = std::filesystem;

void ensure_directory(const Config cfg); 
fs::path make_run_folder(const Config& cfg, int run_id, MPI_Comm comm = MPI_COMM_WORLD);
fs::path make_run_folder(const YAML::Node yaml_config, int run_id, MPI_Comm comm = MPI_COMM_WORLD);
std::vector<fs::path> targets_from_save_root(const Config& cfg);