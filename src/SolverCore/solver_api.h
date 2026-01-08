#pragma once

#include <memory>
#include <filesystem>
#include <string>
#include <vector>
#include <utility>
#include <cstddef>


#include "Config.h"
#include "mfem.hpp"

struct SimulationResult
{
    bool success = false;
    std::string error_message;

    std::unique_ptr<mfem::ParMesh>                 mesh;
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    std::unique_ptr<mfem::ParFiniteElementSpace>   pfes;

    std::unique_ptr<mfem::ParGridFunction> V;
    std::unique_ptr<mfem::GridFunction> E;
    std::unique_ptr<mfem::GridFunction> Emag;

};

// Core solver entry point (no CLI, no file I/O)
SimulationResult run_simulation(
    std::shared_ptr<Config> cfg,
    const std::filesystem::path& model_path
);

void save_results(const SimulationResult &result, const std::filesystem::path &root_path);

SimulationResult load_results(const Config &cfg, const std::filesystem::path &root_path);

std::string run_one(const Config& cfg, const std::vector<std::pair<std::string, std::string>>& active_params, std::size_t run_index);