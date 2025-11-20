#pragma once

#include <memory>
#include <filesystem>
#include <string>

#include "Config.h"
#include "mfem.hpp"

struct SimulationResult
{
    bool success = false;
    std::string error_message;

    // Numerical outputs you want to reuse in Sweep
    std::unique_ptr<mfem::Mesh> mesh;
    std::unique_ptr<mfem::GridFunction> V;     // potential
    std::unique_ptr<mfem::GridFunction> E;     // electric field (vector GridFunction)
    std::unique_ptr<mfem::GridFunction> Emag;  // |E| as scalar GridFunction

};

// Core solver entry point (no CLI, no file I/O)
SimulationResult run_simulation(
    const Config& cfg,
    const std::filesystem::path& model_path
);