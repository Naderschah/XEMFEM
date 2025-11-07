// src/config/Config.cpp
#include "Config.h"
#include <yaml-cpp/yaml.h>
#include <stdexcept>
#include <iostream>

// read int or "auto"
static std::pair<int,bool> read_int_or_auto(const YAML::Node& n, int dflt_val, bool dflt_auto=true)
{
    if (!n) return {dflt_val, dflt_auto};
    if (n.IsScalar()) {
        const std::string s = n.as<std::string>();
        if (s == "auto" || s == "AUTO" || s == "Auto") return {dflt_val, true};
        try {
            int v = std::stoi(s);
            return {v, false};
        } catch (...) {
            // If scalar but not parseable, fall back to defaults
            return {dflt_val, dflt_auto};
        }
    }
    if (n.IsSequence() || n.IsMap()) {
        // unexpected type; keep defaults
        return {dflt_val, dflt_auto};
    }
    // Try as int directly
    try { return { n.as<int>(), false }; }
    catch (...) { return { dflt_val, dflt_auto }; }
}


// apply preset defaults first; explicit fields will override later
static void apply_preset_defaults(Config &cfg)
{
  const std::string p = cfg.compute.preset;

  if (p.empty()) return;
  if (p=="null") return; 

  // start from conservative defaults
  cfg.compute.mpi.enabled = false;
  cfg.compute.mpi.ranks_auto = true;
  cfg.compute.mpi.ranks = 1;
  cfg.compute.threads.enabled = true;
  cfg.compute.threads.num_auto = true;
  cfg.compute.threads.num = 0;
  cfg.compute.threads.affinity = "compact";
  cfg.compute.device.type = "none";
  cfg.compute.device.id_auto = true;
  cfg.compute.device.id = 0;
  cfg.compute.device.per_rank = 1;

  if (p == "serial") {
    cfg.compute.mpi.enabled = false;
    cfg.compute.threads.enabled = false;
    cfg.compute.device.type = "none";
    cfg.solver.assembly_mode = "full";
  }
  else if (p == "threads") {
    cfg.compute.mpi.enabled        = false;
    cfg.compute.threads.enabled    = true;
    cfg.compute.threads.num        = -1;          // auto → we'll cap to available CPUs
    cfg.compute.threads.affinity   = "compact";   // compact | scatter | none
    cfg.compute.device.type        = "none";
    cfg.solver.assembly_mode       = "partial";
  }
  else if (p == "mpi") {
    cfg.compute.mpi.enabled        = true;
    cfg.compute.mpi.ranks_auto     = true;
    cfg.compute.threads.enabled    = false;
    cfg.compute.device.type        = "none";
    cfg.solver.assembly_mode       = "full";
  }
  else if (p == "gpu") { 
    cfg.compute.mpi.enabled = false;
    cfg.compute.threads.enabled = true;
    cfg.compute.device.type = "cuda"; // default; user can override to hip/occa/cpu
    if (cfg.solver.assembly_mode.empty()) cfg.solver.assembly_mode = "partial";
  }
  else if (p == "mpi+gpu") {
    cfg.compute.mpi.enabled = true;
    cfg.compute.threads.enabled = false;
    cfg.compute.device.type = "cuda";
    cfg.compute.device.per_rank = 1;
    if (cfg.solver.assembly_mode.empty()) cfg.solver.assembly_mode = "partial";
  }
  else if (p == "mpi+threads") {
    cfg.compute.mpi.enabled       = true;
    cfg.compute.threads.enabled   = true;
    cfg.compute.threads.num       = -1;
    cfg.compute.threads.affinity  = "scatter";
    cfg.compute.device.type       = "none";
    cfg.solver.assembly_mode      = "partial";
  }
  else {
    std::cout << "\033[0;33mWARNING Preset " << p << " not recognized\033[0m" << std::endl;
  }
}



// parse compute.* blocks
static void parse_compute(Config &cfg, const YAML::Node& root)
{
    // preset
    if (root["preset"])
        cfg.compute.preset = root["preset"].as<std::string>("");

    // Apply preset defaults first (explicit fields below can override)
    apply_preset_defaults(cfg);

    if (root["compute"]) {
        const auto C = root["compute"];

        // mpi
        if (C["mpi"]) {
            const auto M = C["mpi"];
            cfg.compute.mpi.enabled = M["enabled"].as<bool>(cfg.compute.mpi.enabled);

            // ranks: int or "auto"
            auto [rk, rk_auto] = read_int_or_auto(M["ranks"], cfg.compute.mpi.ranks, cfg.compute.mpi.ranks_auto);
            cfg.compute.mpi.ranks = rk;
            cfg.compute.mpi.ranks_auto = rk_auto;

            cfg.compute.mpi.hostfile = M["hostfile"].as<std::string>(cfg.compute.mpi.hostfile);
            cfg.compute.mpi.repartition_after_refine =
                M["repartition_after_refine"].as<bool>(cfg.compute.mpi.repartition_after_refine);
        }

        // threads
        if (C["threads"]) {
            const auto T = C["threads"];
            cfg.compute.threads.enabled = T["enabled"].as<bool>(cfg.compute.threads.enabled);
            auto [tnum, t_auto] = read_int_or_auto(T["num"], cfg.compute.threads.num, cfg.compute.threads.num_auto);
            cfg.compute.threads.num = tnum;
            cfg.compute.threads.num_auto = t_auto;
            cfg.compute.threads.affinity = T["affinity"].as<std::string>(cfg.compute.threads.affinity);
        }

        // device
        if (C["device"]) {
            const auto D = C["device"];
            cfg.compute.device.type = D["type"].as<std::string>(cfg.compute.device.type);
            auto [did, d_auto] = read_int_or_auto(D["id"], cfg.compute.device.id, cfg.compute.device.id_auto);
            cfg.compute.device.id = did;
            cfg.compute.device.id_auto = d_auto;
            cfg.compute.device.per_rank = D["per_rank"].as<int>(cfg.compute.device.per_rank);
        }
    }

    // Convenience: if GPU selected and no assembly_mode explicitly set → partial
    if (cfg.compute.device.type != "none" && cfg.solver.assembly_mode.empty())
        cfg.solver.assembly_mode = "partial";

    // If threads disabled → force single-thread behavior
    if (!cfg.compute.threads.enabled) {
        cfg.compute.threads.num_auto = false;
        if (cfg.compute.threads.num <= 0) cfg.compute.threads.num = 1;
    }
}



// -----------------------------------------------------------------------------
// Internal common parser
// -----------------------------------------------------------------------------
namespace {
Config LoadFromNode(const YAML::Node &root) {
  Config cfg;

  cfg.schema_version = root["schema_version"].as<int>(1);
  cfg.geometry_id    = root["geometry_id"].as<std::string>("");

  // --- Mesh path
  if (root["mesh"] && root["mesh"]["path"])
    cfg.mesh.path = root["mesh"]["path"].as<std::string>("geometry.msh");


  // --- Debug
  if (root["debug"]) {
    cfg.debug.debug      = root["debug"]["debug"].as<bool>(false);
    cfg.debug.quick_mesh = root["debug"]["quick_mesh"].as<bool>(false);
  }

  // --- Solver
  if (root["solver"]) {
    const auto s = root["solver"];
    cfg.solver.atol        = s["atol"].as<double>(1.0);
    cfg.solver.rtol        = s["rtol"].as<double>(0.0);
    cfg.solver.maxiter     = s["maxiter"].as<int>(100000);
    cfg.solver.printlevel  = s["printlevel"].as<int>(1);

    cfg.solver.mesh_save_path     = s["mesh_save_path"].as<std::string>("simulation_mesh.msh");
    cfg.solver.V_solution_path    = s["V_solution_path"].as<std::string>("solution_V.gf");
    cfg.solver.Emag_solution_path = s["Emag_solution_path"].as<std::string>("solution_Emag.gf");

    // NEW: assembly/solver/precond (optional)
    if (s["assembly_mode"]) cfg.solver.assembly_mode = s["assembly_mode"].as<std::string>(cfg.solver.assembly_mode);
    if (s["solver"])        cfg.solver.solver        = s["solver"].as<std::string>(cfg.solver.solver);
    if (s["precond"])       cfg.solver.precond       = s["precond"].as<std::string>(cfg.solver.precond);
  }

  // --- Materials
  if (root["materials"]) {
    for (const auto &it : root["materials"]) {
      std::string name = it.first.as<std::string>();
      auto node = it.second;
      Material m;
      m.id         = node["attr_id"]   ? node["attr_id"].as<int>(-1) : -1;
      m.epsilon_r  = node["epsilon_r"] ? node["epsilon_r"].as<double>(1.0) : 1.0;
      cfg.materials[name] = m;
    }
  }

  // --- Boundaries
  if (root["boundaries"]) {
    for (const auto &it : root["boundaries"]) {
      std::string name = it.first.as<std::string>();
      auto node = it.second;
      Boundary b;
      b.bdr_id = node["bdr_id"].as<int>(-1);
      b.type   = node["type"].as<std::string>("dirichlet");
      b.value  = node["value"].as<double>(0.0);
      if (b.bdr_id <= 0)
          throw std::runtime_error("Boundary '" + name + "' is missing a valid bdr_id");
      cfg.boundaries[name] = b;
    }
  }

  // --- Compute (preset + mpi/threads/device) with precedence & back-compat
  parse_compute(cfg, root);


  return cfg;
}
} // namespace
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Public loaders
// -----------------------------------------------------------------------------
Config Config::Load(const std::string &path) {
    YAML::Node root = YAML::LoadFile(path);
    return LoadFromNode(root);
}

Config Config::LoadFromString(const std::string &yaml_str) {
    YAML::Node root = YAML::Load(yaml_str);
    return LoadFromNode(root);
}