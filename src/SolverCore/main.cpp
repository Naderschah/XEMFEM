#define MFEM_DEBUG
#include "mfem.hpp"
#include "load_mesh.h"
#include "boundary_conditions.h"
#include "solver.h"
#include "ComputeElectricField.h"
#include "Config.h"
#include "cmdLineParser.h"
#include "solver_api.h"
#include "parallelization.h"

// Standard Libararies 
#include <iostream>
#include <memory>
#include <cstdlib> 
#include <thread> 

// For Multithreading 
#include <omp.h>

int main(int argc, char *argv[])
{
  // Extract cmd line options 
  cli::InputParser args(argc, argv);
  if (args.has("-h") || args.has("--help")) {
      cli::print_usage(argv[0]);
      return 0;
  }

  //----------------------   Read options ---------------------------
  // config path 
  auto config_str_opt = args.get("-c");
  if (!config_str_opt) config_str_opt = args.get("--config");
  if (!config_str_opt) {
      std::cerr << "Error: missing required argument -c/--config\n";
      cli::print_usage(argv[0]);
      return 1;
  }
  auto config_path = cli::to_absolute(*config_str_opt);
  if (!std::filesystem::exists(config_path)) {
      std::cerr << "Error: config file not found: " << config_path << "\n";
      return 1;
  }
  // Load yaml config containing geometry and solver parameters
  auto cfg = std::make_shared<const Config>(
      Config::Load(config_path)
  );
  // mesh path 
  // check CLI (-m / --model)
  std::filesystem::path model_path;
    
  // First check CLI (-m / --model)
  auto model_str_opt = args.get("-m");
  if (!model_str_opt) {
      model_str_opt = args.get("--model");
  }
  if (model_str_opt)
  {
      // CLI override provided
      model_path = cli::to_absolute(*model_str_opt);

      if (!std::filesystem::exists(model_path))
      {
          std::cerr << "Error: model/mesh file not found: " << model_path << "\n";
          return 1;
      }
  }
  else
  {
      // No CLI override → fall back to config
      if (cfg->mesh.path.empty())
      {
          std::cerr << "Error: missing mesh path (config.mesh.path and -m/--model are both absent).\n";
          return 1;
      }

      model_path = cli::to_absolute(cfg->mesh.path);

      if (!std::filesystem::exists(model_path))
      {
          std::cerr << "Error: model/mesh file not found: " << model_path << "\n";
          return 1;
      }
  }

  // Log and continue 
  std::cout << "[Config] " << config_path << "\n";
  std::cout << "[Mesh]   " << model_path  << "\n";
  
  parallel::init_environment(*cfg, argc, argv)
  // Apply Voltage network
  apply_fieldcage_network(cfg);
  SimulationResult result = run_simulation(cfg, model_path);
    if (!result.success) {
        std::cerr << "Simulation failed: " << result.error_message << "\n";
        return 1;
    }
  // TODO Paths 
  // 4) Save components and magnitude
  SaveEComponents(*result.E, "field");          // field_ex.gf, field_ey.gf, (field_ez.gf)
  { std::ofstream ofs("field_mag.gf"); result.Emag->Save(ofs); }

  // 5. Save Data
  result.mesh->Save(cfg->solver.mesh_save_path.c_str());
  result.V->Save(cfg->solver.V_solution_path.c_str());

  if (!result.success) {
    std::cerr << "Simulation failed: " << result.error_message << "\n";
    return 1;
  }
  return 0;
}

