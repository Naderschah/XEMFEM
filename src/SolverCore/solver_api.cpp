
#include "mfem.hpp"
#include "load_mesh.h"
#include "boundary_conditions.h"
#include "solver.h"
#include "ComputeElectricField.h"
#include "Config.h"
#include "cmdLineParser.h"
#include "solver_api.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <omp.h>
#include <filesystem>



using namespace mfem;


SimulationResult run_simulation(std::shared_ptr<Config> cfg,
                                const std::filesystem::path& model_path)
{
  // What will this default to if I do nothing with MPI
  MPI_Comm comm = MPI_COMM_WORLD;

  // Solve voltages in circuit
  apply_fieldcage_network(*cfg);

  // 1. Create the mesh
  auto mesh = CreateSimulationDomain(model_path, comm); 
  if (cfg->solver.axisymmetric) {
    if (mesh->Dimension() != 2) { std::cerr << "Axisymmetric Simulation but Geometry is not 2D" << std::endl;  }
    // Check r axis starts at null
    CheckAxisymmetricMesh(*mesh,0);
  }

  // 2. Create finite element collection and space
  std::unique_ptr<mfem::FiniteElementCollection> fec =
    std::make_unique<mfem::H1_FECollection>(cfg->solver.order,
                                            mesh->Dimension());
  auto pfes = std::make_unique<mfem::ParFiniteElementSpace>(mesh.get(), fec.get());

  // 3. Get Dirichlet boundary attributes
  Array<int> dirichlet_arr = GetDirichletAttributes(mesh.get(), cfg);
  // FIXME: For Neumann and Robin and axisymmetric we still need to supply boundary markers


  // 4. Solve Poisson
  auto V = SolvePoisson(*pfes, dirichlet_arr, cfg);

  std::unique_ptr<mfem::FiniteElementSpace> vec_fes, scalar_fes;
  mfem::ParMesh *pmesh = mesh.get();  // mesh is std::unique_ptr<ParMesh>

  vec_fes    = std::make_unique<mfem::ParFiniteElementSpace>(pmesh,
                                                            fec.get(),
                                                            pmesh->Dimension());
  scalar_fes = std::make_unique<mfem::ParFiniteElementSpace>(pmesh,
                                                            fec.get(),
                                                            1);

  // Initialize postprocessor 
  InitFieldPostprocessor(*pfes);

  // 3) Allocate output fields on the right spaces
  std::unique_ptr<mfem::GridFunction> E    = CreateE();  
  std::unique_ptr<mfem::GridFunction> Emag = CreateEmag();  

  // 4) Compute E and |E|
  ComputeElectricField(*V, *E, /*scale=*/-1.0); // V/m
  ComputeFieldMagnitude(*E, *Emag);

  // 5) Pack result and transfer ownership of all dependent objects
  SimulationResult result;
  result.mesh = std::move(mesh);
  result.fec  = std::move(fec);
  result.pfes = std::move(pfes);
  result.V    = std::move(V);
  result.E    = std::move(E);
  result.Emag = std::move(Emag);

  result.success = true;

  return result;
}


void save_results(const SimulationResult &result, const std::filesystem::path &root_path)
{
  // TODO Make some config entry for which to save

  {
  // Save E field components 
  std::filesystem::path EComponentPath = root_path / "E";
  SaveEComponents(*result.E, EComponentPath.string());  // E_ex.gf, E_ey.gf, ...
  std::ofstream ofs(root_path / "Emag.gf");
  if (!ofs) {
      std::cerr << "Warning: could not open Emag.gf in " << root_path << "\n";
  } else {
      result.Emag->Save(ofs);
  }
  // Save V and mesh 
  std::filesystem::path mesh_path = root_path / "simulation_mesh.msh";
  result.mesh->Save(mesh_path.c_str());
  std::filesystem::path V_path = root_path / "V.gf";
  result.V->Save(V_path.c_str());
  }

  {
  // ---------- ParaView output ----------
  // Collection name appears in the .pvd file
  std::string collection_name = "Simulation";
  auto *pvdc = new ParaViewDataCollection(collection_name, result.mesh.get());

  // Write under root_path (directory will be created if needed)
  pvdc->SetPrefixPath(root_path.string());

  // Register fields (they can be H1, L2, etc., as long as they live on result.mesh)
  pvdc->RegisterField("V", result.V.get());       // H1 scalar
  pvdc->RegisterField("E", result.E.get());       // L2 vector (components in the GF)
  pvdc->RegisterField("Emag", result.Emag.get()); // scalar magnitude

  // Optional: control output quality
  int order = result.V->FESpace()->GetOrder(0);
  pvdc->SetLevelsOfDetail(order);
  pvdc->SetDataFormat(VTKFormat::BINARY);
  pvdc->SetHighOrderOutput(true);

  pvdc->Save();
  }

}
namespace
{
std::optional<std::filesystem::path>
find_first_with_prefix(const std::filesystem::path &dir,
                       const std::string &prefix)
{
    std::error_code ec;
    if (!std::filesystem::exists(dir, ec)) {
        return std::nullopt;
    }

    for (const auto &entry : std::filesystem::directory_iterator(dir, ec)) {
        if (ec) break;
        if (!entry.is_regular_file()) continue;

        const std::string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) == 0) {  // starts with prefix
            return entry.path();
        }
    }
    return std::nullopt;
}
} // namespace

SimulationResult load_results(const Config &cfg,
                              const std::filesystem::path &root_path)
{
    using namespace mfem;
    std::cerr << "[LOAD RESULT] This will error later something about the way we create the mesh TODO Fix\n";

    SimulationResult result;
    result.success = false;
    result.error_message.clear();

    std::error_code ec;
    if (!std::filesystem::exists(root_path, ec)) {
        result.error_message =
            "load_results: directory does not exist: " + root_path.string();
        std::cerr << result.error_message << "\n";
        return result;
    }

    // ----------------- 1) Load serial mesh -----------------
    auto mesh_file = find_first_with_prefix(root_path, "simulation_mesh.msh");
    if (!mesh_file) {
        result.error_message =
            "load_results: mesh file with prefix 'simulation_mesh.msh' "
            "not found in " + root_path.string();
        std::cerr << result.error_message << "\n";
        return result;
    }

    std::ifstream mesh_in(*mesh_file);
    if (!mesh_in) {
        result.error_message =
            "load_results: cannot open mesh file: " + mesh_file->string();
        std::cerr << result.error_message << "\n";
        return result;
    }

    Mesh serial_mesh;
    try {
        serial_mesh.Load(mesh_in, 1, 1);
    } catch (const std::exception &e) {
        result.error_message =
            std::string("load_results: failed to load serial Mesh: ") + e.what();
        std::cerr << result.error_message << "\n";
        return result;
    }

    // ----------------- 2) Build ParMesh from serial Mesh -----------------
    try {
        result.mesh = std::make_unique<ParMesh>(MPI_COMM_WORLD, serial_mesh);
    } catch (const std::exception &e) {
        result.error_message =
            std::string("load_results: failed to build ParMesh: ") + e.what();
        std::cerr << result.error_message << "\n";
        return result;
    }

    const int dim = result.mesh->Dimension();

    // ----------------- 3) Rebuild H1 FE space for V -----------------
    // Use the same order you used in the original run.
    const int order_h1 = cfg.solver.order;  // adjust field name if different

    auto fec_h1  = std::make_unique<H1_FECollection>(order_h1, dim);
    auto pfes_h1 = std::make_unique<ParFiniteElementSpace>(
        result.mesh.get(), fec_h1.get(), /*vdim=*/1);

    // ----------------- 4) Load V.gf into ParGridFunction -----------------
    auto V_file = find_first_with_prefix(root_path, "V.gf");
    if (!V_file) {
        result.error_message =
            "load_results: file with prefix 'V.gf' not found in " +
            root_path.string();
        std::cerr << result.error_message << "\n";
        return result;
    }

    {
        std::ifstream in(*V_file);
        if (!in) {
            result.error_message =
                "load_results: cannot open V.gf file: " + V_file->string();
            std::cerr << result.error_message << "\n";
            return result;
        }

        auto V = std::make_unique<ParGridFunction>(pfes_h1.get());
        V->Load(in);
        result.V = std::move(V);
    }

    // store FE data if SimulationResult has these members
    result.fec  = std::move(fec_h1);
    result.pfes = std::move(pfes_h1);

    // ----------------- 5) Rebuild E and Emag via ElectricFieldPostprocessor ----
    try {
        // ElectricFieldPostprocessor takes a FiniteElementSpace&, and
        // ParFiniteElementSpace derives from FiniteElementSpace.
        ElectricFieldPostprocessor efpp(*result.pfes);

        // E: vector L2 field reconstructed from saved components
        auto E = efpp.MakeE();  // std::unique_ptr<mfem::GridFunction>
        std::string prefix = (root_path / "E").string(); // same as SaveComponents prefix
        efpp.LoadE(*E, prefix);
        result.E = std::move(E);

        // Emag: recompute from E
        auto Emag = efpp.MakeEmag();  // std::unique_ptr<mfem::GridFunction>
        efpp.ComputeFieldMagnitude(*result.E, *Emag);
        result.Emag = std::move(Emag);
    } catch (const std::exception &e) {
        std::cerr << "load_results: failed to reconstruct E/Emag from components: "
                  << e.what() << "\n";
        // leave E / Emag null; V + mesh are still valid
    }

    result.success = true;
    result.error_message.clear();
    return result;
}


// Single Simulation with file path handling 

std::string run_one(const Config &cfg,
                           const std::vector<std::pair<std::string, std::string>> &active_params,
                           std::size_t run_index) // 0-based
{
    std::filesystem::path model_path = cfg.mesh.path;
    // Determine root save directory from Config
    std::filesystem::path save_root(cfg.save_path);
    std::error_code ec;

    // Check save path is empty 
    if (std::filesystem::exists(save_root)) {
        bool empty = std::filesystem::is_empty(save_root, ec);
        if (ec) {
            throw std::runtime_error("Failed to check save_path directory: " + ec.message());
        }
        if (!empty) {
            if (!cfg.delete_files_present) {
                throw std::runtime_error(
                    "Directory '" + save_root.string() + "' is not empty.");
            } else {
                // Delete all contents but not the directory itself
                for (const auto& entry : std::filesystem::directory_iterator(save_root)) {
                    std::filesystem::remove_all(entry.path(), ec);
                    if (ec) {
                        throw std::runtime_error(
                            "Failed to delete '" + entry.path().string() +
                            "': " + ec.message());
                    }
                }
            }
        }
    }

    // Ensure directory exists
    std::filesystem::create_directories(save_root, ec);
    if (ec)
    {
        std::cerr << "Warning: could not create save_root directory "
                  << save_root << " : " << ec.message() << "\n";
    }

    // Per-run directory: run_0001, run_0002, ...
    std::size_t display_index = run_index + 1;  // make it 1-based for naming
    std::ostringstream run_name;
    run_name << "run_" << std::setw(4) << std::setfill('0') << display_index;

    std::filesystem::path run_dir = save_root / run_name.str();
    std::filesystem::create_directories(run_dir, ec);
    if (ec)
    {
        std::cerr << "Warning: could not create run directory "
                  << run_dir << " : " << ec.message() << "\n";
    }

    // Copy config and override solver output paths so they land in run_dir
    Config cfg_copy = cfg;

    auto cfg_ptr = std::make_shared<Config>(cfg_copy);
    if (!cfg.debug.dry_run)
    {
        SimulationResult result = run_simulation(cfg_ptr, model_path);
        if (!result.success)
        {
            std::cerr << "Simulation failed for " << run_name.str()
                      << ": " << result.error_message << "\n";
            return {};  // signal failure to caller
        }

        save_results(result, save_root);
    }

    if (cfg.debug.debug){
        std::cout << "[DEBUG] Boundary Condition values" << std::endl;
        for (const auto& [name, b] : cfg_copy.boundaries) {
            std::cout << name << ": " << b.value << '\n';
        }
    }
    return run_name.str();
}