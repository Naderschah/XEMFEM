
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
  BoundaryConditionGroups BCs = GetBoundaryConditionGroups(mesh.get(), cfg);

  // 4. Solve Poisson
  auto V = SolvePoisson(*pfes, BCs, cfg);

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
    // --- Create directory structure ---
    std::error_code ec;
    std::filesystem::create_directories(root_path, ec);
    if (ec) {
        std::cerr << "Warning: could not create output directory "
                  << root_path << " : " << ec.message() << "\n";
    }

    // --- Save E field TODO Add flag in config ---
    {
        std::filesystem::path EComponentPath = root_path / "E";
        SaveEComponents(*result.E, EComponentPath.string());  // E_ex.gf, E_ey.gf, ...

        std::ofstream ofs(root_path / "Emag.gf");
        if (!ofs) {
            std::cerr << "Warning: could not open Emag.gf in " << root_path << "\n";
        } else {
            result.Emag->Save(ofs);
        }
    }

    // --- Save mesh + V as serial outputs ---
    std::filesystem::path mesh_path = root_path / "simulation_mesh.msh";
    std::ofstream mesh_out(mesh_path);
    if (!mesh_out) {
        std::cerr << "Warning: could not open mesh file " << mesh_path << "\n";
    } else {
        result.mesh->PrintAsSerial(mesh_out);
    }
    std::filesystem::path V_path = root_path / "V.gf";
    result.V->SaveAsSerial(V_path.c_str(), /*precision=*/16, /*save_rank=*/0);

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

    SimulationResult result;
    result.success = false;

    // --- Explicit filenames (no prefix search) ---
    const auto mesh_path = root_path / "simulation_mesh.msh";
    const auto V_path    = root_path / "V.gf";

    if (!std::filesystem::exists(mesh_path))
    {
        result.error_message = "load_results: missing " + mesh_path.string();
        return result;
    }
    if (!std::filesystem::exists(V_path))
    {
        result.error_message = "load_results: missing " + V_path.string();
        return result;
    }

    // --- 1) Load SERIAL mesh ---
    Mesh serial_mesh;
    {
        std::ifstream mesh_in(mesh_path);
        if (!mesh_in)
        {
            result.error_message = "load_results: cannot open " + mesh_path.string();
            return result;
        }
        serial_mesh.Load(mesh_in, /*generate_edges=*/1, /*refine=*/0);
    }

    // --- 2) Build SERIAL space + load SERIAL V ---
    auto fec_h1_s  = std::make_unique<H1_FECollection>(cfg.solver.order, serial_mesh.Dimension());
    auto fes_h1_s  = std::make_unique<FiniteElementSpace>(&serial_mesh, fec_h1_s.get(), 1);

    std::unique_ptr<mfem::GridFunction> V_s;
    std::ifstream vin(V_path);
    if (!vin)
    {
        result.error_message = "load_results: cannot open " + V_path.string();
        return result;
    }
    V_s = std::make_unique<mfem::GridFunction>(&serial_mesh, vin);

    // --- 3) Convert to PAR objects required by SimulationResult ---
    result.mesh = std::make_unique<ParMesh>(MPI_COMM_WORLD, serial_mesh);

    const int dim = result.mesh->Dimension();
    result.fec  = std::make_unique<H1_FECollection>(cfg.solver.order, dim);
    result.pfes = std::make_unique<ParFiniteElementSpace>(result.mesh.get(), result.fec.get(), 1);

    result.V = std::make_unique<ParGridFunction>(result.pfes.get());
    {
        GridFunctionCoefficient Vcoeff(V_s.get());
        result.V->ProjectCoefficient(Vcoeff);
    }

    // --- 4) Rebuild E/Emag from components as you already do ---
    try
    {
        ElectricFieldPostprocessor efpp(*result.pfes);
        auto E = efpp.MakeE();
        efpp.LoadE(*E, (root_path / "E").string());
        result.E = std::move(E);

        auto Emag = efpp.MakeEmag();
        efpp.ComputeFieldMagnitude(*result.E, *Emag);
        result.Emag = std::move(Emag);
    }
    catch (const std::exception &e)
    {
        // leave E/Emag null; V + mesh valid
        std::cerr << "load_results: failed to reconstruct E/Emag: " << e.what() << "\n";
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

    if (cfg.debug.printBoundaryConditions){
        std::cout << "[DEBUG] Boundary Condition values" << std::endl;
        for (const auto& [name, b] : cfg_copy.boundaries) {
            std::cout << name << ": " << b.value << '\n';
        }
    }
    return run_name.str();
}