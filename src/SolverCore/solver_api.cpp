
#include "mfem.hpp"
#include "load_mesh.h"
#include "boundary_conditions.h"
#include "solver.h"
#include "ComputeElectricField.h"
#include "Config.h"
#include "cmdLineParser.h"
#include "solver_api.h"

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
    if (mesh->Dimension() != 2) { std::cerr << "Axisymmetric Simulation Geometry 3D" << std::endl;  }
    // Check r axis starts at null
    CheckAxisymmetricMesh(*mesh,
                          0,
                          cfg->solver.axisymmetric_r0_bd_attribute);
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
  InitFieldPostprocessor(*pfes, /*smooth_output=*/false);

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