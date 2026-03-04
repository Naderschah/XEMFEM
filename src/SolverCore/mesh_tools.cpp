#include "mesh_tools.h"

#include <fstream>
#include <iomanip>
#include <sstream>

static std::string RankSuffix6(int rank)
{
    std::ostringstream os;
    os << "." << std::setw(6) << std::setfill('0') << rank;
    return os.str();
}

void SaveParallelMesh(const mfem::ParMesh &pmesh,
                      const std::filesystem::path &out_dir,
                      const std::string &prefix,
                      int precision)
{
    namespace fs = std::filesystem;

    const MPI_Comm comm = pmesh.GetComm();
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0)
    {
        std::error_code ec;
        fs::create_directories(out_dir, ec);
        MFEM_VERIFY(!ec,
                    "SaveParallelMesh: could not create directory '" << out_dir.string()
                    << "': " << ec.message());
    }
    MPI_Barrier(comm);

    const fs::path file_path = out_dir / (prefix + RankSuffix6(rank));
    std::ofstream os(file_path);
    MFEM_VERIFY(os.good(),
                "SaveParallelMesh: could not open '" << file_path.string() << "'");
    os.precision(precision);
    os.setf(std::ios::scientific);
    pmesh.ParPrint(os);
}

void SaveSerialMesh(const mfem::ParMesh &pmesh,
                    const std::filesystem::path &mesh_path,
                    int precision,
                    const std::string &comments)
{
    namespace fs = std::filesystem;

    const MPI_Comm comm = pmesh.GetComm();
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0)
    {
        std::error_code ec;
        const fs::path parent = mesh_path.parent_path();
        if (!parent.empty())
        {
            fs::create_directories(parent, ec);
            MFEM_VERIFY(!ec,
                        "SaveSerialMesh: could not create directory '" << parent.string()
                        << "': " << ec.message());
        }
    }
    MPI_Barrier(comm);

    std::ofstream mesh_out;
    if (rank == 0)
    {
        mesh_out.open(mesh_path);
    }

    int ok = 1;
    if (rank == 0)
    {
        ok = mesh_out.is_open() ? 1 : 0;
        if (ok == 1)
        {
            mesh_out.precision(precision);
            mesh_out.setf(std::ios::scientific);
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, comm);
    MFEM_VERIFY(ok == 1, "SaveSerialMesh: could not open serial mesh file on rank 0.");

    pmesh.PrintAsSerial(mesh_out, comments.c_str());
}

void SaveParallelVTU(mfem::ParMesh &pmesh,
                     const std::filesystem::path &out_prefix,
                     mfem::VTKFormat format,
                     bool high_order_output,
                     int compression_level,
                     bool bdr_elements)
{
    namespace fs = std::filesystem;
    const MPI_Comm comm = pmesh.GetComm();
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0)
    {
        std::error_code ec;
        const fs::path parent = out_prefix.parent_path();
        if (!parent.empty())
        {
            fs::create_directories(parent, ec);
        }
        MFEM_VERIFY(!ec,
                    "SaveParallelVTU: could not create directory '" << parent.string()
                    << "': " << ec.message());
    }
    MPI_Barrier(comm);

    pmesh.PrintVTU(out_prefix.string(), format, high_order_output, compression_level, bdr_elements);
}

void SaveAMRMeshArtifacts(mfem::ParMesh &pmesh,
                          const std::filesystem::path &mesh_dir,
                          const AMRSettings &amr_io,
                          int precision)
{
    SaveParallelMesh(pmesh, mesh_dir, "amr_mesh", precision);
    if (amr_io.save_serial)
    {
        SaveSerialMesh(pmesh, mesh_dir / "amr_mesh_serial.mesh", precision, "AMR serial mesh export");
    }
    if (amr_io.save_vtu_parallel)
    {
        SaveParallelVTU(pmesh, mesh_dir / "amr_mesh_vtu", mfem::VTKFormat::BINARY, true, 0, false);
    }
}

std::unique_ptr<mfem::ParMesh> CreateSimulationDomain(const std::string &path, MPI_Comm comm)
{
    int rank = 0, nranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    namespace fs = std::filesystem;
    const fs::path p(path);

    // ---------------------------------------------------------------------
    // Case A: Parallel restart mesh (path is a directory)
    // Convention: directory contains files "amr_mesh.000000", ...
    // ---------------------------------------------------------------------
    if (fs::exists(p) && fs::is_directory(p))
    {
        const fs::path prefix = p / "amr_mesh";               // fixed name per your saver
        const std::string fname = (prefix.string() + RankSuffix6(rank));

        if (!fs::exists(fname))
        {
            std::ostringstream msg;
            msg << "CreateSimulationDomain: missing parallel mesh piece '"
                << fname << "'. (Wrong directory, wrong prefix, or MPI size mismatch?)";
            throw std::runtime_error(msg.str());
        }

        std::ifstream mesh_in(fname);
        if (!mesh_in)
        {
            throw std::runtime_error("CreateSimulationDomain: cannot open " + fname);
        }

        auto pmesh = std::make_unique<mfem::ParMesh>(comm, mesh_in);

        if (pmesh->bdr_attributes.Size() == 0)
        {
            throw std::runtime_error("CreateSimulationDomain: no boundary attributes in parallel mesh!");
        }

        return pmesh;
    }

    // ---------------------------------------------------------------------
    // Case B: Serial mesh file (path is a file)
    // ---------------------------------------------------------------------
    if (!fs::exists(p) || !fs::is_regular_file(p))
    {
        throw std::runtime_error("CreateSimulationDomain: path does not exist or is not a file/directory: " + path);
    }

    // Load serial mesh
    auto serial = std::make_unique<mfem::Mesh>(path.c_str(), 1, 1, false);

    if (serial->bdr_attributes.Size() == 0)
    {
        throw std::runtime_error("CreateSimulationDomain: no boundary attributes in serial mesh!");
    }

    // If you want NC meshes for quad/hex local AMR later:
    // serial->EnsureNCMesh();

    // Partition to parallel mesh (explicit in-memory prepartition).
    // This keeps startup behavior consistent across execution paths without
    // writing/reading intermediate mesh files.
    int *part_raw = serial->GeneratePartitioning(nranks, 1);
    if (part_raw == nullptr)
    {
        throw std::runtime_error("CreateSimulationDomain: GeneratePartitioning returned null.");
    }
    std::unique_ptr<int[]> part(part_raw);

    auto pmesh = std::make_unique<mfem::ParMesh>(comm, *serial, part.get());
    return pmesh;
}

void CheckAxisymmetricMesh(const mfem::ParMesh &mesh,
                           int radial_coord_index)
{
    MPI_Comm comm = mesh.GetComm();
    CheckAxisymmetricMesh(mesh,
                          radial_coord_index,
                          comm);
}
void CheckAxisymmetricMesh(const mfem::ParMesh &mesh,
                           int radial_coord_index,
                           MPI_Comm comm)
{
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // 1) Must be a 2D meridian mesh (local check is sufficient)
    if (mesh.Dimension() != 2)
    {
        if (rank == 0)
        {
            std::cerr << "[Axisym] ERROR: Axisymmetric mode requires a 2D (r,z) mesh.\n";
        }
        MPI_Abort(comm, 1);
    }

    // 2) Local extrema of r
    double rmin_loc =  1e300;
    double rmax_loc = -1e300;

    for (int v = 0; v < mesh.GetNV(); ++v)
    {
        const double *X = mesh.GetVertex(v);
        const double r  = X[radial_coord_index];
        rmin_loc = std::min(rmin_loc, r);
        rmax_loc = std::max(rmax_loc, r);
    }

    // 3) Global extrema
    double rmin = 0.0, rmax = 0.0;
    MPI_Allreduce(&rmin_loc, &rmin, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&rmax_loc, &rmax, 1, MPI_DOUBLE, MPI_MAX, comm);

    // 4) Global validation
    const double tol = 1e-12 * std::max(1.0, rmax - rmin);
    if (rmin < -tol)
    {
        if (rank == 0)
        {
            std::cerr << "[Axisym] ERROR: Mesh contains vertices with r < 0.\n"
                      << "         Axisymmetric meshes must lie entirely in r >= 0.\n"
                      << "         Global minimum r detected: " << rmin << "\n";
        }
        MPI_Abort(comm, 1);
    }

    // 5) Informational output (rank 0 only)
    if (rank == 0)
    {
        std::cout << "[Axisym] Mesh OK. Global r-range: ["
                  << rmin << ", " << rmax << "]\n";
    }
}

// ------------------------------ AMR Specifics ----------------------------------------------
static void UpdateTransferAndRebalanceIfNeeded(mfem::ParMesh &pmesh,
                                              mfem::ParFiniteElementSpace &pfes,
                                              mfem::ParGridFunction &V)
{
  // 1) Build/update transfer operator for the current mesh change (refine/derefine/rebalance).
  pfes.Update();
  V.Update();

  // 2) In parallel nonconforming AMR, rebalance to avoid rank imbalance.
  if (pmesh.Nonconforming())
  {
    pmesh.Rebalance();
    pfes.Update();
    V.Update();
  }

  // 3) Free update matrices to save memory (must be called after all GF updates).
  pfes.UpdatesFinished();
}

bool ApplyAMRRefineDerefineStep(mfem::ParMesh &pmesh,
                                mfem::ParFiniteElementSpace &pfes,
                                mfem::ParGridFunction &V,
                                const Config &cfg)
{
  if (!cfg.mesh.amr.enable) { return false; }
  const bool amr_verbose = cfg.mesh.amr.verbose;
  const MPI_Comm comm = pmesh.GetComm();
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  auto global_counts = [&](long long &ne_global, long long &tdof_global)
  {
    const long long ne_local = static_cast<long long>(pmesh.GetNE());
    const long long tdof_local = static_cast<long long>(pfes.GetTrueVSize());
    MPI_Allreduce(&ne_local, &ne_global, 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&tdof_local, &tdof_global, 1, MPI_LONG_LONG, MPI_SUM, comm);
  };

  long long ne_before = 0, tdof_before = 0;
  global_counts(ne_before, tdof_before);

  // --- 1) Build an element-wise error estimator (ex15p-style).
  //
  // KellyErrorEstimator requires an integrator with ComputeElementFlux().
  // For Poisson with (Q grad u, grad v), DiffusionIntegrator is appropriate.
  // If you have a spatially varying coefficient, you can wire it here instead of "one".
  mfem::ConstantCoefficient one(1.0);
  mfem::DiffusionIntegrator integ(one);

  const int dim  = pmesh.Dimension();
  const int sdim = pmesh.SpaceDimension();

  // ex15p uses an L2 space for discontinuous flux storage.
  mfem::L2_FECollection flux_fec(cfg.solver.order, dim);

  // IMPORTANT: In ex15p they allocate flux_fes with new and give it to the estimator.
  // MFEM estimators take ownership of the flux_fes pointer they are constructed with.
  auto *flux_fes = new mfem::ParFiniteElementSpace(&pmesh, &flux_fec, sdim);

  std::unique_ptr<mfem::ErrorEstimator> estimator =
      std::make_unique<mfem::KellyErrorEstimator>(integ, V, flux_fes);

  // --- 2) Configure refiner (purely local threshold driven by LocalErrorGoal).
  mfem::ThresholdRefiner refiner(*estimator);
  refiner.SetTotalErrorFraction(0.0);                 // purely local threshold
  refiner.SetLocalErrorGoal(cfg.mesh.amr.local_error_goal);
  refiner.SetNCLimit(cfg.mesh.amr.nc_limit);

  if (cfg.mesh.amr.max_elements > 0) { refiner.SetMaxElements(cfg.mesh.amr.max_elements); }

  if (cfg.mesh.amr.prefer_conforming_refinement)
  {
    refiner.PreferConformingRefinement();
  }
  else
  {
    refiner.PreferNonconformingRefinement();
  }

  // Force recomputation of errors for this call (matches ex15p pattern).
  refiner.Reset();

  bool changed = false;
  bool refined = false;
  bool stop_after_refine = false;
  bool did_derefine = false;
  bool derefined = false;

  // --- 3) Apply refinement.
  refiner.Apply(pmesh);
  stop_after_refine = refiner.Stop();
  refined = refiner.Refined();

  // STOP is satisfied if: stopping criterion met OR nothing was marked.
  // Either way, no mesh modification that we should continue from.
  if (stop_after_refine)
  {
    // Even if Stop() is true, it can be because of "no elements marked".
    // In either case, there's nothing to update/transfer.
  }
  else if (refined)
  {
    UpdateTransferAndRebalanceIfNeeded(pmesh, pfes, V);
    changed = true;
  }

  // --- 4) Optional derefinement (coarsening) using hysteresis*goal.
  if (cfg.mesh.amr.enable_derefine)
  {
    mfem::ThresholdDerefiner derefiner(*estimator);
    derefiner.SetNCLimit(cfg.mesh.amr.nc_limit);

    const double hysteresis = (cfg.mesh.amr.derefine_hysteresis > 0.0)
                                ? cfg.mesh.amr.derefine_hysteresis
                                : 0.5; // safe default
    derefiner.SetThreshold(hysteresis * cfg.mesh.amr.local_error_goal);

    derefiner.Reset();

    // MeshOperator::Apply returns true if something happened (derefined + continue).
    did_derefine = derefiner.Apply(pmesh);
    derefined = derefiner.Derefined();
    if (did_derefine && derefined)
    {
      UpdateTransferAndRebalanceIfNeeded(pmesh, pfes, V);
      changed = true;
    }
  }

  if (amr_verbose)
  {
    long long ne_after = 0, tdof_after = 0;
    global_counts(ne_after, tdof_after);
    if (rank == 0)
    {
      std::cout << "[AMR] refine: stop=" << (stop_after_refine ? "true" : "false")
                << ", refined=" << (refined ? "true" : "false")
                << " | derefine: attempted=" << (cfg.mesh.amr.enable_derefine ? "true" : "false")
                << ", applied=" << (did_derefine ? "true" : "false")
                << ", derefined=" << (derefined ? "true" : "false")
                << " | changed=" << (changed ? "true" : "false")
                << "\n";
      std::cout << "[AMR] global size: NE " << ne_before << " -> " << ne_after
                << ", true_dofs " << tdof_before << " -> " << tdof_after << "\n";
    }
  }

  return changed;
}
