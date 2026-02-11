#include "mesh_tools.h"

static std::string RankSuffix6(int rank)
{
    std::ostringstream os;
    os << "." << std::setw(6) << std::setfill('0') << rank;
    return os.str();
}

std::unique_ptr<mfem::ParMesh> CreateSimulationDomain(const std::string &path, MPI_Comm comm)
{
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

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

    // Partition to parallel mesh
    auto pmesh = std::make_unique<mfem::ParMesh>(comm, *serial);
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

  // --- 3) Apply refinement.
  refiner.Apply(pmesh);

  // STOP is satisfied if: stopping criterion met OR nothing was marked.
  // Either way, no mesh modification that we should continue from.
  if (refiner.Stop())
  {
    // Even if Stop() is true, it can be because of "no elements marked".
    // In either case, there's nothing to update/transfer.
  }
  else if (refiner.Refined())
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
    const bool did_derefine = derefiner.Apply(pmesh);
    if (did_derefine && derefiner.Derefined())
    {
      UpdateTransferAndRebalanceIfNeeded(pmesh, pfes, V);
      changed = true;
    }
  }

  return changed;
}
