#include "interpolator.h"

// Lazy TODO move to utils module same function in optimzation.cpp
static inline std::string rtrim(const std::string &s)
{
    std::string::size_type end = s.find_last_not_of(" \t\r\n");
    if (end == std::string::npos) return "";
    return s.substr(0, end + 1);
}
static std::vector<std::string> read_root_level_runs_from_meta(const std::filesystem::path &save_root)
{
    std::vector<std::string> runs;

    std::filesystem::path meta_path = save_root / "meta.txt";
    if (!std::filesystem::exists(meta_path)) {
        return runs; // empty → no meta.txt
    }

    std::ifstream meta(meta_path);
    if (!meta) {
        std::cerr << "Warning: could not open meta.txt in " << save_root << "\n";
        return runs;
    }

    std::string line;
    while (std::getline(meta, line)) {
        // Ignore empty lines
        if (line.empty()) continue;

        // Ignore lines starting with whitespace (indented blocks)
        if (!line.empty() && (line[0] == ' ' || line[0] == '\t')) {
            continue;
        }

        // We only care about root-level lines like "run_0001:"
        std::string trimmed = rtrim(line);
        if (trimmed.size() < 2) continue;

        if (trimmed.back() == ':') {
            std::string name = trimmed.substr(0, trimmed.size() - 1);
            if (!name.empty()) {
                runs.push_back(name);
            }
        }
    }

    return runs;
}

// ------------------------- H1 Projection --------------------------
struct ProjectedH1VectorField
{
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    std::unique_ptr<mfem::ParFiniteElementSpace>   fes;
    std::unique_ptr<mfem::ParGridFunction>         E;
};

static ProjectedH1VectorField BuildEFieldH1Projection(mfem::ParMesh &pmesh,
                                                      const mfem::ParGridFunction &V,
                                                      Config cfg)
{
    // Solve minimization of grad V - E = 0 over FESpace where E is H1 and gradV is L2
    using namespace mfem;

    const int field_dim = pmesh.Dimension();
    const ParFiniteElementSpace *fesV = V.ParFESpace();
    const int order = fesV->GetMaxElementOrder();

    ProjectedH1VectorField out;
    out.fec = std::make_unique<H1_FECollection>(order, field_dim);
    out.fes = std::make_unique<ParFiniteElementSpace>(&pmesh, out.fec.get(), field_dim, Ordering::byVDIM);
    out.E   = std::make_unique<ParGridFunction>(out.fes.get());
    *(out.E) = 0.0;

    ParBilinearForm m(out.fes.get());
    m.AddDomainIntegrator(new VectorMassIntegrator());
    m.Assemble();
    m.Finalize();

    ParLinearForm b(out.fes.get());
    GradientGridFunctionCoefficient gradV(&V); // +grad(V)
    b.AddDomainIntegrator(new VectorDomainLFIntegrator(gradV));
    b.Assemble();

    std::unique_ptr<mfem::HypreParMatrix> M(m.ParallelAssemble());
    std::unique_ptr<mfem::HypreParVector> B(b.ParallelAssemble());
    mfem::HypreParVector X(*B);  // same partitioning/layout, owns its memory
    X = 0.0;

    mfem::HypreSmoother prec;
    prec.SetType(mfem::HypreSmoother::l1GS);   // good default on CPU for mass matrices
    prec.SetOperator(*M);

    mfem::CGSolver cg(pmesh.GetComm());
    cg.SetPrintLevel(0);
    cg.SetRelTol(cfg.solver.rtol);
    cg.SetAbsTol(cfg.solver.atol);
    cg.SetMaxIter(cfg.solver.maxiter);
    cg.SetPreconditioner(prec);
    cg.SetOperator(*M);
    cg.Mult(*B, X);

    out.E->Distribute(X);

    // Convert +grad(V) -> -grad(V)
    *(out.E) *= -1.0;

    return out;
}
// ------------------------- End of H1 Projection ---------------------------
// --------------------------- Sampling Helpers ---------------------------
static mfem::Vector BuildPointPositionsByNodes(const int space_dim,
                                               const int Nx,
                                               const int Ny,
                                               const int Nz,
                                               const mfem::Vector &bbmin3,
                                               const mfem::Vector &spacing)
{
    using namespace mfem;
    const int N = Nx * Ny * Nz;
    Vector point_pos(space_dim * N);
    int p = 0;
    for (int k = 0; k < Nz; ++k)
    for (int j = 0; j < Ny; ++j)
    for (int i = 0; i < Nx; ++i)
    {
        point_pos[0 * N + p] = bbmin3[0] + i * spacing[0];
        point_pos[1 * N + p] = bbmin3[1] + j * spacing[1];
        if (space_dim == 3)
        {
            const double z = bbmin3[2] + k * spacing[2];
            point_pos[2 * N + p] = z;
        }
        ++p;
    }
    return point_pos;
}

static void SampleGradV(mfem::ParMesh &pmesh,
                        const mfem::ParGridFunction &V,
                        const mfem::Array<unsigned int> &code,
                        const mfem::Array<unsigned int> &elid,
                        const mfem::Vector &ref_by_nodes, // component-major: ref[d*N + p]
                        const int N,
                        const bool accept_surface_projection,
                        GridSample &out)
{
    using namespace mfem;

    const int field_dim = pmesh.Dimension();

    GradientGridFunctionCoefficient gradV(&V);
    Vector g(field_dim);
    IntegrationPoint ip;

    for (int p = 0; p < N; ++p)
    {
        const unsigned int c = code[p];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf)) { continue; }

        const double r0 = ref_by_nodes[0 * N + p];
        const double r1 = ref_by_nodes[1 * N + p];
        const double r2 = (field_dim == 3) ? ref_by_nodes[2 * N + p] : 0.0;

        if (field_dim == 2) { ip.Set2(r0, r1); }
        else                { ip.Set3(r0, r1, r2); }

        const int e = static_cast<int>(elid[p]);
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        T->SetIntPoint(&ip);

        gradV.Eval(g, *T, ip);
        g *= -1.0;

        out.valid[p] = 1;
        out.Ex[p] = g[0];
        out.Ey[p] = g[1];
        if (field_dim == 3) { out.Ez[p] = g[2]; }
        out.Emag[p] = std::sqrt(g * g);
    }
}

static void SampleH1Projection(mfem::ParMesh &pmesh,
                              const mfem::ParGridFunction &E_h, // vdim == field_dim
                              const mfem::Array<unsigned int> &code,
                              const mfem::Array<unsigned int> &elid,
                              const mfem::Vector &ref_by_nodes, // component-major: ref[d*N + p]
                              const int N,
                              const bool accept_surface_projection,
                              GridSample &out)
{
    using namespace mfem;

    const int field_dim = pmesh.Dimension();
    Vector g(field_dim);
    IntegrationPoint ip;

    for (int p = 0; p < N; ++p)
    {
        const unsigned int c = code[p];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf)) { continue; }

        const double r0 = ref_by_nodes[0 * N + p];
        const double r1 = ref_by_nodes[1 * N + p];
        const double r2 = (field_dim == 3) ? ref_by_nodes[2 * N + p] : 0.0;

        if (field_dim == 2) { ip.Set2(r0, r1); }
        else                { ip.Set3(r0, r1, r2); }

        const int e = static_cast<int>(elid[p]);
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        T->SetIntPoint(&ip);

        // Evaluate the projected continuous vector field at this point
        E_h.GetVectorValue(*T, ip, g);

        out.valid[p] = 1;
        out.Ex[p] = g[0];
        out.Ey[p] = g[1];
        if (field_dim == 3) { out.Ez[p] = g[2]; }
        out.Emag[p] = std::sqrt(g * g);
    }
}


void SampleEFieldOnCartesianGrid(mfem::ParMesh &pmesh,
                                 const mfem::ParGridFunction &V,
                                 const int Nx, const int Ny, const int Nz,
                                 GridSample &out,
                                 Config cfg,
                                 const bool H1_project,
                                 const bool accept_surface_projection)
{
    using namespace mfem;

    const int field_dim = pmesh.Dimension();      // 2 or 3: size of grad(V)
    const int space_dim = pmesh.SpaceDimension(); // 2 or 3: physical coordinate dimension for FindPoints
    MFEM_VERIFY(field_dim == 2 || field_dim == 3, "Only 2D/3D elements supported");
    MFEM_VERIFY(space_dim == 2 || space_dim == 3, "Only 2D/3D physical space supported");

    // Bounding box in physical space, then pad to 3D for export
    Vector bbmin_s(space_dim), bbmax_s(space_dim);
    pmesh.GetBoundingBox(bbmin_s, bbmax_s);

    mfem::Vector bbmin3(3), bbmax3(3);
    bbmin3 = 0.0; bbmax3 = 0.0;
    for (int d = 0; d < space_dim; ++d) { bbmin3[d] = bbmin_s[d]; bbmax3[d] = bbmax_s[d]; }

    // Output grid metadata: ALWAYS 3D coordinates (x,y,z)
    out.Nx = Nx; out.Ny = Ny; out.Nz = Nz;

    out.origin.SetSize(3);
    out.spacing.SetSize(3);
    out.origin = bbmin3;
    out.spacing = 0.0;

    MFEM_VERIFY(Nx > 1 && Ny > 1, "Nx and Ny must be > 1");
    out.spacing[0] = (bbmax3[0] - bbmin3[0]) / (Nx - 1);
    out.spacing[1] = (bbmax3[1] - bbmin3[1]) / (Ny - 1);

    // For space_dim==3, z spacing comes from bbox; for space_dim==2 it's a dummy axis.
    if (Nz > 1 && space_dim == 3)
    {
        out.spacing[2] = (bbmax3[2] - bbmin3[2]) / (Nz - 1);
    }
    else
    {
        out.spacing[2] = 0.0; // single slice OR 2D physical space => repeated slices
    }

    const int N = Nx * Ny * Nz;

    out.Ex.assign(N, std::numeric_limits<double>::quiet_NaN());
    out.Ey.assign(N, std::numeric_limits<double>::quiet_NaN());
    out.Ez.assign(N, 0.0); // kept even for 2D field_dim; writer can ignore/omit based on availability
    out.Emag.assign(N, std::numeric_limits<double>::quiet_NaN());
    out.valid.assign(N, 0);

    // Ensure Nodes exist for FindPointsGSLIB
    pmesh.EnsureNodes();
    MFEM_VERIFY(pmesh.GetNodes() != nullptr, "Mesh nodes are required for FindPointsGSLIB");
    // Build point samples  
    mfem::Vector point_pos = BuildPointPositionsByNodes(space_dim, Nx, Ny, Nz, bbmin3, out.spacing);

    // Locate points in the mesh
    mfem::FindPointsGSLIB finder(pmesh.GetComm());
    finder.Setup(pmesh);
    finder.FindPoints(point_pos, mfem::Ordering::byNODES);

    const auto &code = finder.GetCode();
    const auto &elid = finder.GetElem();
    const mfem::Vector &ref = finder.GetReferencePosition(); // component-major: ref[d*N + p]

    // Sample either raw -grad(V) or projected H1 E_h
    if (!H1_project)
    {
        SampleGradV(pmesh, V, code, elid, ref, N, accept_surface_projection, out);
    }
    else
    {
        const ProjectedH1VectorField proj = BuildEFieldH1Projection(pmesh, V, cfg);
        SampleH1Projection(pmesh, *(proj.E), code, elid, ref, N, accept_surface_projection, out);
    }
}




static inline void h5_check(const herr_t st, const char *what)
{
  if (st < 0) { throw std::runtime_error(std::string("HDF5 error: ") + what); }
}

static void WriteGridSampleBinary(const GridSample &g,
                                  const std::filesystem::path &out_dir,
                                  const std::string &filename = "E_interpolated.h5")
{
  const hsize_t Nx  = static_cast<hsize_t>(g.Nx);
  const hsize_t Ny  = static_cast<hsize_t>(g.Ny);
  const hsize_t Nz  = static_cast<hsize_t>(g.Nz);
  const hsize_t dim = static_cast<hsize_t>(g.dim);

  const size_t N = static_cast<size_t>(g.Nx) * static_cast<size_t>(g.Ny) * static_cast<size_t>(g.Nz);

  if (g.dim != 2 && g.dim != 3)
    throw std::runtime_error("WriteGridSampleBinary: g.dim must be 2 or 3.");

  if (g.Ex.size() != N || g.Ey.size() != N || (g.dim == 3 && g.Ez.size() != N))
    throw std::runtime_error("WriteGridSampleBinary: Ex/Ey/Ez sizes do not match grid size.");

  // Coordinates are ALWAYS treated as 3D (x,y,z), regardless of g.dim.
  // For 2D geometry you still must provide origin/spacing as (x,y,z),
  // with Nz possibly 1 if the geometry is truly 2D-in-(x,y).
  mfem::Vector origin(3), spacing(3);
  if (g.origin.Size() == g.spacing.Size())
  {
    const int n = g.origin.Size();

    if (n < 1 || n > 3)
      throw std::runtime_error("WriteGridSampleBinary: origin/spacing must have length 1, 2, or 3.");

    origin = 0.0;
    spacing = 0.0;

    for (int d = 0; d < n; ++d)
    {
      origin[d]  = g.origin[d];
      spacing[d] = g.spacing[d];
    }
  }
  else
  {
    throw std::runtime_error("WriteGridSampleBinary: origin and spacing must have the same length.");
}
  std::filesystem::create_directories(out_dir);
  const auto out_path = out_dir / filename;

  // ------------------------------------------------------------------
  // Pack E into contiguous buffer for /E/field with shape [Nz, Ny, Nx, dim]
  // p = i + Nx*(j + Ny*k)
  // buffer offset = p*dim + c
  // ------------------------------------------------------------------
  std::vector<double> E(N * static_cast<size_t>(g.dim), std::numeric_limits<double>::quiet_NaN());
  for (size_t p = 0; p < N; ++p)
  {
    E[p * static_cast<size_t>(g.dim) + 0] = g.Ex[p];
    E[p * static_cast<size_t>(g.dim) + 1] = g.Ey[p];
    if (g.dim == 3) { E[p * static_cast<size_t>(g.dim) + 2] = g.Ez[p]; }
  }

  // ------------------------------------------------------------------
  // Create file
  // ------------------------------------------------------------------
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (fapl < 0) throw std::runtime_error("HDF5 error: H5Pcreate(FILE_ACCESS)");

  h5_check(H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG), "H5Pset_fclose_degree");

  hid_t file = H5Fcreate(out_path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  if (file < 0) { H5Pclose(fapl); throw std::runtime_error("HDF5 error: H5Fcreate"); }
  H5Pclose(fapl);

  // ------------------------------------------------------------------
  // Helpers
  // ------------------------------------------------------------------
  auto write_attr_i32 = [&](hid_t obj, const char *name, int value)
  {
    hsize_t one[1] = {1};
    hid_t as = H5Screate_simple(1, one, nullptr);
    if (as < 0) throw std::runtime_error("HDF5 error: H5Screate_simple(attr i32)");

    hid_t at = H5Acreate2(obj, name, H5T_STD_I32LE, as, H5P_DEFAULT, H5P_DEFAULT);
    if (at < 0) { H5Sclose(as); throw std::runtime_error("HDF5 error: H5Acreate2(attr i32)"); }

    h5_check(H5Awrite(at, H5T_NATIVE_INT, &value), "H5Awrite(i32)");
    H5Aclose(at);
    H5Sclose(as);
  };

  auto write_attr_f64_vec = [&](hid_t obj, const char *name, const mfem::Vector &v)
  {
    hsize_t n[1] = {static_cast<hsize_t>(v.Size())};
    hid_t as = H5Screate_simple(1, n, nullptr);
    if (as < 0) throw std::runtime_error("HDF5 error: H5Screate_simple(attr f64 vec)");

    hid_t at = H5Acreate2(obj, name, H5T_IEEE_F64LE, as, H5P_DEFAULT, H5P_DEFAULT);
    if (at < 0) { H5Sclose(as); throw std::runtime_error("HDF5 error: H5Acreate2(attr f64 vec)"); }

    h5_check(H5Awrite(at, H5T_NATIVE_DOUBLE, v.GetData()), "H5Awrite(f64vec)");
    H5Aclose(at);
    H5Sclose(as);
  };

  auto write_1d_f64 = [&](hid_t parent, const char *name, const double *data, hsize_t n)
  {
    hsize_t d1[1] = {n};
    hid_t sp = H5Screate_simple(1, d1, nullptr);
    if (sp < 0) throw std::runtime_error("HDF5 error: H5Screate_simple(1D)");

    hid_t ds = H5Dcreate2(parent, name, H5T_IEEE_F64LE, sp,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(sp); throw std::runtime_error(std::string("HDF5 error: H5Dcreate2(/grid/") + name + ")"); }

    const std::string what = std::string("H5Dwrite(/grid/") + name + ")";
    h5_check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
             what.c_str());

    H5Dclose(ds);
    H5Sclose(sp);
  };

  try
  {
    // ------------------------------------------------------------------
    // /E/field
    // ------------------------------------------------------------------
    hid_t grpE = H5Gcreate2(file, "/E", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (grpE < 0) throw std::runtime_error("HDF5 error: H5Gcreate2(/E)");

    hsize_t d4[4] = {Nz, Ny, Nx, dim};
    hid_t spaceE = H5Screate_simple(4, d4, nullptr);
    if (spaceE < 0) { H5Gclose(grpE); throw std::runtime_error("HDF5 error: H5Screate_simple(/E/field)"); }

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0) { H5Sclose(spaceE); H5Gclose(grpE); throw std::runtime_error("HDF5 error: H5Pcreate(DCPL)"); }

    hsize_t chunk4[4] = {
      std::min<hsize_t>(Nz, 16),
      std::min<hsize_t>(Ny, 64),
      std::min<hsize_t>(Nx, 64),
      dim
    };
    h5_check(H5Pset_chunk(dcpl, 4, chunk4), "H5Pset_chunk");

    hid_t dsetE = H5Dcreate2(grpE, "field", H5T_IEEE_F64LE, spaceE,
                            H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (dsetE < 0)
    {
      H5Pclose(dcpl);
      H5Sclose(spaceE);
      H5Gclose(grpE);
      throw std::runtime_error("HDF5 error: H5Dcreate2(/E/field)");
    }

    h5_check(H5Dwrite(dsetE, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, E.data()),
             "H5Dwrite(/E/field)");

    write_attr_i32(dsetE, "dim", g.dim);
    write_attr_i32(dsetE, "Nx",  g.Nx);
    write_attr_i32(dsetE, "Ny",  g.Ny);
    write_attr_i32(dsetE, "Nz",  g.Nz);
    write_attr_f64_vec(dsetE, "origin",  origin);   // MUST be (x,y,z)
    write_attr_f64_vec(dsetE, "spacing", spacing);  // MUST be (dx,dy,dz)

    // Optional: record dataset axis meaning explicitly: field[z,y,x,comp]
    {
      int axis_order[4] = {2, 1, 0, 3};
      hsize_t n[1] = {4};
      hid_t as = H5Screate_simple(1, n, nullptr);
      hid_t at = H5Acreate2(dsetE, "axis_order", H5T_STD_I32LE, as, H5P_DEFAULT, H5P_DEFAULT);
      h5_check(H5Awrite(at, H5T_NATIVE_INT, axis_order), "H5Awrite(axis_order)");
      H5Aclose(at); H5Sclose(as);
    }

    H5Dclose(dsetE);
    H5Pclose(dcpl);
    H5Sclose(spaceE);
    H5Gclose(grpE);

    // ------------------------------------------------------------------
    // /grid (ALWAYS x,y,z)
    // ------------------------------------------------------------------
    std::vector<double> x(static_cast<size_t>(g.Nx));
    std::vector<double> y(static_cast<size_t>(g.Ny));
    std::vector<double> z(static_cast<size_t>(g.Nz));

    const double x0 = origin[0], dx = spacing[0];
    const double y0 = origin[1], dy = spacing[1];
    const double z0 = origin[2], dz = spacing[2];

    for (int i = 0; i < g.Nx; ++i) x[static_cast<size_t>(i)] = x0 + dx * static_cast<double>(i);
    for (int j = 0; j < g.Ny; ++j) y[static_cast<size_t>(j)] = y0 + dy * static_cast<double>(j);
    for (int k = 0; k < g.Nz; ++k) z[static_cast<size_t>(k)] = z0 + dz * static_cast<double>(k);

    hid_t grpG = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (grpG < 0) throw std::runtime_error("HDF5 error: H5Gcreate2(/grid)");

    write_1d_f64(grpG, "x", x.data(), Nx);
    write_1d_f64(grpG, "y", y.data(), Ny);
    write_1d_f64(grpG, "z", z.data(), Nz);

    // Optional: coordinate kind (x,y,z)
    {
      int coord_kind[3] = {0, 1, 2}; // 0=x,1=y,2=z
      hsize_t n[1] = {3};
      hid_t as = H5Screate_simple(1, n, nullptr);
      hid_t at = H5Acreate2(grpG, "coord_kind", H5T_STD_I32LE, as, H5P_DEFAULT, H5P_DEFAULT);
      h5_check(H5Awrite(at, H5T_NATIVE_INT, coord_kind), "H5Awrite(coord_kind)");
      H5Aclose(at); H5Sclose(as);
    }

    H5Gclose(grpG);

    H5Fclose(file);
  }
  catch (...)
  {
    H5Fclose(file);
    throw;
  }
}


int do_interpolate(Config cfg)
{
  // Mess with mpi so only Rank 0 machine does interpolation
  int mpi_inited = 0;
  MPI_Initialized(&mpi_inited);
  int world_rank = 0, world_size = 1;
  if (mpi_inited)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  }
  // Rank-0-only execution
  if (world_rank != 0) { return 0; }
  MPI_Comm comm_for_single_rank = MPI_COMM_SELF;

  // Get filepath
  std::filesystem::path save_root(cfg.save_path);
  std::filesystem::path run_dir;
  // TODOO Make multiple runs possible same TODO in optimization move this bit to utils module
  auto runs = read_root_level_runs_from_meta(save_root);
  if (runs.empty()) {
      run_dir = save_root;
      std::cout << "[METRICS] meta.txt not found or empty; "
                << "using save_root directly: " << run_dir << "\n";
  } else if (runs.size() == 1) {
      run_dir = save_root / runs.front();
      std::cout << "[METRICS] Found single run in meta.txt: "
                << runs.front() << " → " << run_dir << "\n";
  } else {
      std::cerr << "[METRICS] meta.txt contains multiple runs under "
                << save_root << ":\n";
      for (const auto &r : runs) {
          std::cerr << "  - " << r << "\n";
      }
      std::cerr << "Please specify which run to use (CLI/config).\n";
      std::exit(1);
  }
  // Load data 
  SimulationResult result = load_results(cfg, run_dir);

  // Interpolation grid settings
  const int Nx = cfg.interp.Nx;
  const int Ny = cfg.interp.Ny;
  const int Nz = cfg.interp.Nz; // ignored if dim==2
  const bool accept_surface = false;//cfg.interp_accept_surface_projection;

  GridSample grid;
  grid.dim = result.mesh->Dimension();

  // Interpolate 
  SampleEFieldOnCartesianGrid(*result.mesh, *result.V, Nx, Ny, Nz, grid, cfg, cfg.interp.H1_project, accept_surface);

  // Write output
  const auto out_dir = run_dir / "interpolated";
  std::filesystem::create_directories(out_dir);
  WriteGridSampleBinary(grid, out_dir); 

  std::cout << "[INTERP] wrote interpolated grid to: " << out_dir << "\n";
  return 0;
}