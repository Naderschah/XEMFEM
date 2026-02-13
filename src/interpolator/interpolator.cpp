#include "interpolator.h"

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
                        mfem::FindPointsGSLIB &finder,
                        const mfem::Array<unsigned int> &code, // size N, on all ranks
                        const int N,
                        const bool accept_surface_projection,
                        GridSample &out)
{
    using namespace mfem;

    const int field_dim = pmesh.Dimension();
    const int vdim = field_dim;

    // 1) Redistribute point info to owning ranks (local elem ids + local ref coords + local codes)
    Array<unsigned int> recv_elem;
    Array<unsigned int> recv_code;
    Vector recv_ref_by_vdim; // component-major: recv_ref[d*nrecv + i]
    finder.DistributePointInfoToOwningMPIRanks(recv_elem, recv_ref_by_vdim, recv_code);

    const int nrecv = recv_elem.Size();
    MFEM_VERIFY(recv_code.Size() == nrecv, "recv_code size mismatch");
    MFEM_VERIFY(recv_ref_by_vdim.Size() == nrecv * field_dim, "recv_ref size mismatch");

    // 2) Evaluate locally on owning ranks
    GradientGridFunctionCoefficient gradV(&V);
    Vector g(field_dim);
    IntegrationPoint ip;

    // Local interpolation buffer for the received points.
    // Use byNODES (point-major): [p0: x y z][p1: x y z]...
    Vector vals_local;
    vals_local.SetSize(nrecv * vdim);
    vals_local = 0.0;

    for (int i = 0; i < nrecv; ++i)
    {
        const unsigned int c = recv_code[i];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf)) { continue; }

        const double r0 = recv_ref_by_vdim[0 * nrecv + i];
        const double r1 = recv_ref_by_vdim[1 * nrecv + i];
        const double r2 = (field_dim == 3) ? recv_ref_by_vdim[2 * nrecv + i] : 0.0;

        if (field_dim == 2) { ip.Set2(r0, r1); }
        else                { ip.Set3(r0, r1, r2); }

        const int e = static_cast<int>(recv_elem[i]);
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_VERIFY(T != nullptr, "Null ElementTransformation");

        T->SetIntPoint(&ip);
        gradV.Eval(g, *T, ip);
        g *= -1.0;

        // byNODES layout
        vals_local[i * vdim + 0] = g[0];
        vals_local[i * vdim + 1] = g[1];
        if (field_dim == 3) { vals_local[i * vdim + 2] = g[2]; }
    }

    // 3) Gather interpolated values back to the ranks that originated the query
    Vector vals_global;
    vals_global.SetSize(N * vdim);
    vals_global = 0.0;

    finder.DistributeInterpolatedValues(vals_local, vdim, Ordering::byNODES, vals_global);

    // 4) Fill outputs on each rank (use original 'code' to set validity consistently)
    //    Assumes out.* arrays are sized to N and initialized (or we set them here).
    for (int p = 0; p < N; ++p)
    {
        const unsigned int c = code[p];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf))
        {
            out.valid[p] = 0;
            continue;
        }

        const double ex = vals_global[p * vdim + 0];
        const double ey = vals_global[p * vdim + 1];
        const double ez = (field_dim == 3) ? vals_global[p * vdim + 2] : 0.0;

        out.valid[p] = 1;
        out.Ex[p] = ex;
        out.Ey[p] = ey;
        if (field_dim == 3) { out.Ez[p] = ez; }
        out.Emag[p] = std::sqrt(ex*ex + ey*ey + ez*ez);
    }
}


static void SampleH1Projection(mfem::ParMesh &pmesh,
                               const mfem::ParGridFunction &E_h, // vdim == field_dim
                               mfem::FindPointsGSLIB &finder,
                               const mfem::Array<unsigned int> &code, // size N (original query ordering)
                               const int N,
                               const bool accept_surface_projection,
                               GridSample &out)
{
    using namespace mfem;

    const int field_dim = pmesh.Dimension();
    const int vdim = field_dim;

    // 1) Redistribute point info to owning ranks
    Array<unsigned int> recv_elem;
    Array<unsigned int> recv_code;
    Vector recv_ref_by_vdim; // component-major: ref[d*nrecv + i]
    finder.DistributePointInfoToOwningMPIRanks(recv_elem, recv_ref_by_vdim, recv_code);

    const int nrecv = recv_elem.Size();
    MFEM_VERIFY(recv_code.Size() == nrecv, "recv_code size mismatch");
    MFEM_VERIFY(recv_ref_by_vdim.Size() == nrecv * field_dim, "recv_ref size mismatch");

    // 2) Evaluate locally on owning ranks
    Vector g(field_dim);
    IntegrationPoint ip;

    // Local interpolation buffer for received points, byNODES (point-major)
    Vector vals_local(nrecv * vdim);
    vals_local = 0.0;

    for (int i = 0; i < nrecv; ++i)
    {
        const unsigned int c = recv_code[i];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf)) { continue; }

        const double r0 = recv_ref_by_vdim[0 * nrecv + i];
        const double r1 = recv_ref_by_vdim[1 * nrecv + i];
        const double r2 = (field_dim == 3) ? recv_ref_by_vdim[2 * nrecv + i] : 0.0;

        if (field_dim == 2) { ip.Set2(r0, r1); }
        else                { ip.Set3(r0, r1, r2); }

        const int e = static_cast<int>(recv_elem[i]);
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_VERIFY(T != nullptr, "Null ElementTransformation");

        T->SetIntPoint(&ip);

        // Evaluate the projected continuous vector field at this point
        E_h.GetVectorValue(*T, ip, g);

        // byNODES layout
        vals_local[i * vdim + 0] = g[0];
        vals_local[i * vdim + 1] = g[1];
        if (field_dim == 3) { vals_local[i * vdim + 2] = g[2]; }
    }

    // 3) Gather interpolated values back to original querying ranks/order
    Vector vals_global(N * vdim);
    vals_global = 0.0;

    finder.DistributeInterpolatedValues(vals_local, vdim, Ordering::byNODES, vals_global);

    // 4) Fill outputs on each rank using original 'code' for validity
    for (int p = 0; p < N; ++p)
    {
        const unsigned int c = code[p];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf))
        {
            out.valid[p] = 0;
            continue;
        }

        const double ex = vals_global[p * vdim + 0];
        const double ey = vals_global[p * vdim + 1];
        const double ez = (field_dim == 3) ? vals_global[p * vdim + 2] : 0.0;

        out.valid[p] = 1;
        out.Ex[p] = ex;
        out.Ey[p] = ey;
        if (field_dim == 3) { out.Ez[p] = ez; }
        out.Emag[p] = std::sqrt(ex*ex + ey*ey + ez*ez);
    }
}

static mfem::Vector BuildPointPositionsByNodesRange(const int space_dim,
                                                    const int Nx,
                                                    const int Ny,
                                                    const int Nz,
                                                    const mfem::Vector &bbmin3,
                                                    const mfem::Vector &spacing,
                                                    const int p0,
                                                    const int nloc)
{
    using namespace mfem;

    Vector point_pos(space_dim * nloc);

    for (int lp = 0; lp < nloc; ++lp)
    {
        const int p = p0 + lp;

        const int i = p % Nx;
        const int t = p / Nx;
        const int j = t % Ny;
        const int k = t / Ny;

        point_pos[0 * nloc + lp] = bbmin3[0] + i * spacing[0];
        point_pos[1 * nloc + lp] = bbmin3[1] + j * spacing[1];
        if (space_dim == 3)
        {
            point_pos[2 * nloc + lp] = bbmin3[2] + k * spacing[2];
        }
    }

    return point_pos;
}
void SampleEFieldOnCartesianGrid(mfem::ParMesh &pmesh,
                                 const mfem::ParGridFunction &V,
                                 const int Nx, const int Ny, const int Nz,
                                 GridSample &out,
                                 Config cfg,
                                 const bool H1_project,    // IGNORED – always direct gradient
                                 const bool accept_surface_projection) // IGNORED
{
    using namespace mfem;

    MPI_Comm comm = pmesh.GetComm();
    int myid = 0, np = 1;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &np);

    const int space_dim = pmesh.SpaceDimension();
    const int field_dim = pmesh.Dimension();
    MFEM_VERIFY(space_dim == 2 || space_dim == 3, "Only 2D/3D physical space supported");

    // ------------------------------------------------------------------
    // Global bounding box
    // ------------------------------------------------------------------
    Vector bbmin_s(space_dim), bbmax_s(space_dim);
    pmesh.GetBoundingBox(bbmin_s, bbmax_s);

    Vector bbmin3(3), bbmax3(3);
    bbmin3 = 0.0; bbmax3 = 0.0;
    for (int d = 0; d < space_dim; ++d)
    {
        bbmin3[d] = bbmin_s[d];
        bbmax3[d] = bbmax_s[d];
    }

    Vector origin(3), spacing(3);
    origin = bbmin3;
    spacing = 0.0;

    MFEM_VERIFY(Nx > 1 && Ny > 1, "Nx and Ny must be > 1");
    spacing[0] = (bbmax3[0] - bbmin3[0]) / (Nx - 1);
    spacing[1] = (bbmax3[1] - bbmin3[1]) / (Ny - 1);
    if (Nz > 1 && space_dim == 3)
    {
        spacing[2] = (bbmax3[2] - bbmin3[2]) / (Nz - 1);
    }
    else
    {
        spacing[2] = 0.0;
    }

    const int N = Nx * Ny * Nz;

    // ------------------------------------------------------------------
    // Partition global points
    // ------------------------------------------------------------------
    const int base = (np > 0) ? (N / np) : 0;
    const int rem  = (np > 0) ? (N % np) : 0;
    const int nloc = base + (myid < rem ? 1 : 0);
    const int p0   = myid * base + (myid < rem ? myid : rem);

    // ------------------------------------------------------------------
    // Local query points – dummy if nloc == 0
    // ------------------------------------------------------------------
    const bool need_dummy = (nloc == 0);
    const int nq = need_dummy ? 1 : nloc;

    // Create vector of coordinates in "byNODES" layout (all x, then all y, then all z)
    Vector point_pos(space_dim * nq);   // FIXED: no ordering enum in constructor

    if (!need_dummy)
    {
        // Build local subset deterministically in byNODES order
        for (int lp = 0; lp < nq; ++lp)
        {
            const int p = p0 + lp;
            const int i = p % Nx;
            const int t = p / Nx;
            const int j = t % Ny;
            const int k = t / Ny;

            point_pos[0 * nq + lp] = bbmin3[0] + i * spacing[0];
            point_pos[1 * nq + lp] = bbmin3[1] + j * spacing[1];
            if (space_dim == 3)
                point_pos[2 * nq + lp] = bbmin3[2] + k * spacing[2];
        }
    }
    else
    {
        // Dummy point inside bbox, still byNODES layout
        point_pos = 0.0;
        point_pos[0] = bbmin3[0];
        if (space_dim >= 2) point_pos[1] = bbmin3[1];
        if (space_dim == 3) point_pos[2] = bbmin3[2];
    }

    // ------------------------------------------------------------------
    // Local output containers (size nq, but only first nloc are meaningful)
    // ------------------------------------------------------------------
    GridSample local;
    local.dim = field_dim;
    local.Nx = Nx; local.Ny = Ny; local.Nz = Nz;
    local.origin.SetSize(3); local.spacing.SetSize(3);
    local.origin = origin; local.spacing = spacing;

    local.Ex.assign(nq, std::numeric_limits<double>::quiet_NaN());
    local.Ey.assign(nq, std::numeric_limits<double>::quiet_NaN());
    local.Ez.assign(nq, 0.0);
    local.Emag.assign(nq, std::numeric_limits<double>::quiet_NaN());
    local.valid.assign(nq, static_cast<unsigned char>(0));

    // ------------------------------------------------------------------
    // Point location (collective) – explicitly specify byNODES layout
    // ------------------------------------------------------------------
    pmesh.EnsureNodes();
    MFEM_VERIFY(pmesh.GetNodes() != nullptr, "Mesh nodes are required for FindPointsGSLIB");

    FindPointsGSLIB finder(comm);
    finder.Setup(pmesh);
    finder.FindPoints(point_pos, Ordering::byNODES);   // FIXED: ordering passed here

    // ------------------------------------------------------------------
    // --- PORTED LOGIC: distribute to owning ranks, evaluate E = -∇V ---
    // ------------------------------------------------------------------
    Array<unsigned int> recv_elem, recv_code;
    Vector recv_ref;   // reference coordinates, returned in byVDIM layout
    finder.DistributePointInfoToOwningMPIRanks(recv_elem, recv_ref, recv_code);

    const int nrecv = recv_elem.Size();
    const int sdim = space_dim;

    // Buffer for E field on the points we received (byVDIM)
    Vector recv_E(nrecv * sdim);
    Vector grad(sdim);

    for (int j = 0; j < nrecv; ++j)
    {
        if (recv_code[j] == 2)
        {
            recv_E(j * sdim + 0) = std::numeric_limits<double>::quiet_NaN();
            recv_E(j * sdim + 1) = std::numeric_limits<double>::quiet_NaN();
            if (sdim == 3)
                recv_E(j * sdim + 2) = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        const int elem_id = (int)recv_elem[j];
        ElementTransformation *T = pmesh.GetElementTransformation(elem_id);
        MFEM_VERIFY(T != nullptr, "Null ElementTransformation");

        IntegrationPoint ip;
        if (sdim == 2)
            ip.Set2(recv_ref(j * sdim + 0), recv_ref(j * sdim + 1));
        else // sdim == 3
            ip.Set3(recv_ref(j * sdim + 0), recv_ref(j * sdim + 1), recv_ref(j * sdim + 2));

        T->SetIntPoint(&ip);
        V.GetGradient(*T, grad);
        grad *= -1.0;   // E = -∇φ

        recv_E(j * sdim + 0) = grad[0];
        recv_E(j * sdim + 1) = grad[1];
        if (sdim == 3)
            recv_E(j * sdim + 2) = grad[2];
    }

    // Scatter the computed E field back to the original ranks / original point order
    // DistributeInterpolatedValues always returns in byVDIM layout.
    Vector Evals_flat(nq * sdim);
    finder.DistributeInterpolatedValues(recv_E, sdim, Ordering::byVDIM, Evals_flat);

    // ------------------------------------------------------------------
    // Fill local arrays for the first nloc real points
    // ------------------------------------------------------------------
    const Array<unsigned int> &code_local = finder.GetCode(); // original point status
    for (int lp = 0; lp < nloc; ++lp)
    {
        const double ex = Evals_flat(lp * sdim + 0);
        const double ey = Evals_flat(lp * sdim + 1);
        const double ez = (sdim == 3) ? Evals_flat(lp * sdim + 2) : 0.0;

        local.Ex[lp] = ex;
        local.Ey[lp] = ey;
        local.Ez[lp] = ez;

        double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
        local.Emag[lp] = std::isfinite(emag) ? emag : std::numeric_limits<double>::quiet_NaN();

        // Valid if the point was successfully located (code == 0)
        local.valid[lp] = (code_local[lp] == 0) ? 1 : 0;
    }

    // ------------------------------------------------------------------
    // Gather to rank 0 (same as original)
    // ------------------------------------------------------------------
    const int send_n = nloc;

    std::vector<int> recvcounts(np, 0), displs(np, 0);
    for (int r = 0; r < np; ++r)
    {
        const int rn  = base + (r < rem ? 1 : 0);
        const int rp0 = r * base + (r < rem ? r : rem);
        recvcounts[r] = rn;
        displs[r]     = rp0;
    }

    if (myid == 0)
    {
        out.dim = field_dim;
        out.Nx = Nx; out.Ny = Ny; out.Nz = Nz;
        out.origin.SetSize(3); out.spacing.SetSize(3);
        out.origin = origin; out.spacing = spacing;

        out.Ex.assign(N, std::numeric_limits<double>::quiet_NaN());
        out.Ey.assign(N, std::numeric_limits<double>::quiet_NaN());
        out.Ez.assign(N, 0.0);
        out.Emag.assign(N, std::numeric_limits<double>::quiet_NaN());
        out.valid.assign(N, static_cast<unsigned char>(0));
    }

    double dummy_d = 0.0;
    unsigned char dummy_uc = 0;

    const double *sendEx   = (send_n > 0) ? local.Ex.data()   : nullptr;
    const double *sendEy   = (send_n > 0) ? local.Ey.data()   : nullptr;
    const double *sendEz   = (send_n > 0) ? local.Ez.data()   : nullptr;
    const double *sendEmag = (send_n > 0) ? local.Emag.data() : nullptr;
    const unsigned char *sendVal = (send_n > 0) ? local.valid.data() : nullptr;

    double *recvEx   = (myid == 0) ? out.Ex.data()   : &dummy_d;
    double *recvEy   = (myid == 0) ? out.Ey.data()   : &dummy_d;
    double *recvEz   = (myid == 0) ? out.Ez.data()   : &dummy_d;
    double *recvEmag = (myid == 0) ? out.Emag.data() : &dummy_d;
    unsigned char *recvVal = (myid == 0) ? out.valid.data() : &dummy_uc;

    MPI_Gatherv(sendEx, send_n, MPI_DOUBLE,
                recvEx, recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, comm);
    MPI_Gatherv(sendEy, send_n, MPI_DOUBLE,
                recvEy, recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, comm);
    if (field_dim == 3)
    {
        MPI_Gatherv(sendEz, send_n, MPI_DOUBLE,
                    recvEz, recvcounts.data(), displs.data(), MPI_DOUBLE,
                    0, comm);
    }
    MPI_Gatherv(sendEmag, send_n, MPI_DOUBLE,
                recvEmag, recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, comm);
    MPI_Gatherv(sendVal, send_n, MPI_UNSIGNED_CHAR,
                recvVal, recvcounts.data(), displs.data(), MPI_UNSIGNED_CHAR,
                0, comm);
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
  int world_rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  const auto targets = targets_from_save_root(cfg);

  for (const auto& run_dir : targets)
  {
    // Load data (must be called on all ranks if it builds ParMesh/ParGridFunction on world comm)
    SimulationResult result = load_results(cfg, run_dir);

    const int Nx = cfg.interp.Nx;
    const int Ny = cfg.interp.Ny;
    const int Nz = cfg.interp.Nz;
    const bool accept_surface = false;

    GridSample grid;
    grid.dim = result.mesh->Dimension();

    // Interpolate (collective across result.mesh->GetComm())
    SampleEFieldOnCartesianGrid(*result.mesh, *result.V,
                                Nx, Ny, Nz, grid,
                                cfg, cfg.interp.H1_project, accept_surface);

    // Rank-0-only write
    if (world_rank == 0)
    {
      const auto out_dir = run_dir / "interpolated";
      std::filesystem::create_directories(out_dir);
      WriteGridSampleBinary(grid, out_dir);
      std::cout << "[INTERP] wrote interpolated grid to: " << out_dir << "\n";
    }
  }
  return 0;
}
