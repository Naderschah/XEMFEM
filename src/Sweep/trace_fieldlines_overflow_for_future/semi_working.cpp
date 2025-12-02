//Headers 


struct DriftTimeEstimate
{
    double H            = 0.0;  // drift height
    double E_rms        = 0.0;  // RMS field in tracking region
    double a_typ        = 0.0;  // |q/m| * E_rms
    double t_drift      = 0.0;  // time to cross H under a_typ
    double dt_traverse  = 0.0;  // base dt giving T_trav drift crossings in max_steps
    double crossings_est= 0.0;  // estimated crossings with dt_traverse (= T_trav)
};
// Holds data for the binary CIV search
struct ElementGeometryData
{
    std::vector<double> r_centroid;
    std::vector<double> z_centroid;
    std::vector<double> volume;
};
struct ColumnRefineResult
{
    bool   found;
    double z_ci;              // refined highest CI height
    bool   is_cathode_ci;
    bool   is_wall_ci;
};
struct CIVolumeResult
{
    double V_total; // total volume in the considered region
    double V_CI;    // charge-insensitive volume (cathode + wall)
};
struct InitialColumnCIResult
{
    bool   found;
    double r;
    double z_ci;          // highest CI height at r
    int    elem;          // element containing (r,z_ci)
    bool   cathode_ci;    // seen cathode CI
    bool   wall_ci;       // seen wall CI
};
struct ColumnCISearchResult
{
    bool   found;             // whether any CI was detected
    double z_ci;              // highest z with CI along this column
    bool   is_cathode_ci;     // true if CI dominated by cathode at this column
    bool   is_wall_ci;        // true if CI dominated by wall at this column
};

struct ProbeRequest
{
    double r;
    double z;
    int    elem_guess;   // optional (-1 means unknown)
};
struct ProbeResult
{
    bool   ok         = false; // mapping + tracing succeeded
    bool   is_ci      = false; // IsChargeInsensitive(...)
    bool   is_cathode = false; // IsCathodeInsensitive(...)
    bool   is_wall    = false; // IsWallInsensitive(...)

    double r          = 0.0;   // probed (r,z)
    double z          = 0.0;
    int    elem       = -1;    // element id returned/used by tracer
    int    exit_code  = -1;    // int(ElectronExitCode) for debugging
};

// COde
static std::vector<std::vector<int>> BuildVertexAdjacency(mfem::ParMesh &mesh)
{
    using namespace mfem;

    const int ne = mesh.GetNE();
    const int nv = mesh.GetNV();

    // vertex -> elements containing that vertex
    std::vector<std::vector<int>> vertex_to_elems(nv);

    Array<int> verts;
    for (int e = 0; e < ne; ++e)
    {
        mesh.GetElementVertices(e, verts);
        for (int k = 0; k < verts.Size(); ++k)
        {
            int v = verts[k];
            MFEM_ASSERT(v >= 0 && v < nv, "Invalid vertex index");
            vertex_to_elems[v].push_back(e);
        }
    }

    // element -> vertex neighbors (all elems sharing any vertex)
    std::vector<std::vector<int>> elem_vertex_neighbors(ne);

    for (int e = 0; e < ne; ++e)
    {
        mesh.GetElementVertices(e, verts);
        auto &nbrs = elem_vertex_neighbors[e];

        for (int k = 0; k < verts.Size(); ++k)
        {
            int v = verts[k];
            for (int ee : vertex_to_elems[v])
            {
                if (ee != e)
                {
                    nbrs.push_back(ee);
                }
            }
        }

        // optional: deduplicate
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    return elem_vertex_neighbors;
}

// Return true if this exit is cathode-based charge insensitivity
static bool IsCathodeInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitCathode);
}

// Return true if this exit is wall-based charge insensitivity
static bool IsWallInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitWall);
}
static std::vector<std::vector<int>>
BuildUnifiedAdjacency(const ElementAdjacency &adj,
                      const std::vector<std::vector<int>> &vertex_adj,
                      int ne)
{
    std::vector<std::vector<int>> unified(ne);

    for (int e = 0; e < ne; ++e)
    {
        auto &u = unified[e];

        // insert face neighbors
        const auto &fn = adj.neighbors[e];
        u.insert(u.end(), fn.begin(), fn.end());

        // insert vertex neighbors
        const auto &vn = vertex_adj[e];
        u.insert(u.end(), vn.begin(), vn.end());

        // remove duplicates & self
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), u.end()), u.end());

        // safety: remove the element itself if present
        u.erase(std::remove(u.begin(), u.end(), e), u.end());
    }

    return unified;
}

// Axisymmetric element centroids and volumes (no 2π factor)
static ElementGeometryData ComputeElementGeometry(mfem::ParMesh &mesh)
{
    using namespace mfem;

    const int ne  = mesh.GetNE();
    const int dim = mesh.SpaceDimension();

    ElementGeometryData geom;
    geom.r_centroid.resize(ne);
    geom.z_centroid.resize(ne);
    geom.volume.resize(ne);

    IntegrationRules irs;

    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = mesh.GetElementTransformation(e);
        MFEM_ASSERT(T, "ComputeElementGeometry: ElementTransformation is null");

        Geometry::Type geom_type = T->GetGeometryType();
        // Order 2 is enough for centroid and volume here
        const IntegrationRule &ir = irs.Get(geom_type, 2);
        const int n_ip            = ir.GetNPoints();

        double V_e   = 0.0;
        double r_acc = 0.0;
        double z_acc = 0.0;

        for (int i = 0; i < n_ip; ++i)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            Vector x(dim);
            T->Transform(ip, x);
            const double r = x(0);
            const double z = x(1);

            const double detJ  = T->Weight();
            const double w_ref = ip.weight;
            const double dV    = r * detJ * w_ref; // axisymmetric, no 2π

            V_e   += dV;
            r_acc += r * dV;
            z_acc += z * dV;
        }

        if (V_e > 0.0)
        {
            geom.r_centroid[e] = r_acc / V_e;
            geom.z_centroid[e] = z_acc / V_e;
        }
        else
        {
            geom.r_centroid[e] = 0.0;
            geom.z_centroid[e] = 0.0;
        }

        geom.volume[e] = V_e;
    }

    return geom;
}

// Evaluate CI boundary height z_b at radius r using linear interpolation
// over the given interface points (assumed sorted by r).
// If r is outside the sampled range, we extrapolate using the nearest point.
static double EvaluateInterfaceZ(const std::vector<InterfacePoint> &iface,
                                 double                             r)
{
    const std::size_t n = iface.size();
    if (n == 0)
    {
        return 0.0;
    }

    if (n == 1)
    {
        return iface[0].z_ci;
    }

    // If r is below the smallest sampled radius, clamp to first point
    if (r <= iface.front().r)
    {
        return iface.front().z_ci;
    }

    // If r is above the largest sampled radius, clamp to last point
    if (r >= iface.back().r)
    {
        return iface.back().z_ci;
    }

    // Find segment [i, i+1] s.t. iface[i].r <= r <= iface[i+1].r
    for (std::size_t i = 0; i + 1 < n; ++i)
    {
        const double r0 = iface[i].r;
        const double r1 = iface[i+1].r;
        if (r >= r0 && r <= r1)
        {
            const double z0 = iface[i].z_ci;
            const double z1 = iface[i+1].z_ci;

            const double t = (r - r0) / (r1 - r0);
            return (1.0 - t) * z0 + t * z1;
        }
    }

    // Fallback (should not reach here if above logic is correct)
    return iface.back().z_ci;
}


// Parallel interface
static void ProbeBatch(const std::vector<ProbeRequest>       &requests,
                       const std::vector<TracerContext*>     &tctx_pool,
                       std::vector<ProbeResult>              &out_results)
{
    const int N = (int)requests.size();
    out_results.clear();
    out_results.resize(N);

    if (N == 0 || tctx_pool.empty())
    {
        return;
    }

    const int n_ctx = (int)tctx_pool.size();
    const int max_threads = std::max(1, n_ctx);

    #pragma omp parallel for num_threads(max_threads)
    for (int i = 0; i < N; ++i)
    {
        int tid = omp_get_thread_num();
        // Safety in case OpenMP chooses more threads than tctx_pool size
        if (tid >= n_ctx) tid = tid % n_ctx;

        TracerContext &ctx = *tctx_pool[tid];

        const ProbeRequest &req = requests[i];
        ProbeResult        &out = out_results[i];

        int elem_guess = req.elem_guess;
        ElectronTraceResult res;

        bool ok = TraceSingleCIProbe(ctx,
                                     req.r,
                                     req.z,
                                     elem_guess,
                                     res);

        out.ok         = ok;
        out.r          = req.r;
        out.z          = req.z;
        out.elem       = elem_guess;
        out.exit_code  = ok ? int(res.exit_code) : -1;
        out.is_ci      = ok && IsChargeInsensitive(res);
        out.is_cathode = ok && IsCathodeInsensitive(res);
        out.is_wall    = ok && IsWallInsensitive(res);
    }
}




// Build CivSeeds from a list of physical points (r,z)
// Any point that cannot be mapped to an element is skipped.
static CivSeeds BuildSeedsFromPhysicalPoints(mfem::ParMesh                         &mesh,
                                             const std::vector<mfem::Vector>       &points)
{
    using namespace mfem;

    CivSeeds seeds;

    const int dim = mesh.SpaceDimension();
    MFEM_VERIFY(dim == 2, "BuildSeedsFromPhysicalPoints: only 2D axisymmetric meshes supported.");

    const std::size_t n_pts = points.size();
    seeds.positions.reserve(n_pts);
    seeds.elements.reserve(n_pts);
    seeds.ips.reserve(n_pts);
    seeds.volumes.reserve(n_pts); // will typically be dummy in interface search

    for (const auto &x_phys : points)
    {
        MFEM_VERIFY(x_phys.Size() == dim,
                    "BuildSeedsFromPhysicalPoints: point dimension mismatch.");

        int elem = -1;
        IntegrationPoint ip;

        if (!FindElementForPointGlobal(mesh, x_phys, elem, ip))
        {
            // Point lies outside the mesh; skip
            continue;
        }

        seeds.positions.push_back(x_phys);
        seeds.elements.push_back(elem);
        seeds.ips.push_back(ip);

        // In the interface search we don't use per-seed volumes for CIV;
        // set a dummy value (e.g. 1.0) or 0.0 as a placeholder.
        seeds.volumes.push_back(1.0);
    }

    return seeds;
}



// Trace a batch of electrons starting from physical (r,z) positions.
// Points that cannot be mapped to an element are silently skipped;
// out_results will have the same length as the resulting seeds.
static void TraceBatchFromPhysicalPoints(const SimulationResult         &sim,
                                         const Config                   &cfg,
                                         const std::vector<mfem::Vector> &points,
                                         std::vector<ElectronTraceResult> &out_results)
{
    using namespace mfem;

    ParMesh &mesh = *sim.mesh;

    // Build seeds (this will skip points outside the mesh)
    CivSeeds seeds = BuildSeedsFromPhysicalPoints(mesh, points);

    out_results.clear();

    const std::size_t n_seeds = seeds.positions.size();
    if (n_seeds == 0)
    {
        return; // nothing to trace
    }

    out_results.resize(n_seeds);

    if (cfg.debug.debug) {std::cout << "Tracing Batch of electorns n_elems=" << n_seeds << std::endl;}

    TraceElectronFieldLines(sim, cfg, seeds, out_results);
}



static InitialColumnCIResult FindInitialColumnCI(TracerContext &base_ctx)
{
    using namespace mfem;

    InitialColumnCIResult out;
    out.found      = false;
    out.r          = 0.0;
    out.z_ci       = 0.0;
    out.elem       = -1;
    out.cathode_ci = false;
    out.wall_ci    = false;

    const TpcGeometry        &tpc    = base_ctx.tpc;
    const ElectronTraceParams &params = base_ctx.params;
    const Config &cfg = base_ctx.cfg;
    int n_threads = cfg.compute.threads.num;
    bool debug = cfg.debug.debug;

    const double r_min    = tpc.r_min;
    const double r_max    = tpc.r_max;
    const double z_min    = tpc.z_min;
    const double z_max    = tpc.z_max;
    const double geom_tol = params.geom_tol;

    const double r0 = r_max - geom_tol;
    if (r0 <= r_min + geom_tol) {
        return out;
    }

    const double z_mid = 0.5 * (z_min + z_max);

    // 1) Build a batch of candidate z’s between z_mid and z_min
    const int max_samples = std::max(4, n_threads);  // at least a few, at most a few dozen
    std::vector<ProbeRequest> reqs;
    reqs.reserve(max_samples);

    for (int k = 0; k < max_samples; ++k) {
        // Geometric-ish spacing from mid downwards
        double frac = std::pow(0.5, k); // 1, 1/2, 1/4, ...
        double z    = z_mid - (z_mid - z_min) * frac;
        if (z < z_min) z = z_min;

        ProbeRequest pr;
        pr.r         = r0;
        pr.z         = z;
        pr.elem_guess = -1;
        reqs.push_back(pr);
    }

    // 2) Build per-thread contexts as in 1.3
    int n_used_threads = std::max(1, n_threads);
    std::vector<std::unique_ptr<TracerContext>> extra;
    extra.reserve(n_used_threads - 1);

    std::vector<TracerContext*> ctx_pool(n_used_threads);
    ctx_pool[0] = &base_ctx;

    const FiniteElementCollection *fec = base_ctx.fes_phi.FEColl();
    const int                      ord = base_ctx.fes_phi.GetOrdering();

    for (int t = 1; t < n_used_threads; ++t) {
        extra.emplace_back(std::make_unique<TracerContext>(
            base_ctx.mesh, fec, ord, base_ctx.phi,
            base_ctx.adj, base_ctx.geom, base_ctx.unified_adj, cfg));
        ctx_pool[t] = extra.back().get();
    }

    // 3) Run a single parallel ProbeBatch
    std::vector<ProbeResult> results;
    ProbeBatch(reqs, ctx_pool, results);
    

    // 4) Select the highest CI z
    bool   found_ci    = false;
    double best_z      = z_min;
    int    best_elem   = -1;
    bool   saw_cathode = false;
    bool   saw_wall    = false;

    for (std::size_t i = 0; i < results.size(); ++i) {
        const ProbeResult &pr = results[i];
        if (!pr.ok || !pr.is_ci) continue;

        if (!found_ci || pr.z > best_z) {
            found_ci  = true;
            best_z    = pr.z;
            best_elem = pr.elem;
        }

        saw_cathode = saw_cathode || pr.is_cathode;
        saw_wall    = saw_wall    || pr.is_wall;
    }

    if (!found_ci) {
        return out;
    }

    // Optional: small serial refinement around best_z (reuse old bisection if you like),
    // but now centered near best_z and with reduced z-interval.

    out.found      = true;
    out.r          = r0;
    out.z_ci       = best_z;
    out.elem       = best_elem;
    out.cathode_ci = saw_cathode;
    out.wall_ci    = saw_wall;

    if (debug) {
        std::cout << "[DEBUG:MARCH] InitialColumnCI_Parallel: r=" << out.r
                  << " z_ci=" << out.z_ci
                  << " elem=" << out.elem
                  << " cath=" << out.cathode_ci
                  << " wall=" << out.wall_ci << "\n";
    }

    return out;
}




// Classify elements as charge-insensitive by comparing their centroids
// to the interface curve, then sum volumes.
static CIVolumeResult ComputeCIVFromInterface(mfem::ParMesh                 &mesh,
                                              const Config                  &cfg,
                                              const ElementGeometryData     &geom,
                                              const std::vector<InterfacePoint> &iface)
{
    using namespace mfem;

    CIVolumeResult res;
    res.V_total = 0.0;
    res.V_CI    = 0.0;

    const int ne = mesh.GetNE();
    if (ne <= 0 || iface.empty())
    {
        return res;
    }

    TpcGeometry tpc(cfg);
    const double z_check_max = cfg.tracing_params.z_max; // or a dedicated z_check_max if you add it

    for (int e = 0; e < ne; ++e)
    {
        const double r_c = geom.r_centroid[e];
        const double z_c = geom.z_centroid[e];
        const double V_e = geom.volume[e];

        if (V_e <= 0.0)
        {
            continue;
        }

        // Optionally restrict CIV computation to lower part of TPC
        if (z_c > z_check_max)
        {
            continue;
        }

        // Ignore elements outside the TPC cylinder
        if (!tpc.Inside(r_c, z_c))
        {
            continue;
        }

        res.V_total += V_e;

        // Interface-based CI decision: below or on the curve => CI
        const double z_boundary = EvaluateInterfaceZ(iface, r_c);
        if (z_c <= z_boundary)
        {
            res.V_CI += V_e;
        }
    }

    return res;
}

static ColumnRefineResult RefineColumnCI(const SimulationResult &sim,
                                         const Config           &cfg,
                                         double                  r_fixed,
                                         double                  z_lo,
                                         double                  z_hi,
                                         int                     n_samples)
{
    using namespace mfem;

    ColumnRefineResult ref;
    ref.found         = false;
    ref.z_ci          = z_lo;
    ref.is_cathode_ci = false;
    ref.is_wall_ci    = false;

    const int dim = sim.mesh->SpaceDimension();
    MFEM_VERIFY(dim == 2, "RefineColumnCI: only 2D axisymmetric meshes supported.");

    if (n_samples <= 0 || z_hi <= z_lo)
    {
        return ref;
    }

    std::vector<Vector> points;
    points.reserve(n_samples);

    const double dz = (z_hi - z_lo) / static_cast<double>(n_samples + 1);

    for (int k = 1; k <= n_samples; ++k)
    {
        double z = z_lo + k * dz;
        Vector x_phys(dim);
        x_phys[0] = r_fixed;
        x_phys[1] = z;
        points.push_back(x_phys);
    }

    std::vector<ElectronTraceResult> results;
    TraceBatchFromPhysicalPoints(sim, cfg, points, results);

    const std::size_t n_res = results.size();
    const std::size_t n_pts = points.size();
    const std::size_t n_eval = std::min(n_res, n_pts);

    double z_ci_max       = z_lo;
    bool   found_any_ci   = false;
    bool   any_cathode_ci = false;
    bool   any_wall_ci    = false;

    for (std::size_t i = 0; i < n_eval; ++i)
    {
        const ElectronTraceResult &res = results[i];
        const Vector              &x0  = points[i];
        const double z              = x0[1];

        if (IsChargeInsensitive(res))
        {
            found_any_ci = true;
            if (z > z_ci_max) { z_ci_max = z; }

            if (IsCathodeInsensitive(res)) { any_cathode_ci = true; }
            if (IsWallInsensitive(res))    { any_wall_ci    = true; }
        }
    }

    if (!found_any_ci)
    {
        return ref;
    }

    ref.found         = true;
    ref.z_ci          = z_ci_max;
    ref.is_cathode_ci = any_cathode_ci;
    ref.is_wall_ci    = any_wall_ci;

    return ref;
}


/*
Not referenced anywhere 
*/
static ColumnCISearchResult SearchColumnForCI(const SimulationResult &sim,
                                              const Config           &cfg,
                                              double                  r_fixed,
                                              double                  z_min,
                                              double                  z_max)
{
    using namespace mfem;

    ColumnCISearchResult col_res;
    col_res.found         = false;
    col_res.z_ci          = z_min;
    col_res.is_cathode_ci = false;
    col_res.is_wall_ci    = false;

    const int dim = sim.mesh->SpaceDimension();
    MFEM_VERIFY(dim == 2, "SearchColumnForCI: only 2D axisymmetric meshes supported.");

    const double z_mid = 0.5 * (z_min + z_max);

    // Build a sequence of test heights: z_mid, then halfway towards z_min,
    // then halfway again, etc. We keep them all and trace in one batch.
    const int max_levels = 8; // enough to get close to cathode

    std::vector<double> z_tests;
    z_tests.reserve(max_levels);

    double z_hi = z_mid;
    double z_lo = z_min;

    // z0 = mid, then mid - (mid - z_min)/2, etc
    for (int k = 0; k < max_levels; ++k)
    {
        double zk;
        if (k == 0)
        {
            zk = z_mid;
        }
        else
        {
            double frac = std::pow(0.5, k); // 1/2, 1/4, 1/8, ...
            zk = z_mid - (z_mid - z_min) * frac;
        }
        if (zk < z_min) { zk = z_min; }
        if (zk > z_max) { zk = z_max; }

        // Avoid duplicates if rounding collapses
        if (!z_tests.empty() && std::fabs(zk - z_tests.back()) < 1e-12)
        {
            continue;
        }

        z_tests.push_back(zk);
    }

    if (z_tests.empty())
    {
        return col_res;
    }

    // Build physical points along this column
    std::vector<Vector> points;
    points.reserve(z_tests.size());

    for (double z : z_tests)
    {
        Vector x_phys(dim);
        x_phys[0] = r_fixed;
        x_phys[1] = z;
        points.push_back(x_phys);
    }

    // Trace all in one batch
    std::vector<ElectronTraceResult> results;
    TraceBatchFromPhysicalPoints(sim, cfg, points, results);

    // Note: some points may have been skipped by BuildSeedsFromPhysicalPoints
    // if they were outside the mesh. So we need to be careful: we assume
    // that if a point could not be mapped, it simply won't appear in results.
    // For a well-chosen r_fixed,z-range this should be rare; if it happens,
    // we effectively ignore that z in the search.

    // Walk through results to find highest z that is CI
    double z_ci_max         = z_min;
    bool   found_any_ci     = false;
    bool   any_cathode_ci   = false;
    bool   any_wall_ci      = false;

    // We assume results are in the same order as 'points'
    const std::size_t n_res = results.size();
    const std::size_t n_pts = points.size();
    const std::size_t n_eval = std::min(n_res, n_pts);

    for (std::size_t i = 0; i < n_eval; ++i)
    {
        const ElectronTraceResult &res = results[i];
        const Vector              &x0  = points[i];
        const double z              = x0[1];

        if (IsChargeInsensitive(res))
        {
            found_any_ci = true;
            if (z > z_ci_max) { z_ci_max = z; }

            if (IsCathodeInsensitive(res)) { any_cathode_ci = true; }
            if (IsWallInsensitive(res))    { any_wall_ci    = true; }
        }
    }

    if (!found_any_ci)
    {
        return col_res; // no CI in this column
    }

    col_res.found         = true;
    col_res.z_ci          = z_ci_max;
    col_res.is_cathode_ci = any_cathode_ci;
    col_res.is_wall_ci    = any_wall_ci;

    return col_res;
}

/*
Main Logic fails on nearest neighbor search, I think the base assumption here is wrong 

Best bet 
GPU offload start at the bottom use bracketing logic use cathode periodicity as base
*/
static std::vector<InterfacePoint>
BuildInterfaceCurve_Marching(TracerContext &ctx)
{
    using namespace mfem;

    ParMesh                   &mesh    = ctx.mesh;
    const ElementAdjacency    &adj     = ctx.adj;
    const ElementGeometryData &geom    = ctx.geom;
    const TpcGeometry         &tpc     = ctx.tpc;
    const ElectronTraceParams &params  = ctx.params;
    const auto                &unified_adj = ctx.unified_adj; // <--- use unified adjacency
    bool debug                         = ctx.cfg.debug.debug;

    if (debug)
    {
        std::cout << "[DEBUG:MARCH] BuildInterfaceCurve_Marching\n";
        std::cout << "[DEBUG:MARCH] Mesh NE=" << mesh.GetNE() << "\n";
        std::cout << "[DEBUG:MARCH] r_min=" << tpc.r_min
                  << " r_max=" << tpc.r_max
                  << " z_min=" << tpc.z_min
                  << " z_max=" << tpc.z_max
                  << " geom_tol=" << params.geom_tol << "\n";
    }

    std::vector<InterfacePoint> iface;

    const double r_min    = tpc.r_min;
    const double r_max    = tpc.r_max;
    const double z_min    = tpc.z_min;
    const double z_max    = tpc.z_max;
    const double geom_tol = params.geom_tol;

    // ---------------------------------------------------------------------
    // 1) Initial column CI (at r ~ r_max)
    // ---------------------------------------------------------------------
    InitialColumnCIResult init = FindInitialColumnCI(ctx);

    if (debug)
    {
        std::cout << "[DEBUG:MARCH] InitialColumnCI: found=" << init.found
                  << " r=" << init.r
                  << " z_ci=" << init.z_ci
                  << " elem=" << init.elem
                  << " cathode_ci=" << init.cathode_ci
                  << " wall_ci=" << init.wall_ci << "\n";
    }

    if (!init.found || init.elem < 0)
    {
        if (debug)
            std::cout << "[DEBUG:MARCH] No initial CI found. Exiting.\n";
        return iface;
    }

    // First interface point from initial column search
    iface.push_back({init.r, init.z_ci, init.elem,
                     init.cathode_ci, init.wall_ci});

    // Marching state
    int    elem_i = init.elem;
    double r_i    = init.r;
    double z_i    = init.z_ci;

    double smallest_r_reached = r_i;

    // Cathode detection state
    bool   cathode_seen   = init.cathode_ci;
    double z_first_cath   = init.cathode_ci ? init.z_ci
                                            : std::numeric_limits<double>::infinity();

    const int   max_steps = 200;
    const double eps_r    = 1e-9;
    const double eps_z    = 1e-9;

    // ---------------------------------------------------------------------
    // 2) Main adjacency-based march: move to CI neighbors that are down-left.
    //    Use unified adjacency (face + vertex neighbors).
    //    Stop as soon as we see the first cathode hit.
    // ---------------------------------------------------------------------
    if (!cathode_seen) // if initial column already saw cathode, skip marching
    {
        for (int step = 0; step < max_steps; ++step)
        {
            if (debug)
            {
                std::cout << "\n[DEBUG:MARCH] --- Step " << step
                          << " at elem=" << elem_i
                          << " r_i=" << r_i
                          << " z_i=" << z_i << "\n";
            }

            // Stop near inner radius for wall-branch
            if (r_i <= r_min + geom_tol)
            {
                if (debug)
                    std::cout << "[DEBUG:MARCH] Reached inner radius stop (wall marching).\n";
                break;
            }

            // Use unified adjacency instead of face-only neighbors
            const auto &nbrs = unified_adj[elem_i];

            if (debug)
            {
                std::cout << "[DEBUG:MARCH] Unified neighbors of elem " << elem_i
                          << ": " << nbrs.size() << "\n";
            }

            bool   found_ci_here   = false;
            double best_z_ci       = -std::numeric_limits<double>::infinity();
            int    best_elem       = -1;
            bool   cathode_here    = false;
            bool   wall_here       = false;

            bool   trigger_cathode = false; // becomes true on first cathode hit

            // Scan neighbors; pick those whose centroids lie "down-left"
            for (int nb : nbrs)
            {
                if (nb < 0 || nb >= mesh.GetNE())
                    continue;

                const double r_n = geom.r_centroid[nb];
                const double z_n = geom.z_centroid[nb];

                // down-left filter (allow tiny tolerance)
                if (r_n > r_i + eps_r) continue;
                if (z_n > z_i + eps_z) continue;

                if (debug)
                {
                    std::cout << "[DEBUG:MARCH]   neighbor elem=" << nb
                              << " r_n=" << r_n
                              << " z_n=" << z_n << " (candidate)\n";
                }

                double rr = r_n;
                double zz = z_n;

                ElectronTraceResult res;
                int elem_guess = nb;

                bool ok = TraceSingleCIProbe(ctx, rr, zz, elem_guess, res);

                if (debug)
                {
                    std::cout << "[DEBUG:MARCH]      probe ok=" << ok
                              << " exit=" << int(res.exit_code)
                              << " cathode=" << IsCathodeInsensitive(res)
                              << " wall=" << IsWallInsensitive(res)
                              << "\n";
                }

                if (!ok) continue;
                if (!IsChargeInsensitive(res)) continue;

                // First time we see a cathode hit anywhere → trigger cathode phase
                if (IsCathodeInsensitive(res) && !cathode_seen)
                {
                    cathode_seen   = true;
                    z_first_cath   = z_n;
                    trigger_cathode = true;

                    if (debug)
                    {
                        std::cout << "[DEBUG:MARCH] First cathode hit detected at "
                                  << "r=" << r_n << " z=" << z_n
                                  << " → stopping wall marching and entering cathode phase.\n";
                    }
                    break; // break neighbor loop
                }

                // Wall-based CI candidate
                // Keep highest z among CI neighbors
                if (!found_ci_here || z_n > best_z_ci)
                {
                    found_ci_here = true;
                    best_z_ci     = z_n;
                    best_elem     = elem_guess;
                    cathode_here  = IsCathodeInsensitive(res);
                    wall_here     = IsWallInsensitive(res);
                }
                else
                {
                    // OR-aggregate CI type flags
                    cathode_here = cathode_here || IsCathodeInsensitive(res);
                    wall_here    = wall_here    || IsWallInsensitive(res);
                }
            } // neighbor loop

            if (trigger_cathode)
            {
                // We just detected the first cathode hit; stop wall marching.
                break;
            }

            if (!found_ci_here || best_elem < 0)
            {
                if (debug)
                    std::cout << "[DEBUG:MARCH] No CI neighbor found → terminating wall march\n";
                break;
            }

            // Accept new wall interface point
            elem_i = best_elem;
            r_i    = geom.r_centroid[best_elem];
            z_i    = best_z_ci;

            iface.push_back({r_i, z_i, elem_i, cathode_here, wall_here});

            smallest_r_reached = std::min(smallest_r_reached, r_i);
        } // wall-marching loop
    } // if (!cathode_seen)

    // ---------------------------------------------------------------------
    // 3) Cathode plane construction (triggered as soon as a cathode hit is seen)
    // ---------------------------------------------------------------------
    if (cathode_seen)
    {
        if (debug)
        {
            std::cout << "[DEBUG:MARCH] Entering cathode plane construction.\n";
            std::cout << "[DEBUG:MARCH] z_first_cath (initial guess)=" << z_first_cath << "\n";
        }

        double z_cath = z_first_cath;
        if (!std::isfinite(z_cath))
        {
            // Fallback: if for some reason we didn't record z_first_cath,
            // take the last interface z or fall back to z_min.
            z_cath = (iface.empty() ? z_min : iface.back().z_ci);
        }

        double r_plane_start = smallest_r_reached;

        if (!iface.empty())
        {
            // Find index of smallest-r point
            int idx_min = 0;
            double r_min_iface = iface[0].r;
            for (int i = 1; i < (int)iface.size(); ++i)
            {
                if (iface[i].r < r_min_iface)
                {
                    r_min_iface = iface[i].r;
                    idx_min     = i;
                }
            }

            // Update that point to lie on the cathode plane
            iface[idx_min].z_ci       = z_cath;
            iface[idx_min].cathode_ci = true;
            iface[idx_min].wall_ci    = false;

            // Add an additional plane point at inner radius (if distinct in r)
            const double r_inner = r_min + geom_tol;
            if (r_inner < r_min_iface - 1e-12)
            {
                InterfacePoint p;
                p.r          = r_inner;
                p.z_ci       = z_cath;
                p.elem       = iface[idx_min].elem;
                p.cathode_ci = true;
                p.wall_ci    = false;
                iface.push_back(p);
            }
        }
        else
        {
            // If somehow the interface has no points, add just a cathode plane segment
            const double r_inner = r_min + geom_tol;
            InterfacePoint p1{r_plane_start, z_cath, -1, true, false};
            InterfacePoint p2{r_inner,       z_cath, -1, true, false};
            iface.push_back(p1);
            iface.push_back(p2);
        }
    }

    // ---------------------------------------------------------------------
    // 4) Sort by radius and print summary
    // ---------------------------------------------------------------------
    std::sort(iface.begin(), iface.end(),
              [](const InterfacePoint &a, const InterfacePoint &b)
              {
                  return a.r < b.r;
              });

    if (debug)
    {
        std::cout << "[DEBUG:MARCH] Completed marching. Points = "
                  << iface.size() << "\n";
        for (const auto &p : iface)
        {
            std::cout << "  r=" << p.r
                      << " z_ci=" << p.z_ci
                      << " elem=" << p.elem
                      << " cath=" << p.cathode_ci
                      << " wall=" << p.wall_ci << "\n";
        }
    }

    return iface;
}



// Main Wrapper
double ComputeCIV_Marching(const SimulationResult &sim,
                           const Config           &cfg)
{
    using namespace mfem;

    MFEM_VERIFY(sim.mesh, "ComputeCIV_Marching: mesh is null");
    ParMesh &global_mesh = *sim.mesh;

    MFEM_VERIFY(sim.V, "ComputeCIV_Marching: potential V is null");
    GridFunction &global_phi = *sim.V;

    const FiniteElementSpace      *fes_phi = global_phi.FESpace();
    const FiniteElementCollection *fec     = fes_phi->FEColl();
    const int                      ordering = fes_phi->GetOrdering();

    // 1) Build shared geometry/connectivity ONCE
    ElementAdjacency           adj         = BuildAdjacency(global_mesh);
    ElementGeometryData        geom        = ComputeElementGeometry(global_mesh);
    std::vector<std::vector<int>> vertex_adj = BuildVertexAdjacency(global_mesh);
    std::vector<std::vector<int>> unified_adj =
        BuildUnifiedAdjacency(adj, vertex_adj, global_mesh.GetNE());

    // 2) Build a single base TracerContext
    TracerContext base_ctx(global_mesh,
                           fec,
                           ordering,
                           global_phi,
                           adj,
                           geom,
                           unified_adj,
                           cfg);

    // 3) Build interface using base_ctx (it will internally clone per-thread)
    std::vector<InterfacePoint> iface =
        BuildInterfaceCurve_Marching(base_ctx);

    if (iface.empty()) {
        if (cfg.debug.debug)
            std::cout << "[DEBUG:CIV] No interface points found; CIV = 0.\n";
        return 0.0;
    }

    std::sort(iface.begin(), iface.end(),
              [](const InterfacePoint &a, const InterfacePoint &b)
              { return a.r < b.r; });

    if (cfg.tracing_params.dump_civ_boundary) {
        namespace fs = std::filesystem;
        fs::path outdir(cfg.save_path);
        fs::path fname = outdir / "civ_interface_boundary.csv";
        std::ofstream ofs(fname);
        ofs << "# r, z_ci\n";
        for (const auto &ip : iface)
            ofs << ip.r << "," << ip.z_ci << "\n";
        if (cfg.debug.debug)
            std::cout << "[DEBUG:CIV] Dumped CIV boundary to: "
                      << fname.string() << "\n";
    }

    // 4) Volume classification uses the same geom we already have
    CIVolumeResult vres = ComputeCIVFromInterface(global_mesh, cfg, geom, iface);
    double V_total = vres.V_total;
    double V_CI    = vres.V_CI;
    double frac    = (V_total > 0.0) ? (V_CI / V_total) : 0.0;

    if (cfg.debug.debug) {
        std::cout << "\n----- Charge Insensitive Volume (Marching) -----\n";
        std::cout << "Total Volume considered: " << V_total << "\n";
        std::cout << "CIV Volume:              " << V_CI    << "\n";
        std::cout << "CIV Fraction:            " << frac    << "\n";
        std::cout << "------------------------------------------------\n";
    }

    return frac;
}