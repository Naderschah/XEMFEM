#include "interpolator.h"
#include <cmath>
#include <cstdint>
#include <utility>

namespace
{
static void ComputeGlobalBoundingBox(mfem::ParMesh &pmesh,
                                     mfem::Vector &bbmin3,
                                     mfem::Vector &bbmax3)
{
    const int space_dim = pmesh.SpaceDimension();
    mfem::Vector bbmin_s(space_dim), bbmax_s(space_dim);
    pmesh.GetBoundingBox(bbmin_s, bbmax_s);

    double lmin[3] = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()
    };
    double lmax[3] = {
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    };

    for (int d = 0; d < space_dim; ++d)
    {
        lmin[d] = bbmin_s[d];
        lmax[d] = bbmax_s[d];
    }

    double gmin[3] = {0.0, 0.0, 0.0};
    double gmax[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(lmin, gmin, 3, MPI_DOUBLE, MPI_MIN, pmesh.GetComm());
    MPI_Allreduce(lmax, gmax, 3, MPI_DOUBLE, MPI_MAX, pmesh.GetComm());

    bbmin3.SetSize(3);
    bbmax3.SetSize(3);
    for (int d = 0; d < 3; ++d)
    {
        bbmin3[d] = gmin[d];
        bbmax3[d] = gmax[d];
    }
}

static void ApplyInterpolationBoundsFromConfig(const Config &cfg,
                                               const int space_dim,
                                               mfem::Vector &bbmin3,
                                               mfem::Vector &bbmax3)
{
    auto apply_axis = [&](const char *axis_name,
                          const int axis,
                          const double vmin,
                          const double vmax)
    {
        const bool has_min = std::isfinite(vmin);
        const bool has_max = std::isfinite(vmax);
        MFEM_VERIFY(has_min == has_max,
                    std::string("interpolate.") + axis_name +
                    "_min and interpolate." + axis_name +
                    "_max must be set together.");
        if (has_min)
        {
            MFEM_VERIFY(vmax > vmin,
                        std::string("interpolate.") + axis_name +
                        "_max must be greater than interpolate." +
                        axis_name + "_min.");
            bbmin3[axis] = vmin;
            bbmax3[axis] = vmax;
        }
    };

    apply_axis("x", 0, cfg.interp.x_min, cfg.interp.x_max);
    apply_axis("y", 1, cfg.interp.y_min, cfg.interp.y_max);
    if (space_dim == 3)
    {
        apply_axis("z", 2, cfg.interp.z_min, cfg.interp.z_max);
    }
}

static inline double CenteredGridSpacing(const double xmin,
                                         const double xmax,
                                         const int npts)
{
    MFEM_VERIFY(npts > 0, "Grid point count must be > 0.");
    return (xmax - xmin) / static_cast<double>(npts);
}

static inline double CenteredGridCoordinate(const double xmin,
                                            const double spacing,
                                            const int idx)
{
    return xmin + (static_cast<double>(idx) + 0.5) * spacing;
}

static mfem::IntegrationPoint MakeIntegrationPoint(const double *ref, const int sdim)
{
    mfem::IntegrationPoint ip;
    if (sdim == 2) { ip.Set2(ref[0], ref[1]); }
    else           { ip.Set3(ref[0], ref[1], ref[2]); }
    return ip;
}

static inline bool IsLocalElementId(const mfem::ParMesh &mesh, const int elem_id)
{
    return elem_id >= 0 && elem_id < mesh.GetNE();
}

static bool IsInsideElementReferencePoint(mfem::ElementTransformation &T,
                                          const mfem::IntegrationPoint &ip,
                                          const int sdim)
{
    mfem::Vector xphys(sdim);
    T.Transform(ip, xphys);
    mfem::InverseElementTransformation inv(&T);
    mfem::IntegrationPoint ip_back;
    const int status = inv.Transform(xphys, ip_back);
    return status == mfem::InverseElementTransformation::Inside;
}

static bool NudgeReferencePointTowardElementCenter(mfem::ParMesh &mesh,
                                                   const int elem_id,
                                                   const double *ref_in,
                                                   const int sdim,
                                                   mfem::IntegrationPoint &ip_out)
{
    if (!IsLocalElementId(mesh, elem_id)) { return false; }
    mfem::ElementTransformation *T = mesh.GetElementTransformation(elem_id);
    if (T == nullptr) { return false; }

    const mfem::IntegrationPoint ip_ref = MakeIntegrationPoint(ref_in, sdim);
    if (IsInsideElementReferencePoint(*T, ip_ref, sdim))
    {
        ip_out = ip_ref;
        return true;
    }

    const mfem::Geometry::Type geom = mesh.GetElementBaseGeometry(elem_id);
    const mfem::IntegrationPoint &ctr = mfem::Geometries.GetCenter(geom);

    double alpha = 1.0e-8;
    for (int it = 0; it < 14; ++it)
    {
        const double x = (1.0 - alpha) * ip_ref.x + alpha * ctr.x;
        const double y = (1.0 - alpha) * ip_ref.y + alpha * ctr.y;
        const double z = (1.0 - alpha) * ip_ref.z + alpha * ctr.z;

        mfem::IntegrationPoint ip_try;
        if (sdim == 2) { ip_try.Set2(x, y); }
        else           { ip_try.Set3(x, y, z); }

        if (IsInsideElementReferencePoint(*T, ip_try, sdim))
        {
            ip_out = ip_try;
            return true;
        }

        alpha = std::min(1.0, alpha * 4.0);
    }

    // Final fallback: element center should be strictly interior for standard geometries.
    if (sdim == 2) { ip_out.Set2(ctr.x, ctr.y); }
    else           { ip_out.Set3(ctr.x, ctr.y, ctr.z); }
    return IsInsideElementReferencePoint(*T, ip_out, sdim);
}

static bool BuildBoundaryNudgedPhysicalPair(mfem::ParMesh &mesh,
                                            const int elem_id,
                                            const double *ref_in,
                                            const int sdim,
                                            const double eps,
                                            double *x_in,
                                            double *x_out)
{
    if (!IsLocalElementId(mesh, elem_id)) { return false; }
    mfem::ElementTransformation *T = mesh.GetElementTransformation(elem_id);
    if (T == nullptr) { return false; }

    const mfem::IntegrationPoint ip_ref = MakeIntegrationPoint(ref_in, sdim);
    const mfem::Geometry::Type geom = mesh.GetElementBaseGeometry(elem_id);
    const mfem::IntegrationPoint &ctr = mfem::Geometries.GetCenter(geom);
    mfem::IntegrationPoint ip_ctr;
    if (sdim == 2) { ip_ctr.Set2(ctr.x, ctr.y); }
    else           { ip_ctr.Set3(ctr.x, ctr.y, ctr.z); }

    mfem::Vector x0(sdim), xc(sdim);
    T->Transform(ip_ref, x0);
    T->Transform(ip_ctr, xc);

    double n2 = 0.0;
    for (int d = 0; d < sdim; ++d)
    {
        const double dd = x0[d] - xc[d];
        n2 += dd * dd;
    }
    if (!(n2 > 0.0) || !std::isfinite(n2)) { return false; }
    const double inv_n = 1.0 / std::sqrt(n2);

    for (int d = 0; d < sdim; ++d)
    {
        const double dir = (x0[d] - xc[d]) * inv_n;
        x_in[d]  = x0[d] - eps * dir; // toward element center side
        x_out[d] = x0[d] + eps * dir; // opposite side
    }
    return true;
}

static void SetCode1SamplesToNaN(const mfem::Array<unsigned int> &code_initial,
                                 const int nloc,
                                 const int ncomp,
                                 mfem::Vector &vals_flat)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int lp = 0; lp < nloc; ++lp)
    {
        if (code_initial[lp] != 1) { continue; }
        for (int c = 0; c < ncomp; ++c)
        {
            vals_flat[lp * ncomp + c] = nan;
        }
    }
}

template <typename EvalFn>
static void AverageCode1BoundarySamples(mfem::ParMesh &mesh,
                                        const mfem::Vector &spacing,
                                        const mfem::Array<unsigned int> &code_initial,
                                        const mfem::Array<unsigned int> &elem_initial,
                                        const mfem::Vector &ref_initial,
                                        const int nloc,
                                        const int sdim,
                                        const int ncomp,
                                        EvalFn &&eval_fn,
                                        mfem::Vector &vals_flat)
{
    std::vector<int> bnd_lp;
    bnd_lp.reserve(static_cast<std::size_t>(nloc));
    for (int lp = 0; lp < nloc; ++lp)
    {
        if (code_initial[lp] == 1) { bnd_lp.push_back(lp); }
    }
    const int nb = static_cast<int>(bnd_lp.size());
    int nb_global = 0;
    MPI_Allreduce(&nb, &nb_global, 1, MPI_INT, MPI_SUM, mesh.GetComm());
    if (nb_global == 0) { return; }

    double min_h = std::numeric_limits<double>::infinity();
    for (int d = 0; d < sdim; ++d)
    {
        if (spacing[d] > 0.0) { min_h = std::min(min_h, spacing[d]); }
    }
    const double eps = std::isfinite(min_h) ? std::max(1.0e-12, 1.0e-3 * min_h) : 1.0e-9;

    std::vector<unsigned char> pair_ok(static_cast<std::size_t>(nb), static_cast<unsigned char>(0));
    std::vector<int> bnd_to_query(static_cast<std::size_t>(nb), -1);
    std::vector<double> x_in_valid;
    std::vector<double> x_out_valid;
    x_in_valid.reserve(static_cast<std::size_t>(nb * sdim));
    x_out_valid.reserve(static_cast<std::size_t>(nb * sdim));
    int nvalid = 0;

    for (int j = 0; j < nb; ++j)
    {
        const int lp = bnd_lp[static_cast<std::size_t>(j)];
        const int elem_id = static_cast<int>(elem_initial[lp]);
        const double *ref = ref_initial.GetData() + lp * sdim;
        double xin[3] = {0.0, 0.0, 0.0};
        double xout[3] = {0.0, 0.0, 0.0};

        if (BuildBoundaryNudgedPhysicalPair(mesh, elem_id, ref, sdim, eps, xin, xout))
        {
            pair_ok[static_cast<std::size_t>(j)] = 1;
            bnd_to_query[static_cast<std::size_t>(j)] = nvalid++;
            for (int d = 0; d < sdim; ++d)
            {
                x_in_valid.push_back(xin[d]);
                x_out_valid.push_back(xout[d]);
            }
        }
    }

    int nvalid_global = 0;
    MPI_Allreduce(&nvalid, &nvalid_global, 1, MPI_INT, MPI_SUM, mesh.GetComm());
    if (nvalid_global == 0)
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (int j = 0; j < nb; ++j)
        {
            const int lp = bnd_lp[static_cast<std::size_t>(j)];
            for (int c = 0; c < ncomp; ++c)
            {
                vals_flat[lp * ncomp + c] = nan;
            }
        }
        return;
    }

    // eval_fn uses FindPointsGSLIB collectives. All ranks must enter even when
    // this rank has no valid boundary pairs to process.
    const int nquery = (nvalid > 0) ? nvalid : 1;

    mfem::Vector x_in(nquery * sdim);
    mfem::Vector x_out(nquery * sdim);
    x_in = 0.0;
    x_out = 0.0;
    if (nvalid > 0)
    {
        MFEM_VERIFY(static_cast<int>(x_in_valid.size()) == nvalid * sdim,
                    "x_in_valid size mismatch.");
        MFEM_VERIFY(static_cast<int>(x_out_valid.size()) == nvalid * sdim,
                    "x_out_valid size mismatch.");
        for (int i = 0; i < nvalid * sdim; ++i)
        {
            x_in[i] = x_in_valid[static_cast<std::size_t>(i)];
            x_out[i] = x_out_valid[static_cast<std::size_t>(i)];
        }
    }

    mfem::Vector vals_in(nquery * ncomp), vals_out(nquery * ncomp);
    eval_fn(x_in, nquery, vals_in);
    eval_fn(x_out, nquery, vals_out);

    for (int j = 0; j < nb; ++j)
    {
        const int lp = bnd_lp[static_cast<std::size_t>(j)];
        const bool ok_pair = (pair_ok[static_cast<std::size_t>(j)] != 0);
        const int q = ok_pair ? bnd_to_query[static_cast<std::size_t>(j)] : -1;
        bool ok_in = ok_pair;
        bool ok_out = ok_pair;

        for (int c = 0; c < ncomp; ++c)
        {
            if (!ok_pair) { break; }
            ok_in = ok_in && std::isfinite(vals_in[q * ncomp + c]);
            ok_out = ok_out && std::isfinite(vals_out[q * ncomp + c]);
        }

        for (int c = 0; c < ncomp; ++c)
        {
            const double vin = ok_pair ? vals_in[q * ncomp + c]
                                       : std::numeric_limits<double>::quiet_NaN();
            const double vout = ok_pair ? vals_out[q * ncomp + c]
                                        : std::numeric_limits<double>::quiet_NaN();
            if (ok_in && ok_out)      { vals_flat[lp * ncomp + c] = 0.5 * (vin + vout); }
            else if (ok_in)           { vals_flat[lp * ncomp + c] = vin; }
            else if (ok_out)          { vals_flat[lp * ncomp + c] = vout; }
            else                      { vals_flat[lp * ncomp + c] = std::numeric_limits<double>::quiet_NaN(); }
        }
    }
}

template <typename EvalFn>
static void ResolveCode1Samples(mfem::ParMesh &mesh,
                                const mfem::Vector &spacing,
                                const InterpolationCode1Mode mode,
                                const mfem::Array<unsigned int> &code_initial,
                                const mfem::Array<unsigned int> &elem_initial,
                                const mfem::Vector &ref_initial,
                                const int nloc,
                                const int sdim,
                                const int ncomp,
                                EvalFn &&eval_fn,
                                mfem::Vector &vals_flat)
{
    if (mode == InterpolationCode1Mode::NaN)
    {
        SetCode1SamplesToNaN(code_initial, nloc, ncomp, vals_flat);
        return;
    }
    if (mode == InterpolationCode1Mode::AverageElements)
    {
        AverageCode1BoundarySamples(mesh, spacing, code_initial, elem_initial, ref_initial,
                                    nloc, sdim, ncomp, std::forward<EvalFn>(eval_fn), vals_flat);
        return;
    }
    // AcceptValue: keep the value returned by the initial finder/evaluator call.
}

// Duplicated from MPI tracer to keep field sampling behavior aligned.
class MPITracerStyleFieldEvaluator
{
public:
    MPITracerStyleFieldEvaluator(mfem::ParMesh& mesh,
                                 const mfem::ParGridFunction& phi,
                                 mfem::FindPointsGSLIB& finder,
                                 int sdim)
        : mesh_(mesh), phi_(phi), finder_(finder), sdim_(sdim)
    {
    }

    void EvaluateE(const mfem::Vector& points, int n_points, mfem::Vector& E)
    {
        MFEM_VERIFY(points.Size() == n_points * sdim_,
                    "Point coordinate buffer size does not match n_points * sdim.");

        finder_.FindPoints(points, mfem::Ordering::byVDIM);

        finder_.DistributePointInfoToOwningMPIRanks(recv_elem_, recv_ref_, recv_code_);
        const int n_recv = recv_elem_.Size();
        MFEM_VERIFY(recv_ref_.Size() == n_recv * sdim_,
                    "Invalid reference coordinates size.");

        recv_E_.SetSize(n_recv * sdim_);
        mfem::Vector grad(sdim_);

        for (int i = 0; i < n_recv; ++i)
        {
            if (recv_code_[i] == 2) // point not found
            {
                for (int d = 0; d < sdim_; ++d)
                {
                    recv_E_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
                }
                continue;
            }

            const int elem_id = static_cast<int>(recv_elem_[i]);
            if (!IsLocalElementId(mesh_, elem_id))
            {
                for (int d = 0; d < sdim_; ++d)
                {
                    recv_E_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
                }
                continue;
            }
            mfem::ElementTransformation* T = mesh_.GetElementTransformation(elem_id);
            MFEM_VERIFY(T != nullptr, "Null element transformation.");

            mfem::IntegrationPoint ip;
            const double* ref = recv_ref_.GetData() + i * sdim_;

            if (recv_code_[i] == 1) // boundary-projected: nudge inward until inside
            {
                if (!NudgeReferencePointTowardElementCenter(mesh_, elem_id, ref, sdim_, ip))
                {
                    for (int d = 0; d < sdim_; ++d)
                    {
                        recv_E_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
                    }
                    continue;
                }
            }
            else
            {
                ip = MakeIntegrationPoint(ref, sdim_);
            }

            T->SetIntPoint(&ip);
            phi_.GetGradient(*T, grad);
            grad *= -1.0; // E = -grad(phi)

            for (int d = 0; d < sdim_; ++d)
            {
                recv_E_[i * sdim_ + d] = grad[d];
            }
        }

        finder_.DistributeInterpolatedValues(recv_E_, sdim_,
                                             mfem::Ordering::byVDIM, E);
    }

private:
    mfem::ParMesh& mesh_;
    const mfem::ParGridFunction& phi_;
    mfem::FindPointsGSLIB& finder_;
    const int sdim_;

    mfem::Array<unsigned int> recv_elem_;
    mfem::Vector recv_ref_;
    mfem::Array<unsigned int> recv_code_;
    mfem::Vector recv_E_;
};

class MPITracerStylePotentialEvaluator
{
public:
    MPITracerStylePotentialEvaluator(mfem::ParMesh& mesh,
                                     const mfem::ParGridFunction& phi,
                                     mfem::FindPointsGSLIB& finder,
                                     int sdim)
        : mesh_(mesh), phi_(phi), finder_(finder), sdim_(sdim)
    {
    }

    void EvaluateV(const mfem::Vector& points, int n_points, mfem::Vector& Vvals)
    {
        MFEM_VERIFY(points.Size() == n_points * sdim_,
                    "Point coordinate buffer size does not match n_points * sdim.");

        finder_.FindPoints(points, mfem::Ordering::byVDIM);

        finder_.DistributePointInfoToOwningMPIRanks(recv_elem_, recv_ref_, recv_code_);
        const int n_recv = recv_elem_.Size();
        MFEM_VERIFY(recv_ref_.Size() == n_recv * sdim_,
                    "Invalid reference coordinates size.");

        recv_V_.SetSize(n_recv);

        for (int i = 0; i < n_recv; ++i)
        {
            if (recv_code_[i] == 2) // point not found
            {
                recv_V_[i] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            const int elem_id = static_cast<int>(recv_elem_[i]);
            if (!IsLocalElementId(mesh_, elem_id))
            {
                recv_V_[i] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            mfem::ElementTransformation* T = mesh_.GetElementTransformation(elem_id);
            MFEM_VERIFY(T != nullptr, "Null element transformation.");

            mfem::IntegrationPoint ip;
            const double* ref = recv_ref_.GetData() + i * sdim_;

            if (recv_code_[i] == 1) // boundary-projected: nudge inward until inside
            {
                if (!NudgeReferencePointTowardElementCenter(mesh_, elem_id, ref, sdim_, ip))
                {
                    recv_V_[i] = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }
            }
            else
            {
                ip = MakeIntegrationPoint(ref, sdim_);
            }

            recv_V_[i] = phi_.GetValue(elem_id, ip);
        }

        finder_.DistributeInterpolatedValues(recv_V_, 1, mfem::Ordering::byVDIM, Vvals);
    }

private:
    mfem::ParMesh& mesh_;
    const mfem::ParGridFunction& phi_;
    mfem::FindPointsGSLIB& finder_;
    const int sdim_;

    mfem::Array<unsigned int> recv_elem_;
    mfem::Vector recv_ref_;
    mfem::Array<unsigned int> recv_code_;
    mfem::Vector recv_V_;
};

class MPITracerStyleVectorFieldEvaluator
{
public:
    MPITracerStyleVectorFieldEvaluator(mfem::ParMesh& mesh,
                                       const mfem::ParGridFunction& vec_field,
                                       mfem::FindPointsGSLIB& finder,
                                       int sdim)
        : mesh_(mesh), vec_field_(vec_field), finder_(finder), sdim_(sdim)
    {
    }

    void EvaluateVectorField(const mfem::Vector& points, int n_points, mfem::Vector& vec_vals)
    {
        MFEM_VERIFY(points.Size() == n_points * sdim_,
                    "Point coordinate buffer size does not match n_points * sdim.");

        finder_.FindPoints(points, mfem::Ordering::byVDIM);

        finder_.DistributePointInfoToOwningMPIRanks(recv_elem_, recv_ref_, recv_code_);
        const int n_recv = recv_elem_.Size();
        MFEM_VERIFY(recv_ref_.Size() == n_recv * sdim_,
                    "Invalid reference coordinates size.");

        recv_vec_.SetSize(n_recv * sdim_);
        mfem::Vector v(sdim_);

        for (int i = 0; i < n_recv; ++i)
        {
            if (recv_code_[i] == 2) // point not found
            {
                for (int d = 0; d < sdim_; ++d)
                {
                    recv_vec_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
                }
                continue;
            }

            const int elem_id = static_cast<int>(recv_elem_[i]);
            if (!IsLocalElementId(mesh_, elem_id))
            {
                for (int d = 0; d < sdim_; ++d)
                {
                    recv_vec_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
                }
                continue;
            }
            mfem::ElementTransformation* T = mesh_.GetElementTransformation(elem_id);
            MFEM_VERIFY(T != nullptr, "Null element transformation.");

            mfem::IntegrationPoint ip;
            const double* ref = recv_ref_.GetData() + i * sdim_;

            if (recv_code_[i] == 1) // boundary-projected: nudge inward until inside
            {
                if (!NudgeReferencePointTowardElementCenter(mesh_, elem_id, ref, sdim_, ip))
                {
                    for (int d = 0; d < sdim_; ++d)
                    {
                        recv_vec_[i * sdim_ + d] = std::numeric_limits<double>::quiet_NaN();
                    }
                    continue;
                }
            }
            else
            {
                ip = MakeIntegrationPoint(ref, sdim_);
            }

            T->SetIntPoint(&ip);
            vec_field_.GetVectorValue(*T, ip, v);

            for (int d = 0; d < sdim_; ++d)
            {
                recv_vec_[i * sdim_ + d] = v[d];
            }
        }

        finder_.DistributeInterpolatedValues(recv_vec_, sdim_,
                                             mfem::Ordering::byVDIM, vec_vals);
    }

private:
    mfem::ParMesh& mesh_;
    const mfem::ParGridFunction& vec_field_;
    mfem::FindPointsGSLIB& finder_;
    const int sdim_;

    mfem::Array<unsigned int> recv_elem_;
    mfem::Vector recv_ref_;
    mfem::Array<unsigned int> recv_code_;
    mfem::Vector recv_vec_;
};
} // namespace

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

    const int field_dim = pmesh.SpaceDimension();
    MFEM_VERIFY(pmesh.Dimension() == field_dim,
                "BuildEFieldH1Projection expects mesh.Dimension() == mesh.SpaceDimension().");
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
    // Use strict projection-specific tolerances; solver.* tolerances are tuned for
    // the main PDE solve and may be very loose defaults (e.g. atol=1, rtol=0).
    const int proj_maxiter = (cfg.solver.maxiter > 500) ? cfg.solver.maxiter : 500;
    cg.SetRelTol(1e-12);
    cg.SetAbsTol(0.0);
    cg.SetMaxIter(proj_maxiter);
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
        point_pos[p * space_dim + 0] = CenteredGridCoordinate(bbmin3[0], spacing[0], i);
        point_pos[p * space_dim + 1] = CenteredGridCoordinate(bbmin3[1], spacing[1], j);
        if (space_dim == 3)
        {
            const double z = CenteredGridCoordinate(bbmin3[2], spacing[2], k);
            point_pos[p * space_dim + 2] = z;
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

    const int space_dim = pmesh.SpaceDimension();
    const int mesh_dim = pmesh.Dimension();
    const int vdim = space_dim;
    MFEM_VERIFY(mesh_dim == space_dim,
                "SampleGradV expects mesh.Dimension() == mesh.SpaceDimension().");

    // 1) Redistribute point info to owning ranks (local elem ids + local ref coords + local codes)
    Array<unsigned int> recv_elem;
    Array<unsigned int> recv_code;
    Vector recv_ref_by_vdim; // byVDIM point-major: recv_ref[i*space_dim + d]
    finder.DistributePointInfoToOwningMPIRanks(recv_elem, recv_ref_by_vdim, recv_code);

    const int nrecv = recv_elem.Size();
    MFEM_VERIFY(recv_code.Size() == nrecv, "recv_code size mismatch");
    MFEM_VERIFY(recv_ref_by_vdim.Size() == nrecv * space_dim, "recv_ref size mismatch");

    // 2) Evaluate locally on owning ranks
    GradientGridFunctionCoefficient gradV(&V);
    Vector g(space_dim);
    IntegrationPoint ip;

    // Local interpolation buffer for the received points in byVDIM point-major layout.
    Vector vals_local;
    vals_local.SetSize(nrecv * vdim);
    vals_local = 0.0;

    for (int i = 0; i < nrecv; ++i)
    {
        const unsigned int c = recv_code[i];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf)) { continue; }

        const double r0 = recv_ref_by_vdim[i * space_dim + 0];
        const double r1 = recv_ref_by_vdim[i * space_dim + 1];
        const double r2 = (space_dim == 3) ? recv_ref_by_vdim[i * space_dim + 2] : 0.0;

        if (space_dim == 2) { ip.Set2(r0, r1); }
        else                { ip.Set3(r0, r1, r2); }

        const int e = static_cast<int>(recv_elem[i]);
        if (!IsLocalElementId(pmesh, e))
        {
            for (int d = 0; d < vdim; ++d)
            {
                vals_local[i * vdim + d] = std::numeric_limits<double>::quiet_NaN();
            }
            continue;
        }
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_VERIFY(T != nullptr, "Null ElementTransformation");

        T->SetIntPoint(&ip);
        gradV.Eval(g, *T, ip);
        g *= -1.0;

        // byVDIM point-major layout
        vals_local[i * vdim + 0] = g[0];
        vals_local[i * vdim + 1] = g[1];
        if (space_dim == 3) { vals_local[i * vdim + 2] = g[2]; }
    }

    // 3) Gather interpolated values back to the ranks that originated the query
    Vector vals_global;
    vals_global.SetSize(N * vdim);
    vals_global = 0.0;

    finder.DistributeInterpolatedValues(vals_local, vdim, Ordering::byVDIM, vals_global);

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
        const double ez = (space_dim == 3) ? vals_global[p * vdim + 2] : 0.0;
        const bool finite_eval = std::isfinite(ex) && std::isfinite(ey) &&
                                 (space_dim != 3 || std::isfinite(ez));

        if (!finite_eval)
        {
            out.valid[p] = 0;
            continue;
        }

        out.valid[p] = 1;
        out.Ex[p] = ex;
        out.Ey[p] = ey;
        if (space_dim == 3) { out.Ez[p] = ez; }
        out.Emag[p] = std::sqrt(ex*ex + ey*ey + ez*ez);
    }
}


static void SampleH1Projection(mfem::ParMesh &pmesh,
                               const mfem::ParGridFunction &E_h, // vdim == space_dim
                               mfem::FindPointsGSLIB &finder,
                               const mfem::Array<unsigned int> &code, // size N (original query ordering)
                               const int N,
                               const bool accept_surface_projection,
                               GridSample &out)
{
    using namespace mfem;

    const int space_dim = pmesh.SpaceDimension();
    const int mesh_dim = pmesh.Dimension();
    const int vdim = space_dim;
    MFEM_VERIFY(mesh_dim == space_dim,
                "SampleH1Projection expects mesh.Dimension() == mesh.SpaceDimension().");

    // 1) Redistribute point info to owning ranks
    Array<unsigned int> recv_elem;
    Array<unsigned int> recv_code;
    Vector recv_ref_by_vdim; // byVDIM point-major: ref[i*space_dim + d]
    finder.DistributePointInfoToOwningMPIRanks(recv_elem, recv_ref_by_vdim, recv_code);

    const int nrecv = recv_elem.Size();
    MFEM_VERIFY(recv_code.Size() == nrecv, "recv_code size mismatch");
    MFEM_VERIFY(recv_ref_by_vdim.Size() == nrecv * space_dim, "recv_ref size mismatch");

    // 2) Evaluate locally on owning ranks
    Vector g(space_dim);
    IntegrationPoint ip;

    // Local interpolation buffer for received points, byVDIM point-major.
    Vector vals_local(nrecv * vdim);
    vals_local = 0.0;

    for (int i = 0; i < nrecv; ++i)
    {
        const unsigned int c = recv_code[i];
        const bool inside  = (c == 0);
        const bool on_surf = (c == 1);
        if (!inside && !(accept_surface_projection && on_surf)) { continue; }

        const double r0 = recv_ref_by_vdim[i * space_dim + 0];
        const double r1 = recv_ref_by_vdim[i * space_dim + 1];
        const double r2 = (space_dim == 3) ? recv_ref_by_vdim[i * space_dim + 2] : 0.0;

        if (space_dim == 2) { ip.Set2(r0, r1); }
        else                { ip.Set3(r0, r1, r2); }

        const int e = static_cast<int>(recv_elem[i]);
        if (!IsLocalElementId(pmesh, e))
        {
            for (int d = 0; d < vdim; ++d)
            {
                vals_local[i * vdim + d] = std::numeric_limits<double>::quiet_NaN();
            }
            continue;
        }
        ElementTransformation *T = pmesh.GetElementTransformation(e);
        MFEM_VERIFY(T != nullptr, "Null ElementTransformation");

        T->SetIntPoint(&ip);

        // Evaluate the projected continuous vector field at this point
        E_h.GetVectorValue(*T, ip, g);

        // byVDIM point-major layout
        vals_local[i * vdim + 0] = g[0];
        vals_local[i * vdim + 1] = g[1];
        if (space_dim == 3) { vals_local[i * vdim + 2] = g[2]; }
    }

    // 3) Gather interpolated values back to original querying ranks/order
    Vector vals_global(N * vdim);
    vals_global = 0.0;

    finder.DistributeInterpolatedValues(vals_local, vdim, Ordering::byVDIM, vals_global);

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
        const double ez = (space_dim == 3) ? vals_global[p * vdim + 2] : 0.0;
        const bool finite_eval = std::isfinite(ex) && std::isfinite(ey) &&
                                 (space_dim != 3 || std::isfinite(ez));

        if (!finite_eval)
        {
            out.valid[p] = 0;
            continue;
        }

        out.valid[p] = 1;
        out.Ex[p] = ex;
        out.Ey[p] = ey;
        if (space_dim == 3) { out.Ez[p] = ez; }
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

        point_pos[lp * space_dim + 0] = CenteredGridCoordinate(bbmin3[0], spacing[0], i);
        point_pos[lp * space_dim + 1] = CenteredGridCoordinate(bbmin3[1], spacing[1], j);
        if (space_dim == 3)
        {
            point_pos[lp * space_dim + 2] = CenteredGridCoordinate(bbmin3[2], spacing[2], k);
        }
    }

    return point_pos;
}

// Legacy reintegrated E sampler (kept intentionally close to the original implementation).
[[maybe_unused]] static void SampleEFieldOnCartesianGridLegacy(mfem::ParMesh &pmesh,
                                              const mfem::ParGridFunction &V,
                                              const int Nx, const int Ny, const int Nz,
                                              GridSample &out,
                                              const Config &cfg,
                                              const bool H1_project,
                                              const bool accept_surface_projection)
{
    using namespace mfem;
    (void)H1_project; // legacy path: always direct gradient sampling
    (void)accept_surface_projection;

    MPI_Comm comm = pmesh.GetComm();
    int myid = 0, np = 1;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &np);

    const int space_dim = pmesh.SpaceDimension();
    const int field_dim = pmesh.Dimension();
    MFEM_VERIFY(space_dim == 2 || space_dim == 3, "Only 2D/3D physical space supported");

    // Legacy behavior: local bounding box (no global Allreduce).
    Vector bbmin_s(space_dim), bbmax_s(space_dim);
    pmesh.GetBoundingBox(bbmin_s, bbmax_s);

    Vector bbmin3(3), bbmax3(3);
    bbmin3 = 0.0; bbmax3 = 0.0;
    for (int d = 0; d < space_dim; ++d)
    {
        bbmin3[d] = bbmin_s[d];
        bbmax3[d] = bbmax_s[d];
    }
    ApplyInterpolationBoundsFromConfig(cfg, space_dim, bbmin3, bbmax3);

    Vector origin(3), spacing(3);
    origin = bbmin3;
    spacing = 0.0;

    MFEM_VERIFY(Nx > 0 && Ny > 0 && (space_dim != 3 || Nz > 0),
                "Nx, Ny (and Nz for 3D) must be > 0 for centered sampling.");
    spacing[0] = CenteredGridSpacing(bbmin3[0], bbmax3[0], Nx);
    spacing[1] = CenteredGridSpacing(bbmin3[1], bbmax3[1], Ny);
    if (space_dim == 3) { spacing[2] = CenteredGridSpacing(bbmin3[2], bbmax3[2], Nz); }
    else                { spacing[2] = 0.0; }
    origin[0] = CenteredGridCoordinate(bbmin3[0], spacing[0], 0);
    origin[1] = CenteredGridCoordinate(bbmin3[1], spacing[1], 0);
    if (space_dim == 3) { origin[2] = CenteredGridCoordinate(bbmin3[2], spacing[2], 0); }

    const int N = Nx * Ny * Nz;

    const int base = (np > 0) ? (N / np) : 0;
    const int rem  = (np > 0) ? (N % np) : 0;
    const int nloc = base + (myid < rem ? 1 : 0);
    const int p0   = myid * base + (myid < rem ? myid : rem);

    const bool need_dummy = (nloc == 0);
    const int nq = need_dummy ? 1 : nloc;

    // Legacy byNODES query layout: [x-block][y-block][z-block]
    Vector point_pos(space_dim * nq);
    if (!need_dummy)
    {
        for (int lp = 0; lp < nq; ++lp)
        {
            const int p = p0 + lp;
            const int i = p % Nx;
            const int t = p / Nx;
            const int j = t % Ny;
            const int k = t / Ny;

            point_pos[0 * nq + lp] = CenteredGridCoordinate(bbmin3[0], spacing[0], i);
            point_pos[1 * nq + lp] = CenteredGridCoordinate(bbmin3[1], spacing[1], j);
            if (space_dim == 3) { point_pos[2 * nq + lp] = CenteredGridCoordinate(bbmin3[2], spacing[2], k); }
        }
    }
    else
    {
        point_pos = 0.0;
        point_pos[0] = origin[0];
        if (space_dim >= 2) { point_pos[1] = origin[1]; }
        if (space_dim == 3) { point_pos[2] = origin[2]; }
    }

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

    pmesh.EnsureNodes();
    MFEM_VERIFY(pmesh.GetNodes() != nullptr, "Mesh nodes are required for FindPointsGSLIB");

    FindPointsGSLIB finder(comm);
    finder.Setup(pmesh);
    finder.FindPoints(point_pos, Ordering::byNODES);

    Array<unsigned int> recv_elem, recv_code;
    Vector recv_ref; // legacy path keeps original indexing assumptions
    finder.DistributePointInfoToOwningMPIRanks(recv_elem, recv_ref, recv_code);

    const int nrecv = recv_elem.Size();
    const int sdim = space_dim;
    Vector recv_E(nrecv * sdim);
    Vector grad(sdim);

    for (int j = 0; j < nrecv; ++j)
    {
        if (recv_code[j] == 2)
        {
            recv_E(j * sdim + 0) = std::numeric_limits<double>::quiet_NaN();
            recv_E(j * sdim + 1) = std::numeric_limits<double>::quiet_NaN();
            if (sdim == 3) { recv_E(j * sdim + 2) = std::numeric_limits<double>::quiet_NaN(); }
            continue;
        }

        const int elem_id = static_cast<int>(recv_elem[j]);
        if (!IsLocalElementId(pmesh, elem_id))
        {
            recv_E(j * sdim + 0) = std::numeric_limits<double>::quiet_NaN();
            recv_E(j * sdim + 1) = std::numeric_limits<double>::quiet_NaN();
            if (sdim == 3) { recv_E(j * sdim + 2) = std::numeric_limits<double>::quiet_NaN(); }
            continue;
        }
        ElementTransformation *T = pmesh.GetElementTransformation(elem_id);
        MFEM_VERIFY(T != nullptr, "Null ElementTransformation");

        IntegrationPoint ip;
        if (sdim == 2) { ip.Set2(recv_ref(j * sdim + 0), recv_ref(j * sdim + 1)); }
        else           { ip.Set3(recv_ref(j * sdim + 0), recv_ref(j * sdim + 1), recv_ref(j * sdim + 2)); }

        T->SetIntPoint(&ip);
        V.GetGradient(*T, grad);
        grad *= -1.0;

        recv_E(j * sdim + 0) = grad[0];
        recv_E(j * sdim + 1) = grad[1];
        if (sdim == 3) { recv_E(j * sdim + 2) = grad[2]; }
    }

    Vector Evals_flat(nq * sdim);
    finder.DistributeInterpolatedValues(recv_E, sdim, Ordering::byVDIM, Evals_flat);

    const Array<unsigned int> &code_local = finder.GetCode();
    for (int lp = 0; lp < nloc; ++lp)
    {
        const double ex = Evals_flat(lp * sdim + 0);
        const double ey = Evals_flat(lp * sdim + 1);
        const double ez = (sdim == 3) ? Evals_flat(lp * sdim + 2) : 0.0;

        local.Ex[lp] = ex;
        local.Ey[lp] = ey;
        local.Ez[lp] = ez;

        const double emag = std::sqrt(ex * ex + ey * ey + ez * ez);
        local.Emag[lp] = std::isfinite(emag) ? emag : std::numeric_limits<double>::quiet_NaN();

        // Legacy validity behavior.
        local.valid[lp] = (code_local[lp] == 0) ? 1 : 0;
    }

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

    out.has_E = true;
}
void SampleEFieldOnCartesianGrid(mfem::ParMesh &pmesh,
                                 const mfem::ParGridFunction &V,
                                 const int Nx, const int Ny, const int Nz,
                                 GridSample &out,
                                 const Config &cfg,
                                 const bool H1_project,
                                 const bool accept_surface_projection) // currently unused
{
    using namespace mfem;
    (void)accept_surface_projection;

    MPI_Comm comm = pmesh.GetComm();
    int myid = 0, np = 1;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &np);

    const int space_dim = pmesh.SpaceDimension();
    const int mesh_dim = pmesh.Dimension();
    const int vdim = space_dim;
    MFEM_VERIFY(mesh_dim == vdim,
                "SampleEFieldOnCartesianGrid expects mesh.Dimension() == mesh.SpaceDimension().");
    MFEM_VERIFY(space_dim == 2 || space_dim == 3, "Only 2D/3D physical space supported");

    // ------------------------------------------------------------------
    // Global bounding box
    // ------------------------------------------------------------------
    Vector bbmin3(3), bbmax3(3);
    ComputeGlobalBoundingBox(pmesh, bbmin3, bbmax3);
    ApplyInterpolationBoundsFromConfig(cfg, space_dim, bbmin3, bbmax3);

    Vector origin(3), spacing(3);
    origin = bbmin3;
    spacing = 0.0;

    MFEM_VERIFY(Nx > 0 && Ny > 0 && (space_dim != 3 || Nz > 0),
                "Nx, Ny (and Nz for 3D) must be > 0 for centered sampling.");
    spacing[0] = CenteredGridSpacing(bbmin3[0], bbmax3[0], Nx);
    spacing[1] = CenteredGridSpacing(bbmin3[1], bbmax3[1], Ny);
    if (space_dim == 3) { spacing[2] = CenteredGridSpacing(bbmin3[2], bbmax3[2], Nz); }
    else                { spacing[2] = 0.0; }
    origin[0] = CenteredGridCoordinate(bbmin3[0], spacing[0], 0);
    origin[1] = CenteredGridCoordinate(bbmin3[1], spacing[1], 0);
    if (space_dim == 3) { origin[2] = CenteredGridCoordinate(bbmin3[2], spacing[2], 0); }

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

    // Create vector of coordinates in byVDIM point-major layout: [x0,y0,(z0), x1,y1,(z1), ...].
    Vector point_pos(space_dim * nq);   // FIXED: no ordering enum in constructor

    if (!need_dummy)
    {
        // Build local subset deterministically in byVDIM point-major order.
        for (int lp = 0; lp < nq; ++lp)
        {
            const int p = p0 + lp;
            const int i = p % Nx;
            const int t = p / Nx;
            const int j = t % Ny;
            const int k = t / Ny;

            point_pos[lp * space_dim + 0] = CenteredGridCoordinate(bbmin3[0], spacing[0], i);
            point_pos[lp * space_dim + 1] = CenteredGridCoordinate(bbmin3[1], spacing[1], j);
            if (space_dim == 3)
                point_pos[lp * space_dim + 2] = CenteredGridCoordinate(bbmin3[2], spacing[2], k);
        }
    }
    else
    {
        // Dummy point inside bbox, still byVDIM point-major layout.
        point_pos = 0.0;
        point_pos[0] = origin[0];
        if (space_dim >= 2) point_pos[1] = origin[1];
        if (space_dim == 3) point_pos[2] = origin[2];
    }

    // ------------------------------------------------------------------
    // Local output containers (size nq, but only first nloc are meaningful)
    // ------------------------------------------------------------------
    GridSample local;
    local.dim = vdim;
    local.Nx = Nx; local.Ny = Ny; local.Nz = Nz;
    local.origin.SetSize(3); local.spacing.SetSize(3);
    local.origin = origin; local.spacing = spacing;

    local.Ex.assign(nq, std::numeric_limits<double>::quiet_NaN());
    local.Ey.assign(nq, std::numeric_limits<double>::quiet_NaN());
    local.Ez.assign(nq, 0.0);
    local.Emag.assign(nq, std::numeric_limits<double>::quiet_NaN());
    local.valid.assign(nq, static_cast<unsigned char>(0));

    // ------------------------------------------------------------------
    // Point location/evaluation (collective) using MPI tracer field sampler
    // ------------------------------------------------------------------
    pmesh.EnsureNodes();
    MFEM_VERIFY(pmesh.GetNodes() != nullptr, "Mesh nodes are required for FindPointsGSLIB");

    FindPointsGSLIB finder(comm);
    finder.Setup(pmesh);
    finder.SetDefaultInterpolationValue(std::numeric_limits<double>::quiet_NaN());
    finder.SetL2AvgType(mfem::FindPointsGSLIB::AvgType::ARITHMETIC);
    const int sdim = space_dim;
    Vector Evals_flat(nq * sdim);
    Array<unsigned int> code_initial;
    Array<unsigned int> elem_initial;
    Vector ref_initial;

    if (H1_project)
    {
        ProjectedH1VectorField projected = BuildEFieldH1Projection(pmesh, V, cfg);
        MPITracerStyleVectorFieldEvaluator evaluator(pmesh, *(projected.E), finder, sdim);
        evaluator.EvaluateVectorField(point_pos, nq, Evals_flat);

        code_initial = finder.GetCode();
        elem_initial = finder.GetElem();
        ref_initial  = finder.GetReferencePosition();
        MFEM_VERIFY(code_initial.Size() == nq, "Unexpected code array size.");
        MFEM_VERIFY(elem_initial.Size() == nq, "Unexpected element array size.");
        MFEM_VERIFY(ref_initial.Size() == nq * sdim, "Unexpected reference array size.");

        ResolveCode1Samples(pmesh, spacing, cfg.interp.code1_mode,
                            code_initial, elem_initial, ref_initial,
                            nloc, sdim, sdim,
                            [&](const Vector &pts, const int npts, Vector &vals)
                            {
                                evaluator.EvaluateVectorField(pts, npts, vals);
                            },
                            Evals_flat);
    }
    else
    {
        MPITracerStyleFieldEvaluator evaluator(pmesh, V, finder, sdim);
        evaluator.EvaluateE(point_pos, nq, Evals_flat);

        code_initial = finder.GetCode();
        elem_initial = finder.GetElem();
        ref_initial  = finder.GetReferencePosition();
        MFEM_VERIFY(code_initial.Size() == nq, "Unexpected code array size.");
        MFEM_VERIFY(elem_initial.Size() == nq, "Unexpected element array size.");
        MFEM_VERIFY(ref_initial.Size() == nq * sdim, "Unexpected reference array size.");

        ResolveCode1Samples(pmesh, spacing, cfg.interp.code1_mode,
                            code_initial, elem_initial, ref_initial,
                            nloc, sdim, sdim,
                            [&](const Vector &pts, const int npts, Vector &vals)
                            {
                                evaluator.EvaluateE(pts, npts, vals);
                            },
                            Evals_flat);
    }

    // ------------------------------------------------------------------
    // Fill local arrays for the first nloc real points
    // ------------------------------------------------------------------
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

        // Accept found points only if evaluation produced finite values.
        const bool found = (code_initial[lp] != 2);
        const bool finite_eval = std::isfinite(ex) && std::isfinite(ey) &&
                                 (sdim != 3 || std::isfinite(ez));
        local.valid[lp] = (found && finite_eval) ? 1 : 0;
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
        out.dim = vdim;
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
    if (vdim == 3)
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

    out.has_E = true;
}

void SampleVFieldOnCartesianGrid(mfem::ParMesh &pmesh,
                                 const mfem::ParGridFunction &V,
                                 const int Nx, const int Ny, const int Nz,
                                 GridSample &out,
                                 const Config &cfg)
{
    using namespace mfem;

    MPI_Comm comm = pmesh.GetComm();
    int myid = 0, np = 1;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &np);

    const int space_dim = pmesh.SpaceDimension();
    const int mesh_dim = pmesh.Dimension();
    const int vdim = space_dim;
    MFEM_VERIFY(mesh_dim == vdim,
                "SampleVFieldOnCartesianGrid expects mesh.Dimension() == mesh.SpaceDimension().");
    MFEM_VERIFY(space_dim == 2 || space_dim == 3, "Only 2D/3D physical space supported");

    Vector bbmin3(3), bbmax3(3);
    ComputeGlobalBoundingBox(pmesh, bbmin3, bbmax3);
    ApplyInterpolationBoundsFromConfig(cfg, space_dim, bbmin3, bbmax3);

    Vector origin(3), spacing(3);
    origin = bbmin3;
    spacing = 0.0;

    MFEM_VERIFY(Nx > 0 && Ny > 0 && (space_dim != 3 || Nz > 0),
                "Nx, Ny (and Nz for 3D) must be > 0 for centered sampling.");
    spacing[0] = CenteredGridSpacing(bbmin3[0], bbmax3[0], Nx);
    spacing[1] = CenteredGridSpacing(bbmin3[1], bbmax3[1], Ny);
    if (space_dim == 3) { spacing[2] = CenteredGridSpacing(bbmin3[2], bbmax3[2], Nz); }
    else                { spacing[2] = 0.0; }
    origin[0] = CenteredGridCoordinate(bbmin3[0], spacing[0], 0);
    origin[1] = CenteredGridCoordinate(bbmin3[1], spacing[1], 0);
    if (space_dim == 3) { origin[2] = CenteredGridCoordinate(bbmin3[2], spacing[2], 0); }

    const int N = Nx * Ny * Nz;

    const int base = (np > 0) ? (N / np) : 0;
    const int rem  = (np > 0) ? (N % np) : 0;
    const int nloc = base + (myid < rem ? 1 : 0);
    const int p0   = myid * base + (myid < rem ? myid : rem);

    const bool need_dummy = (nloc == 0);
    const int nq = need_dummy ? 1 : nloc;

    Vector point_pos(space_dim * nq);
    if (!need_dummy)
    {
        for (int lp = 0; lp < nq; ++lp)
        {
            const int p = p0 + lp;
            const int i = p % Nx;
            const int t = p / Nx;
            const int j = t % Ny;
            const int k = t / Ny;

            point_pos[lp * space_dim + 0] = CenteredGridCoordinate(bbmin3[0], spacing[0], i);
            point_pos[lp * space_dim + 1] = CenteredGridCoordinate(bbmin3[1], spacing[1], j);
            if (space_dim == 3) { point_pos[lp * space_dim + 2] = CenteredGridCoordinate(bbmin3[2], spacing[2], k); }
        }
    }
    else
    {
        point_pos = 0.0;
        point_pos[0] = origin[0];
        if (space_dim >= 2) { point_pos[1] = origin[1]; }
        if (space_dim == 3) { point_pos[2] = origin[2]; }
    }

    std::vector<double> local_V(static_cast<std::size_t>(nq), std::numeric_limits<double>::quiet_NaN());
    std::vector<unsigned char> local_valid(static_cast<std::size_t>(nq), static_cast<unsigned char>(0));

    pmesh.EnsureNodes();
    MFEM_VERIFY(pmesh.GetNodes() != nullptr, "Mesh nodes are required for FindPointsGSLIB");

    FindPointsGSLIB finder(comm);
    finder.Setup(pmesh);
    finder.SetDefaultInterpolationValue(std::numeric_limits<double>::quiet_NaN());
    finder.SetL2AvgType(mfem::FindPointsGSLIB::AvgType::ARITHMETIC);

    const int sdim = space_dim;
    MPITracerStylePotentialEvaluator evaluator(pmesh, V, finder, sdim);
    Vector Vvals_flat(nq);
    evaluator.EvaluateV(point_pos, nq, Vvals_flat);

    Array<unsigned int> code_initial = finder.GetCode();
    Array<unsigned int> elem_initial = finder.GetElem();
    Vector ref_initial = finder.GetReferencePosition();
    MFEM_VERIFY(code_initial.Size() == nq, "Unexpected code array size.");
    MFEM_VERIFY(elem_initial.Size() == nq, "Unexpected element array size.");
    MFEM_VERIFY(ref_initial.Size() == nq * sdim, "Unexpected reference array size.");

    ResolveCode1Samples(pmesh, spacing, cfg.interp.code1_mode,
                        code_initial, elem_initial, ref_initial,
                        nloc, sdim, 1,
                        [&](const Vector &pts, const int npts, Vector &vals)
                        {
                            evaluator.EvaluateV(pts, npts, vals);
                        },
                        Vvals_flat);

    for (int lp = 0; lp < nloc; ++lp)
    {
        local_V[static_cast<std::size_t>(lp)] = Vvals_flat[lp];
        const bool finite_v = std::isfinite(Vvals_flat[lp]);
        local_valid[static_cast<std::size_t>(lp)] = ((code_initial[lp] != 2) && finite_v) ? 1 : 0;
    }

    const int send_n = nloc;
    std::vector<int> recvcounts(np, 0), displs(np, 0);
    for (int r = 0; r < np; ++r)
    {
        const int rn  = base + (r < rem ? 1 : 0);
        const int rp0 = r * base + (r < rem ? r : rem);
        recvcounts[r] = rn;
        displs[r]     = rp0;
    }

    std::vector<unsigned char> gathered_valid;
    const bool had_valid_in = (myid == 0) && (out.valid.size() == static_cast<std::size_t>(N));
    if (myid == 0)
    {
        out.dim = vdim;
        out.Nx = Nx; out.Ny = Ny; out.Nz = Nz;
        out.origin.SetSize(3); out.spacing.SetSize(3);
        out.origin = origin; out.spacing = spacing;
        out.V.assign(static_cast<std::size_t>(N), std::numeric_limits<double>::quiet_NaN());
        gathered_valid.assign(static_cast<std::size_t>(N), static_cast<unsigned char>(0));
    }

    double dummy_d = 0.0;
    unsigned char dummy_uc = 0;

    const double *sendV = (send_n > 0) ? local_V.data() : nullptr;
    const unsigned char *sendVal = (send_n > 0) ? local_valid.data() : nullptr;

    double *recvV = (myid == 0) ? out.V.data() : &dummy_d;
    unsigned char *recvVal = (myid == 0) ? gathered_valid.data() : &dummy_uc;

    MPI_Gatherv(sendV, send_n, MPI_DOUBLE,
                recvV, recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, comm);
    MPI_Gatherv(sendVal, send_n, MPI_UNSIGNED_CHAR,
                recvVal, recvcounts.data(), displs.data(), MPI_UNSIGNED_CHAR,
                0, comm);

    if (myid == 0)
    {
        if (had_valid_in)
        {
            for (std::size_t p = 0; p < out.valid.size(); ++p)
            {
                out.valid[p] = static_cast<uint8_t>(out.valid[p] && gathered_valid[p]);
            }
        }
        else
        {
            out.valid.assign(gathered_valid.begin(), gathered_valid.end());
        }
    }

    out.has_V = true;
}




static inline void h5_check(const herr_t st, const char *what)
{
  if (st < 0) { throw std::runtime_error(std::string("HDF5 error: ") + what); }
}

static void WriteGridSampleBinary(const GridSample &g,
                                  const std::filesystem::path &out_dir,
                                  const std::string &filename = "interpolated.h5")
{
  if (g.Nx <= 0 || g.Ny <= 0 || g.Nz <= 0)
    throw std::runtime_error("WriteGridSampleBinary: Nx, Ny, Nz must all be > 0.");

  const hsize_t Nx  = static_cast<hsize_t>(g.Nx);
  const hsize_t Ny  = static_cast<hsize_t>(g.Ny);
  const hsize_t Nz  = static_cast<hsize_t>(g.Nz);
  const hsize_t dim = static_cast<hsize_t>(g.dim);
  const bool write_E = g.has_E;
  const bool write_V = g.has_V;

  const size_t N = static_cast<size_t>(g.Nx) * static_cast<size_t>(g.Ny) * static_cast<size_t>(g.Nz);

  if (g.dim != 2 && g.dim != 3)
    throw std::runtime_error("WriteGridSampleBinary: g.dim must be 2 or 3.");

  if (!write_E && !write_V)
    throw std::runtime_error("WriteGridSampleBinary: nothing to write (both has_E and has_V are false).");

  if (write_E && (g.Ex.size() != N || g.Ey.size() != N || (g.dim == 3 && g.Ez.size() != N)))
    throw std::runtime_error("WriteGridSampleBinary: Ex/Ey/Ez sizes do not match grid size.");

  if (write_V && g.V.size() != N)
    throw std::runtime_error("WriteGridSampleBinary: V size does not match grid size.");

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

  std::vector<double> E;
  if (write_E)
  {
    E.assign(N * static_cast<size_t>(g.dim), std::numeric_limits<double>::quiet_NaN());
    for (size_t p = 0; p < N; ++p)
    {
      E[p * static_cast<size_t>(g.dim) + 0] = g.Ex[p];
      E[p * static_cast<size_t>(g.dim) + 1] = g.Ey[p];
      if (g.dim == 3) { E[p * static_cast<size_t>(g.dim) + 2] = g.Ez[p]; }
    }
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
    if (write_E)
    {
      // ----------------------------------------------------------------
      // /E/field
      // ----------------------------------------------------------------
      hid_t grpE = H5Gcreate2(file, "/E", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (grpE < 0) throw std::runtime_error("HDF5 error: H5Gcreate2(/E)");

      hsize_t d4[4] = {Nz, Ny, Nx, dim};
      hid_t spaceE = H5Screate_simple(4, d4, nullptr);
      if (spaceE < 0) { H5Gclose(grpE); throw std::runtime_error("HDF5 error: H5Screate_simple(/E/field)"); }

      hid_t dcplE = H5Pcreate(H5P_DATASET_CREATE);
      if (dcplE < 0) { H5Sclose(spaceE); H5Gclose(grpE); throw std::runtime_error("HDF5 error: H5Pcreate(DCPL E)"); }

      hsize_t chunk4[4] = {
        std::min<hsize_t>(Nz, 16),
        std::min<hsize_t>(Ny, 64),
        std::min<hsize_t>(Nx, 64),
        dim
      };
      h5_check(H5Pset_chunk(dcplE, 4, chunk4), "H5Pset_chunk(E)");

      hid_t dsetE = H5Dcreate2(grpE, "field", H5T_IEEE_F64LE, spaceE,
                               H5P_DEFAULT, dcplE, H5P_DEFAULT);
      if (dsetE < 0)
      {
        H5Pclose(dcplE);
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
      write_attr_f64_vec(dsetE, "origin",  origin);
      write_attr_f64_vec(dsetE, "spacing", spacing);

      {
        int axis_order[4] = {2, 1, 0, 3};
        hsize_t n[1] = {4};
        hid_t as = H5Screate_simple(1, n, nullptr);
        if (as < 0)
        {
          H5Dclose(dsetE);
          H5Pclose(dcplE);
          H5Sclose(spaceE);
          H5Gclose(grpE);
          throw std::runtime_error("HDF5 error: H5Screate_simple(axis_order E)");
        }
        hid_t at = H5Acreate2(dsetE, "axis_order", H5T_STD_I32LE, as, H5P_DEFAULT, H5P_DEFAULT);
        if (at < 0)
        {
          H5Sclose(as);
          H5Dclose(dsetE);
          H5Pclose(dcplE);
          H5Sclose(spaceE);
          H5Gclose(grpE);
          throw std::runtime_error("HDF5 error: H5Acreate2(axis_order E)");
        }
        h5_check(H5Awrite(at, H5T_NATIVE_INT, axis_order), "H5Awrite(axis_order E)");
        H5Aclose(at); H5Sclose(as);
      }

      H5Dclose(dsetE);
      H5Pclose(dcplE);
      H5Sclose(spaceE);
      H5Gclose(grpE);
    }

    if (write_V)
    {
      // ----------------------------------------------------------------
      // /V/field
      // ----------------------------------------------------------------
      hid_t grpV = H5Gcreate2(file, "/V", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (grpV < 0) throw std::runtime_error("HDF5 error: H5Gcreate2(/V)");

      hsize_t d3[3] = {Nz, Ny, Nx};
      hid_t spaceV = H5Screate_simple(3, d3, nullptr);
      if (spaceV < 0) { H5Gclose(grpV); throw std::runtime_error("HDF5 error: H5Screate_simple(/V/field)"); }

      hid_t dcplV = H5Pcreate(H5P_DATASET_CREATE);
      if (dcplV < 0) { H5Sclose(spaceV); H5Gclose(grpV); throw std::runtime_error("HDF5 error: H5Pcreate(DCPL V)"); }

      hsize_t chunk3[3] = {
        std::min<hsize_t>(Nz, 16),
        std::min<hsize_t>(Ny, 64),
        std::min<hsize_t>(Nx, 64)
      };
      h5_check(H5Pset_chunk(dcplV, 3, chunk3), "H5Pset_chunk(V)");

      hid_t dsetV = H5Dcreate2(grpV, "field", H5T_IEEE_F64LE, spaceV,
                               H5P_DEFAULT, dcplV, H5P_DEFAULT);
      if (dsetV < 0)
      {
        H5Pclose(dcplV);
        H5Sclose(spaceV);
        H5Gclose(grpV);
        throw std::runtime_error("HDF5 error: H5Dcreate2(/V/field)");
      }

      h5_check(H5Dwrite(dsetV, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, g.V.data()),
               "H5Dwrite(/V/field)");

      write_attr_i32(dsetV, "Nx",  g.Nx);
      write_attr_i32(dsetV, "Ny",  g.Ny);
      write_attr_i32(dsetV, "Nz",  g.Nz);
      write_attr_f64_vec(dsetV, "origin",  origin);
      write_attr_f64_vec(dsetV, "spacing", spacing);

      {
        int axis_order[3] = {2, 1, 0};
        hsize_t n[1] = {3};
        hid_t as = H5Screate_simple(1, n, nullptr);
        if (as < 0)
        {
          H5Dclose(dsetV);
          H5Pclose(dcplV);
          H5Sclose(spaceV);
          H5Gclose(grpV);
          throw std::runtime_error("HDF5 error: H5Screate_simple(axis_order V)");
        }
        hid_t at = H5Acreate2(dsetV, "axis_order", H5T_STD_I32LE, as, H5P_DEFAULT, H5P_DEFAULT);
        if (at < 0)
        {
          H5Sclose(as);
          H5Dclose(dsetV);
          H5Pclose(dcplV);
          H5Sclose(spaceV);
          H5Gclose(grpV);
          throw std::runtime_error("HDF5 error: H5Acreate2(axis_order V)");
        }
        h5_check(H5Awrite(at, H5T_NATIVE_INT, axis_order), "H5Awrite(axis_order V)");
        H5Aclose(at); H5Sclose(as);
      }

      H5Dclose(dsetV);
      H5Pclose(dcplV);
      H5Sclose(spaceV);
      H5Gclose(grpV);
    }

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
      if (as < 0)
      {
        H5Gclose(grpG);
        throw std::runtime_error("HDF5 error: H5Screate_simple(coord_kind)");
      }
      hid_t at = H5Acreate2(grpG, "coord_kind", H5T_STD_I32LE, as, H5P_DEFAULT, H5P_DEFAULT);
      if (at < 0)
      {
        H5Sclose(as);
        H5Gclose(grpG);
        throw std::runtime_error("HDF5 error: H5Acreate2(coord_kind)");
      }
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
  const bool want_E_interp = cfg.interp.E;
  const bool want_V_interp = cfg.interp.V;

  for (const auto& run_dir : targets)
  {
    // Load data (must be called on all ranks if it builds ParMesh/ParGridFunction on world comm)
    SimulationResult result = load_results(cfg, run_dir);

    if (!result.success || !result.mesh || !result.V)
    {
      if (world_rank == 0)
      {
        std::cerr << "[INTERP] skipping " << run_dir
                  << " (load failed): " << result.error_message << "\n";
      }
      continue;
    }

    const int Nx = cfg.interp.Nx;
    const int Ny = cfg.interp.Ny;
    const int Nz = cfg.interp.Nz;
    const bool accept_surface = false;

    if (!want_E_interp && !want_V_interp)
    {
      if (world_rank == 0)
      {
        std::cout << "[INTERP] skipping " << run_dir
                  << " (both interpolate.E and interpolate.V are false)\n";
      }
      continue;
    }

    const auto out_dir = run_dir / "interpolated";
    if (world_rank == 0) { std::filesystem::create_directories(out_dir); }

    GridSample grid;
    grid.dim = result.mesh->SpaceDimension();
    const std::int64_t npts = static_cast<std::int64_t>(Nx) *
                              static_cast<std::int64_t>(Ny) *
                              static_cast<std::int64_t>(Nz);

    if (world_rank == 0)
    {
      std::cout << "[INTERP] run=" << run_dir
                << " points=" << npts
                << " Nx=" << Nx << " Ny=" << Ny << " Nz=" << Nz
                << " E=" << (want_E_interp ? "true" : "false")
                << " V=" << (want_V_interp ? "true" : "false")
                << " code1_mode=" << static_cast<int>(cfg.interp.code1_mode)
                << "\n";
    }

    if (want_E_interp)
    {
      SampleEFieldOnCartesianGrid(*result.mesh, *result.V,
                                  Nx, Ny, Nz, grid,
                                  cfg, cfg.interp.H1_project, accept_surface);
    }

    if (want_V_interp)
    {
      SampleVFieldOnCartesianGrid(*result.mesh, *result.V, Nx, Ny, Nz, grid, cfg);
    }

    if (world_rank == 0)
    {
      WriteGridSampleBinary(grid, out_dir, "interpolated.h5");
      std::cout << "[INTERP] wrote interpolated fields to: " << out_dir
                << " (E=" << (grid.has_E ? "true" : "false")
                << ", V=" << (grid.has_V ? "true" : "false") << ")\n";
    }
  }
  return 0;
}
