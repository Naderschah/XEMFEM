#include <mfem.hpp>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

#include "field_spread.h"

static bool InROI(const mfem::Vector &x, const Config &cfg)
{
    // Requires dim >= 2
    const double r = x[0];
    const double z = x[1];

    return (r >= cfg.optimize.r_min &&
            r <= cfg.optimize.r_max &&
            z >= cfg.optimize.z_min &&
            z <= cfg.optimize.z_max);
}

static void WeightedQuantilesOnRoot_Gather(
    MPI_Comm comm,
    const std::vector<double> &local_e,
    const std::vector<double> &local_w,
    double global_total_w,
    double lower_frac,
    double upper_frac,
    double &Qlo,
    double &Qhi)
{
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int local_n = static_cast<int>(local_e.size());
    MFEM_VERIFY(local_w.size() == local_e.size(), "local_e/local_w size mismatch");

    std::vector<int> counts(size, 0), displs(size, 0);
    MPI_Gather(&local_n, 1, MPI_INT,
               counts.data(), 1, MPI_INT, 0, comm);

    int global_n = 0;
    if (rank == 0)
    {
        for (int p = 0; p < size; ++p) { displs[p] = global_n; global_n += counts[p]; }
    }

    std::vector<double> all_e, all_w;
    if (rank == 0)
    {
        all_e.resize(global_n);
        all_w.resize(global_n);
    }

    MPI_Gatherv(local_e.data(), local_n, MPI_DOUBLE,
                rank == 0 ? all_e.data() : nullptr,
                counts.data(), displs.data(), MPI_DOUBLE, 0, comm);

    MPI_Gatherv(local_w.data(), local_n, MPI_DOUBLE,
                rank == 0 ? all_w.data() : nullptr,
                counts.data(), displs.data(), MPI_DOUBLE, 0, comm);

    if (rank == 0)
    {
        if (global_total_w <= 0.0 || global_n == 0)
        {
            Qlo = std::numeric_limits<double>::quiet_NaN();
            Qhi = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            std::vector<int> idx(global_n);
            std::iota(idx.begin(), idx.end(), 0);
            std::sort(idx.begin(), idx.end(),
                      [&](int a, int b) { return all_e[a] < all_e[b]; });

            const double target_lo = lower_frac * global_total_w;
            const double target_hi = upper_frac * global_total_w;

            double cum_w = 0.0;
            Qlo = all_e[idx.front()];
            Qhi = all_e[idx.back()];
            bool set_lo = false, set_hi = false;

            for (int k = 0; k < global_n; ++k)
            {
                cum_w += all_w[idx[k]];
                if (!set_lo && cum_w >= target_lo) { Qlo = all_e[idx[k]]; set_lo = true; }
                if (!set_hi && cum_w >= target_hi) { Qhi = all_e[idx[k]]; set_hi = true; break; }
            }
        }
    }

    MPI_Bcast(&Qlo, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&Qhi, 1, MPI_DOUBLE, 0, comm);
}

double computeFieldSpreadMetric(const Config &cfg,
                                const SimulationResult &result)
{
    using namespace mfem;

    static IntegrationRules g_IntRules(0, Quadrature1D::GaussLegendre);

    GridFunction &V         = *result.V;
    FiniteElementSpace &fes = *V.FESpace();
    Mesh &mesh              = *fes.GetMesh();
    const int dim           = mesh.Dimension();

    MFEM_VERIFY(dim >= 2, "InROI expects at least 2D coordinates (r,z)");

    // Detect parallel communicator (works in serial too)
    MPI_Comm comm = MPI_COMM_SELF;
    int rank = 0, size = 1;
    if (auto *pfes = dynamic_cast<ParFiniteElementSpace*>(&fes))
    {
        comm = pfes->GetComm();
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
    }

    GradientGridFunctionCoefficient gradV(&V);

    std::vector<double> e_vals;
    std::vector<double> w_vals;
    e_vals.reserve(mesh.GetNE() * 4);
    w_vals.reserve(mesh.GetNE() * 4);

    double local_total_w = 0.0;
    double local_sum_w_e = 0.0;

    const int order    = fes.GetOrder(0);
    const int ir_order = 2 * order;

    for (int el = 0; el < mesh.GetNE(); el++)
    {
        ElementTransformation *T = mesh.GetElementTransformation(el);
        const FiniteElement *fe  = fes.GetFE(el);
        const Geometry::Type geom = fe->GetGeomType();

        const IntegrationRule &ir = g_IntRules.Get((int)geom, ir_order);

        Vector x(dim), grad(dim);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            T->Transform(ip, x);
            if (!InROI(x, cfg)) { continue; }

            gradV.Eval(grad, *T, ip);
            grad *= -1.0;

            const double e = grad.Norml2();
            const double w = ip.weight * T->Weight();

            e_vals.push_back(e);
            w_vals.push_back(w);

            local_total_w += w;
            local_sum_w_e += w * e;
        }
    }

    // Global reductions for mean(|E|)
    double global_total_w = local_total_w;
    double global_sum_w_e = local_sum_w_e;

    if (size > 1)
    {
        MPI_Allreduce(&local_total_w, &global_total_w, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&local_sum_w_e, &global_sum_w_e, 1, MPI_DOUBLE, MPI_SUM, comm);
    }

    if (global_total_w == 0.0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double mean_E = global_sum_w_e / global_total_w;

    // Global weighted quantiles (gather-to-root for simplicity)
    double Q5  = std::numeric_limits<double>::quiet_NaN();
    double Q95 = std::numeric_limits<double>::quiet_NaN();

    if (size == 1)
    {
        // Serial: compute percentiles locally
        const std::size_t N = e_vals.size();
        if (N == 0) { return std::numeric_limits<double>::quiet_NaN(); }

        std::vector<std::size_t> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
                  [&](std::size_t a, std::size_t b) { return e_vals[a] < e_vals[b]; });

        const double target5  = cfg.field_spread.lower * global_total_w;
        const double target95 = cfg.field_spread.upper * global_total_w;

        double cum_w = 0.0;
        Q5 = e_vals[idx.front()];
        Q95 = e_vals[idx.back()];
        bool set5 = false, set95 = false;

        for (std::size_t k = 0; k < N; k++)
        {
            cum_w += w_vals[idx[k]];
            if (!set5 && cum_w >= target5)  { Q5 = e_vals[idx[k]]; set5 = true; }
            if (!set95 && cum_w >= target95){ Q95 = e_vals[idx[k]]; set95 = true; break; }
        }
    }
    else
    {
        WeightedQuantilesOnRoot_Gather(
            comm, e_vals, w_vals, global_total_w,
            cfg.field_spread.lower, cfg.field_spread.upper,
            Q5, Q95);
    }

    const double spread = (Q95 - Q5) / mean_E;

    if (!std::isfinite(spread) && rank == 0)
    {
        double minE = std::numeric_limits<double>::infinity();
        double maxE = -std::numeric_limits<double>::infinity();

        // local min/max only (optional: Allreduce for global min/max)
        for (double v : e_vals)
        {
            if (std::isfinite(v))
            {
                minE = std::min(minE, v);
                maxE = std::max(maxE, v);
            }
        }

        std::cout
            << "[FieldSpread NaN/Inf]\n"
            << "  spread   = " << spread << "\n"
            << "  Q5       = " << Q5 << "\n"
            << "  Q95      = " << Q95 << "\n"
            << "  mean_E   = " << mean_E << "\n"
            << "  total_w  = " << global_total_w << "\n"
            << "  local_samples(rank0) = " << e_vals.size() << "\n"
            << "  local_min|E|(rank0)  = " << minE << "\n"
            << "  local_max|E|(rank0)  = " << maxE << "\n"
            << std::endl;
    }

    if (cfg.debug.debug && rank == 0)
    {
        std::cout << "Field Spread " << spread << std::endl;
    }

    return spread;
}