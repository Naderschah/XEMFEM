/*
Here we compute the field spread:

F_spread = (E_95 -E_5)/<E>


So we have V which is continuous on H1

But E is piecewise smooth and discontinuous across element faces. 

Discontinuities live on element interfaces, correspond to measure zero sets.
Integrals of the form int g(E) dx are well defined as long as g(E) is in L1 
which is the case for |E| and |E|^2. So integrating over all elements 
should be fine for abs E? 

And for percentiles we need to do sorting + Volume contribution weights



We compute the |E| field pointwise here as this will be removed from the results object.
*/
#include <mfem.hpp>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

#include "field_spread.h"

static bool InROI(const mfem::Vector &x, const Config &cfg)
{
    // x: physical coordinates (dim = 2 or 3)
    double r = x[0];
    double z = x[1];

    return (r >= cfg.optimize.r_min &&
            r <= cfg.optimize.r_max &&
            z >= cfg.optimize.z_min &&
            z <= cfg.optimize.z_max);
}

// Compute spread metric = (Q95 - Q5) / mean(|E|) in ROI
double computeFieldSpreadMetric(const Config &cfg,
                                const SimulationResult &result)
{
    using namespace mfem;
    static mfem::IntegrationRules g_IntRules(0, mfem::Quadrature1D::GaussLegendre);

    GridFunction &V   = *result.V;
    FiniteElementSpace &fes = *V.FESpace();
    Mesh &mesh        = *fes.GetMesh();
    const int dim     = mesh.Dimension();

    // Gradient of potential
    GradientGridFunctionCoefficient gradV(&V);

    std::vector<double> e_vals;
    std::vector<double> w_vals;

    e_vals.reserve(mesh.GetNE() * 4);
    w_vals.reserve(mesh.GetNE() * 4);

    double total_w = 0.0;
    double sum_w_e = 0.0;

    const int order    = fes.GetOrder(0);
    const int ir_order = 2 * order; // conservative

    // Element loop
    for (int el = 0; el < mesh.GetNE(); el++)
    {
        ElementTransformation *T = mesh.GetElementTransformation(el);
        const FiniteElement *fe  = fes.GetFE(el);
        const Geometry::Type geom = fe->GetGeomType();

        const IntegrationRule &ir = g_IntRules.Get((int)geom, ir_order);

        for (int i = 0; i < ir.GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir.IntPoint(i);
            T->SetIntPoint(&ip);

            // Physical coordinates
            Vector x(dim);
            T->Transform(ip, x);

            if (!InROI(x, cfg)) { continue; }

            // E = -∇V
            Vector grad(dim);
            gradV.Eval(grad, *T, ip);
            grad *= -1.0;

            const double e = grad.Norml2();              // |E|
            const double w = ip.weight * T->Weight();    // J * weight

            e_vals.push_back(e);
            w_vals.push_back(w);

            total_w += w;
            sum_w_e += w * e;
        }
    }

    if (total_w == 0.0 || e_vals.empty())
    {
        // ROI has no quadrature points / elements
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double mean_E = sum_w_e / total_w;

    // ------- Weighted 5th and 95th percentiles --------

    const std::size_t N = e_vals.size();

    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b)
              {
                  return e_vals[a] < e_vals[b];
              });

    const double target5  = cfg.field_spread.lower * total_w;
    const double target95 = cfg.field_spread.upper * total_w;

    double cum_w = 0.0;
    double Q5    = e_vals[idx.front()];
    double Q95   = e_vals[idx.back()];
    bool set5 = false, set95 = false;

    for (std::size_t k = 0; k < N; k++)
    {
        cum_w += w_vals[idx[k]];

        if (!set5 && cum_w >= target5)
        {
            Q5 = e_vals[idx[k]];
            set5 = true;
        }
        if (!set95 && cum_w >= target95)
        {
            Q95 = e_vals[idx[k]];
            set95 = true;
            break;
        }
    }

    const double spread = (Q95 - Q5) / mean_E;
    
    // Check its finite 
    if (!std::isfinite(spread))
    {
        double minE = std::numeric_limits<double>::infinity();
        double maxE = -std::numeric_limits<double>::infinity();

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
            << "  total_w  = " << total_w << "\n"
            << "  samples  = " << e_vals.size() << "\n"
            << "  min|E|   = " << minE << "\n"
            << "  max|E|   = " << maxE << "\n"
            << std::endl;
    }

    if (cfg.debug.debug) {std::cout << "Field Spread " << spread << std::endl;}
    return spread;
}