#include <mfem.hpp>
#include <vector>
#include <cmath>
#include <omp.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>

#include "trace_fieldlines.h"

using namespace mfem;

static void clamp_axis(Vector& x, bool axisymmetric, double geom_tol)
{
    if (axisymmetric && x[0] <= geom_tol) {
        x[0] = geom_tol;
    }
}


// Build adjacency and characteristic size h per element
ElementAdjacency BuildAdjacency(mfem::ParMesh &mesh)
{
    ElementAdjacency adj;

    const int ne  = mesh.GetNE();
    const int dim = mesh.SpaceDimension();

    adj.neighbors.resize(ne);
    adj.h.resize(ne);

    mfem::Array<int> vert_ids;

    const mfem::Table &el_to_el = mesh.ElementToElementTable();

    for (int e = 0; e < ne; ++e)
    {
        // neighbors
        {
            mfem::Array<int> nbrs;
            el_to_el.GetRow(e, nbrs);
            adj.neighbors[e].assign(nbrs.begin(), nbrs.end());
        }

        // element "size": shortest edge length
        mesh.GetElementVertices(e, vert_ids);
        const int nv = vert_ids.Size();

        double h_min2 = std::numeric_limits<double>::infinity();

        for (int i = 0; i < nv; ++i)
        {
            const double *xi = mesh.GetVertex(vert_ids[i]);
            for (int j = i+1; j < nv; ++j)
            {
                const double *xj = mesh.GetVertex(vert_ids[j]);
                double dx2 = 0.0;
                for (int d = 0; d < dim; ++d)
                {
                    double diff = xi[d] - xj[d];
                    dx2 += diff * diff;
                }
                if (dx2 < h_min2) { h_min2 = dx2; }
            }
        }

        if (h_min2 == std::numeric_limits<double>::infinity())
        {
            h_min2 = 0.0;
        }

        adj.h[e] = std::sqrt(h_min2);
    }

    return adj;
}

bool FindElementForPointLocal(
    mfem::ParMesh            &mesh,
    const ElementAdjacency   &adj,
    int                       current_elem,
    const mfem::Vector       &x_new,
    int                      &out_elem,
    mfem::IntegrationPoint   &out_ip)
{
    using namespace mfem;

    // --- 1) Try current element ---
    {
        ElementTransformation *T = mesh.GetElementTransformation(current_elem);
        InverseElementTransformation invT(T);

        IntegrationPoint ip;
        int res = invT.Transform(x_new, ip);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = current_elem;
            out_ip   = ip;
            return true;
        }
    }

    // --- 2) Try neighbors ---
    const auto &nbrs = adj.neighbors[current_elem];
    for (int nb : nbrs)
    {
        if (nb < 0 || nb >= mesh.GetNE()) { continue; }

        ElementTransformation *T = mesh.GetElementTransformation(nb);
        InverseElementTransformation invT(T);

        IntegrationPoint ip;
        int res = invT.Transform(x_new, ip);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = nb;
            out_ip   = ip;
            return true;
        }
    }

    return false;
}

static bool FindElementForPointRobust(
    mfem::ParMesh            &mesh,
    const ElementAdjacency   &adj,
    int                       current_elem,
    const mfem::Vector       &x_new,
    int                      &out_elem,
    mfem::IntegrationPoint   &out_ip)
{
    using namespace mfem;
    // Need this in case we hit a vertex

    // First try cheap local search (current element + neighbors)
    if (FindElementForPointLocal(mesh, adj, current_elem, x_new, out_elem, out_ip))
    {
        return true;
    }

    // Fallback: global search over all elements.
    const int ne = mesh.GetNE();
    for (int e = 0; e < ne; ++e)
    {
        ElementTransformation *T = mesh.GetElementTransformation(e);
        InverseElementTransformation invT(T);

        IntegrationPoint ip;
        int res = invT.Transform(x_new, ip);
        if (res == InverseElementTransformation::Inside)
        {
            out_elem = e;
            out_ip   = ip;
            return true;
        }
    }

    // No element contains x_new anywhere in the mesh
    return false;
}

// Evaluate unit drift direction at an already-located (elem, ip).
// Assumes caller has ensured elem is valid and ip is inside that element.
// Uses the sign convention: k = E / |E| (as in your current code).
static void eval_k_in_element(mfem::ParMesh            &mesh,
                              ElectricFieldCoeff       &E_coeff,
                              int                       elem,
                              const mfem::IntegrationPoint &ip,
                              int                       dim,
                              // outputs:
                              mfem::Vector             &k_out,
                              bool                     &terminated,
                              ElectronExitCode         &term_code)
{
    using namespace mfem;

    if (terminated) { return; }

    ElementTransformation *T = mesh.GetElementTransformation(elem);
    T->SetIntPoint(&ip);

    Vector E(dim);
    E_coeff.Eval(E, *T, ip);

    const double Enorm = E.Norml2();
    if (Enorm <= 0.0)
    {
        terminated = true;
        term_code  = ElectronExitCode::DegenerateTimeStep;
        return;
    }

    k_out.SetSize(dim);
    k_out = E;
    k_out *= (1.0 / Enorm);
}

// Evaluate unit drift direction k at a physical point x_query.
// Does point-location first, then calls eval_k_in_element so all methods
// share the exact same field evaluation + normalization logic.
static void eval_k_at_point(mfem::ParMesh                 &mesh,
                            const ElementAdjacency        &adj,
                            ElectricFieldCoeff            &E_coeff,
                            const mfem::Vector            &x_query,
                            int                            search_elem,
                            int                            dim,
                            // outputs:
                            mfem::Vector                  &k_out,
                            int                            &e_out,
                            mfem::IntegrationPoint         &ip_out,
                            bool                           &terminated,
                            bool                           &need_shrink,
                            ElectronExitCode               &term_code)
{
    using namespace mfem;

    if (terminated || need_shrink) { return; }

    int e_tmp = -1;
    IntegrationPoint ip_tmp;

    if (!FindElementForPointRobust(mesh, adj, search_elem, x_query, e_tmp, ip_tmp) ||
        e_tmp < 0 || e_tmp >= mesh.GetNE())
    {
        need_shrink = true;
        return;
    }

    // Field evaluation + normalization centralized here:
    eval_k_in_element(mesh, E_coeff, e_tmp, ip_tmp, dim, k_out, terminated, term_code);
    if (terminated)
    {
        e_out  = e_tmp;
        ip_out = ip_tmp;
        return;
    }

    e_out  = e_tmp;
    ip_out = ip_tmp;
}

// -----------------------------------------------------------------
// Particle propagation algos
// All of it ChatGPT
// -----------------------------------------------------------------
static Vector transform_to_physical(Mesh&                   mesh,
                                    int                     elem_id,
                                    const IntegrationPoint& ip,
                                    int                     dim,
                                    bool                    axisymmetric,
                                    double                  geom_tol)
{
    IntegrationPoint ip_local = ip;

    ElementTransformation* T = mesh.GetElementTransformation(elem_id);
    T->SetIntPoint(&ip_local);

    Vector x(dim);
    T->Transform(ip_local, x);

    if (axisymmetric && x[0] <= geom_tol) {
        x[0] = geom_tol;
    }

    return x;
}

static void LocateEndpointOrShrink(mfem::ParMesh             &mesh,
                                   const ElementAdjacency    &adj,
                                   const mfem::Vector        &x_accept,
                                   int                        search_elem,
                                   // outputs:
                                   int                       &elem_new,
                                   mfem::IntegrationPoint    &ip_new,
                                   bool                      &need_shrink,
                                   bool                      &terminated)
{
    using namespace mfem;

    if (terminated || need_shrink) { return; }

    int e_tmp = -1;
    IntegrationPoint ip_tmp;

    if (!FindElementForPointRobust(mesh, adj, search_elem, x_accept, e_tmp, ip_tmp) ||
        e_tmp < 0 || e_tmp >= mesh.GetNE())
    {
        need_shrink = true;
        return;
    }

    elem_new = e_tmp;
    ip_new   = ip_tmp;
}

static void ClassifyTerminationAt(const TpcGeometry         &geom,
                                  const ElectronTraceParams &params,
                                  bool                       axisymmetric,
                                  const mfem::Vector        &x_accept,
                                  // outputs:
                                  bool                      &terminated,
                                  ElectronExitCode          &term_code)
{
    if (terminated) { return; }

    if (axisymmetric && params.terminate_on_axis && x_accept[0] <= params.geom_tol)
    {
        terminated = true;
        term_code  = ElectronExitCode::HitAxis;
        return;
    }

    // Assumes (r,z) stored in indices (0,1)
    const double r = x_accept[0];
    const double z = x_accept[1];

    const ElectronExitCode bc = geom.ClassifyBoundary(r, z, params.geom_tol);
    if (bc != ElectronExitCode::None)
    {
        terminated = true;
        term_code  = bc;
    }
}

static void AttemptStep_EulerCauchy(mfem::ParMesh                &mesh,
                                   ElectricFieldCoeff           &E_coeff,
                                   const ElementAdjacency       &adj,
                                   const TpcGeometry            &geom,
                                   const ElectronTraceParams    &params,
                                   bool                          axisymmetric,
                                   int                           elem_id,
                                   const mfem::IntegrationPoint &ip,
                                   double                        ds,
                                   // outputs:
                                   mfem::Vector                 &x_new,
                                   int                          &elem_new,
                                   mfem::IntegrationPoint       &ip_new,
                                   double                       &err_metric,
                                   bool                         &need_shrink,
                                   bool                         &terminated,
                                   ElectronExitCode             &term_code)
{
    using namespace mfem;

    const int dim = mesh.SpaceDimension();

    const Vector x0 = transform_to_physical(mesh, elem_id, ip, dim,
                                           axisymmetric, params.geom_tol);

    Vector k1(dim);
    eval_k_in_element(mesh, E_coeff, elem_id, ip, dim, k1, terminated, term_code);
    if (terminated)
    {
        x_new    = x0;
        elem_new = elem_id;
        ip_new   = ip;
        return;
    }

    Vector x_e(x0);
    x_e.Add(ds, k1);
    clamp_axis(x_e, axisymmetric, params.geom_tol);

    Vector k2(dim);
    int e_e = elem_id;
    IntegrationPoint ip_e = ip;

    eval_k_at_point(mesh, adj, E_coeff,
                    x_e, elem_id, dim,
                    k2, e_e, ip_e,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x_e;
        elem_new = e_e;
        ip_new   = ip_e;
        return;
    }

    Vector x_h(x0);
    x_h.Add(0.5 * ds, Vector(k1) += k2);
    clamp_axis(x_h, axisymmetric, params.geom_tol);

    err_metric = (Vector(x_h) -= x_e).Norml2();

    x_new = x_h;

    LocateEndpointOrShrink(mesh, adj, x_new, e_e,
                           elem_new, ip_new,
                           need_shrink, terminated);
    if (need_shrink || terminated)
    {
        if (need_shrink)
        {
            elem_new = e_e;
            ip_new   = ip_e;
        }
        return;
    }

    ClassifyTerminationAt(geom, params, axisymmetric, x_new,
                          terminated, term_code);
}



static void AttemptStep_RK23(mfem::ParMesh                &mesh,
                            ElectricFieldCoeff            &E_coeff,
                            const ElementAdjacency        &adj,
                            const TpcGeometry             &geom,
                            const ElectronTraceParams     &params,
                            bool                           axisymmetric,
                            int                            elem_id,
                            const mfem::IntegrationPoint  &ip,
                            double                         ds,
                            // outputs:
                            mfem::Vector                  &x_new,
                            int                           &elem_new,
                            mfem::IntegrationPoint        &ip_new,
                            double                        &err_metric,
                            bool                          &need_shrink,
                            bool                          &terminated,
                            ElectronExitCode              &term_code)
{
    using namespace mfem;

    const int dim = mesh.SpaceDimension();

    const Vector x0 = transform_to_physical(mesh, elem_id, ip, dim,
                                           axisymmetric, params.geom_tol);

    Vector k1(dim);
    eval_k_in_element(mesh, E_coeff, elem_id, ip, dim, k1, terminated, term_code);
    if (terminated)
    {
        x_new    = x0;
        elem_new = elem_id;
        ip_new   = ip;
        return;
    }

    // Stage 2: k2 at x1 = x0 + (ds/2)*k1
    Vector x1(x0);
    x1.Add(0.5 * ds, k1);
    clamp_axis(x1, axisymmetric, params.geom_tol);

    Vector k2(dim);
    int e1 = elem_id;
    IntegrationPoint ip1 = ip;
    eval_k_at_point(mesh, adj, E_coeff,
                    x1, elem_id, dim,
                    k2, e1, ip1,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x1;
        elem_new = e1;
        ip_new   = ip1;
        return;
    }

    // Stage 3: k3 at x2s = x0 + (3ds/4)*k2
    Vector x2s(x0);
    x2s.Add(0.75 * ds, k2);
    clamp_axis(x2s, axisymmetric, params.geom_tol);

    Vector k3(dim);
    int e2 = e1;
    IntegrationPoint ip2 = ip1;
    eval_k_at_point(mesh, adj, E_coeff,
                    x2s, e1, dim,
                    k3, e2, ip2,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x2s;
        elem_new = e2;
        ip_new   = ip2;
        return;
    }

    // Candidate (3rd order) endpoint: x3
    Vector x3(x0);
    x3.Add(ds * (2.0/9.0), k1);
    x3.Add(ds * (1.0/3.0), k2);
    x3.Add(ds * (4.0/9.0), k3);
    clamp_axis(x3, axisymmetric, params.geom_tol);

    // Stage 4: k4 at x3
    Vector k4(dim);
    int e3 = e2;
    IntegrationPoint ip3 = ip2;
    eval_k_at_point(mesh, adj, E_coeff,
                    x3, e2, dim,
                    k4, e3, ip3,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x3;
        elem_new = e3;
        ip_new   = ip3;
        return;
    }

    // 2nd order embedded estimate: x2
    Vector x2(x0);
    x2.Add(ds * (7.0/24.0), k1);
    x2.Add(ds * (1.0/4.0),  k2);
    x2.Add(ds * (1.0/3.0),  k3);
    x2.Add(ds * (1.0/8.0),  k4);
    clamp_axis(x2, axisymmetric, params.geom_tol);

    // Accept the 3rd-order solution candidate
    x_new = x3;

    // Error metric: ||x3 - x2||
    err_metric = (Vector(x3) -= x2).Norml2();

    // Locate accepted endpoint (common helper)
    LocateEndpointOrShrink(mesh, adj, x_new, e3,
                           elem_new, ip_new,
                           need_shrink, terminated);
    if (need_shrink || terminated)
    {
        if (need_shrink)
        {
            elem_new = e3;
            ip_new   = ip3;
        }
        return;
    }

    // Termination classification (common helper)
    ClassifyTerminationAt(geom, params, axisymmetric, x_new,
                          terminated, term_code);
}

static void AttemptStep_RK45(mfem::ParMesh                &mesh,
                            ElectricFieldCoeff            &E_coeff,
                            const ElementAdjacency        &adj,
                            const TpcGeometry             &geom,
                            const ElectronTraceParams     &params,
                            bool                           axisymmetric,
                            int                            elem_id,
                            const mfem::IntegrationPoint  &ip,
                            double                         ds,
                            // outputs:
                            mfem::Vector                  &x_new,
                            int                           &elem_new,
                            mfem::IntegrationPoint        &ip_new,
                            double                        &err_metric,
                            bool                          &need_shrink,
                            bool                          &terminated,
                            ElectronExitCode              &term_code)
{
    using namespace mfem;

    const int dim = mesh.SpaceDimension();

    const Vector x0 = transform_to_physical(mesh, elem_id, ip, dim,
                                           axisymmetric, params.geom_tol);

    Vector k1(dim);
    eval_k_in_element(mesh, E_coeff, elem_id, ip, dim, k1, terminated, term_code);
    if (terminated)
    {
        x_new    = x0;
        elem_new = elem_id;
        ip_new   = ip;
        return;
    }

    Vector k2(dim), k3(dim), k4(dim), k5(dim), k6(dim);

    // Stage 2: x1 = x0 + ds*(1/5)*k1
    Vector x1(x0);
    x1.Add(ds * (1.0/5.0), k1);
    clamp_axis(x1, axisymmetric, params.geom_tol);

    int e1 = elem_id;
    IntegrationPoint ip1 = ip;
    eval_k_at_point(mesh, adj, E_coeff,
                    x1, elem_id, dim,
                    k2, e1, ip1,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x1;
        elem_new = e1;
        ip_new   = ip1;
        return;
    }

    // Stage 3: x2 = x0 + ds*(3/40*k1 + 9/40*k2)
    Vector x2(x0);
    x2.Add(ds * (3.0/40.0), k1);
    x2.Add(ds * (9.0/40.0), k2);
    clamp_axis(x2, axisymmetric, params.geom_tol);

    int e2 = e1;
    IntegrationPoint ip2 = ip1;
    eval_k_at_point(mesh, adj, E_coeff,
                    x2, e1, dim,
                    k3, e2, ip2,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x2;
        elem_new = e2;
        ip_new   = ip2;
        return;
    }

    // Stage 4: x3s = x0 + ds*(3/10*k1 - 9/10*k2 + 6/5*k3)
    Vector x3s(x0);
    x3s.Add(ds * (3.0/10.0),  k1);
    x3s.Add(ds * (-9.0/10.0), k2);
    x3s.Add(ds * (6.0/5.0),   k3);
    clamp_axis(x3s, axisymmetric, params.geom_tol);

    int e3 = e2;
    IntegrationPoint ip3 = ip2;
    eval_k_at_point(mesh, adj, E_coeff,
                    x3s, e2, dim,
                    k4, e3, ip3,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x3s;
        elem_new = e3;
        ip_new   = ip3;
        return;
    }

    // Stage 5: x4s = x0 + ds*(-11/54*k1 + 5/2*k2 -70/27*k3 +35/27*k4)
    Vector x4s(x0);
    x4s.Add(ds * (-11.0/54.0), k1);
    x4s.Add(ds * (5.0/2.0),    k2);
    x4s.Add(ds * (-70.0/27.0), k3);
    x4s.Add(ds * (35.0/27.0),  k4);
    clamp_axis(x4s, axisymmetric, params.geom_tol);

    int e4 = e3;
    IntegrationPoint ip4 = ip3;
    eval_k_at_point(mesh, adj, E_coeff,
                    x4s, e3, dim,
                    k5, e4, ip4,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x4s;
        elem_new = e4;
        ip_new   = ip4;
        return;
    }

    // Stage 6: x5s = x0 + ds*(1631/55296*k1 + 175/512*k2 + 575/13824*k3
    //                        + 44275/110592*k4 + 253/4096*k5)
    Vector x5s(x0);
    x5s.Add(ds * (1631.0/55296.0),   k1);
    x5s.Add(ds * (175.0/512.0),      k2);
    x5s.Add(ds * (575.0/13824.0),    k3);
    x5s.Add(ds * (44275.0/110592.0), k4);
    x5s.Add(ds * (253.0/4096.0),     k5);
    clamp_axis(x5s, axisymmetric, params.geom_tol);

    int e5 = e4;
    IntegrationPoint ip5 = ip4;
    eval_k_at_point(mesh, adj, E_coeff,
                    x5s, e4, dim,
                    k6, e5, ip5,
                    terminated, need_shrink, term_code);
    if (terminated || need_shrink)
    {
        x_new    = x5s;
        elem_new = e5;
        ip_new   = ip5;
        return;
    }

    // 5th order (accepted): x5
    Vector x5(x0);
    x5.Add(ds * (37.0/378.0),   k1);
    x5.Add(ds * (250.0/621.0),  k3);
    x5.Add(ds * (125.0/594.0),  k4);
    x5.Add(ds * (512.0/1771.0), k6);
    clamp_axis(x5, axisymmetric, params.geom_tol);

    // 4th order (embedded): x4
    Vector x4(x0);
    x4.Add(ds * (2825.0/27648.0),  k1);
    x4.Add(ds * (18575.0/48384.0), k3);
    x4.Add(ds * (13525.0/55296.0), k4);
    x4.Add(ds * (277.0/14336.0),   k5);
    x4.Add(ds * (1.0/4.0),         k6);
    clamp_axis(x4, axisymmetric, params.geom_tol);

    err_metric = (Vector(x5) -= x4).Norml2();

    x_new = x5;

    LocateEndpointOrShrink(mesh, adj, x_new, e5,
                           elem_new, ip_new,
                           need_shrink, terminated);
    if (need_shrink || terminated)
    {
        if (need_shrink)
        {
            elem_new = e5;
            ip_new   = ip5;
        }
        return;
    }

    ClassifyTerminationAt(geom, params, axisymmetric, x_new,
                          terminated, term_code);
}


static void AttemptStep_Debug(mfem::ParMesh                &mesh,
                              ElectricFieldCoeff           &E_coeff,
                              const ElementAdjacency       &adj,
                              const TpcGeometry            &geom,
                              const ElectronTraceParams    &params,
                              bool                          axisymmetric,
                              int                           elem_id,
                              const mfem::IntegrationPoint &ip,
                              double                        ds,
                              // outputs:
                              mfem::Vector                &x_new,
                              int                          &elem_new,
                              mfem::IntegrationPoint       &ip_new,
                              double                       &err_metric,
                              bool                         &need_shrink,
                              bool                         &terminated,
                              ElectronExitCode             &term_code)
{
    using namespace mfem;

    const int dim = mesh.SpaceDimension();

    const Vector x0 = transform_to_physical(mesh, elem_id, ip, dim,
                                           axisymmetric, params.geom_tol);

    Vector k1(dim);
    eval_k_in_element(mesh, E_coeff, elem_id, ip, dim, k1, terminated, term_code);
    if (terminated)
    {
        x_new    = x0;
        elem_new = elem_id;
        ip_new   = ip;
        return;
    }

    x_new = x0;
    x_new.Add(ds, k1);
    clamp_axis(x_new, axisymmetric, params.geom_tol);

    err_metric = 0.0;

    LocateEndpointOrShrink(mesh, adj, x_new, elem_id,
                           elem_new, ip_new,
                           need_shrink, terminated);
    if (need_shrink || terminated)
    {
        if (need_shrink)
        {
            elem_new = elem_id;
            ip_new   = ip;
        }
        return;
    }

    ClassifyTerminationAt(geom, params, axisymmetric, x_new,
                          terminated, term_code);
}



// Trace Electron Line 
ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh                &mesh,
    ElectricFieldCoeff           &E_coeff,
    const ElementAdjacency       &adj,
    const TpcGeometry            &geom,
    const ElectronTraceParams    &params,
    int                           start_elem,
    const mfem::IntegrationPoint &start_ip,
    bool                          axisymmetric,
    bool                          save_pathlines)
{
    using namespace mfem;

    ElectronTraceResult result;
    result.exit_code    = ElectronExitCode::None;
    result.exit_element = start_elem;

    const int dim = mesh.SpaceDimension();

    int              elem_id = start_elem;
    IntegrationPoint ip      = start_ip;

    // Get Step Function
    StepFunction StepFunc = nullptr;
    if      (params.method == "Euler-Cauchy"){ StepFunc = &AttemptStep_EulerCauchy; }
    else if (params.method == "RK23"){         StepFunc = &AttemptStep_RK23; }
    else if (params.method == "RK45"){         StepFunc = &AttemptStep_RK45; }
    else if (params.method == "Debug"){        StepFunc = &AttemptStep_Debug; }
    else    { throw std::invalid_argument( "Integration Method " + params.method + " unknown.\n"); }
    Vector x(dim);

    // Initial point
    ElementTransformation *T0 = mesh.GetElementTransformation(elem_id);
    T0->SetIntPoint(&ip);
    T0->Transform(ip, x);

    if (axisymmetric && x[0] <= params.geom_tol) 
    { x[0] = params.geom_tol; }

    if (save_pathlines) { result.points.push_back(x); }

    int    step      = 0;
    double c_current = 0.0; // adaptive *dimensionless* step for this electron

    while (step < params.max_steps)
    {
        // Local element size
        const double hK = adj.h[elem_id];

        // Per-element bounds in c-space
        const double c_min = params.c_min_factor;
        const double c_max = params.c_max_factor;

        if (c_current <= 0.0) { c_current = params.c_step; }
        c_current = std::min(std::max(c_current, c_min), c_max);

        bool accepted_step = false;

        while (!accepted_step)
        {
            const double ds = c_current * hK;

            // ---------------------------------------------------------------------
            // Trial step (delegated by integrator choice)
            // ---------------------------------------------------------------------
            bool             need_shrink = false;
            bool             terminated  = false;
            ElectronExitCode term_code   = ElectronExitCode::None;

            Vector           x_new(dim);
            int              elem_new = elem_id;
            IntegrationPoint ip_new  = ip;

            double err_metric = 0.0;

            StepFunc(mesh, E_coeff, adj, geom, params, axisymmetric, elem_id, ip, ds,
                    // outputs:
                    x_new, elem_new, ip_new, err_metric, need_shrink, terminated, term_code);

            // ---------------------------------------------------------------------
            // Handle termination / topology shrink
            // ---------------------------------------------------------------------
            if (terminated)
            {
                result.exit_code    = term_code;
                result.exit_element = elem_new; // last known valid in-volume element
                if (save_pathlines) { result.points.push_back(x_new); }
                return result;
            }

            if (need_shrink)
            {
                const double c_new = c_current * params.adapt_shrink;
                if (c_new < c_min)
                {
                    result.exit_code    = ElectronExitCode::LeftVolume;
                    result.exit_element = elem_id;
                    return result;
                }
                c_current = c_new;
                continue; // retry
            }

            // ---------------------------------------------------------------------
            // Error-based acceptance (common policy)
            // ---------------------------------------------------------------------
            const double tol = params.tol_rel;

            if (err_metric > tol && c_current > c_min)
            {
                double c_new = c_current * params.adapt_shrink;
                if (c_new < c_min) { c_new = c_min; }
                c_current = c_new;
                continue; // retry
            }

            // Accept step
            elem_id = elem_new;
            ip      = ip_new;
            x       = x_new;

            if (save_pathlines) { result.points.push_back(x_new); }
            result.exit_element = elem_id;
            accepted_step = true;

            // Optional step growth in smooth regions
            if (err_metric < 0.25 * tol)
            {
                double c_new = c_current * params.adapt_grow;
                if (c_new > c_max) { c_new = c_max; }
                c_current = c_new;
            }
        } // end adaptive trial loop

        ++step;
    }

    if (result.exit_code == ElectronExitCode::None)
    {
        result.exit_code = ElectronExitCode::MaxSteps;
    }

    return result;
}
// -----------------------
// Parallel interface, first a bunch of helpers
// -----------------------
static void DebugFindCentreSeed(mfem::ParMesh            &global_mesh,
                                int                      dim,
                                const ElectronTraceParams &params,
                                const Config             &cfg,
                                int                      &elem_center,
                                mfem::IntegrationPoint   &ip_center)
{
    using namespace mfem;

    // 1) Construct centre point in physical coordinates
    const double r_c = 0.5 * (params.r_min + params.r_max);
    const double z_c = 0.5 * (params.z_min + params.z_max);

    mfem::Vector x_center(dim);
    x_center[0] = r_c;
    x_center[1] = z_c;

    // 2) Find element and reference IntegrationPoint that contain x_center
    elem_center = -1;

    for (int e = 0; e < global_mesh.GetNE(); ++e)
    {
        mfem::ElementTransformation *T = global_mesh.GetElementTransformation(e);
        mfem::InverseElementTransformation invT(T);

        mfem::IntegrationPoint ip_trial;
        int res = invT.Transform(x_center, ip_trial);
        if (res == InverseElementTransformation::Inside)
        {
            elem_center = e;
            ip_center   = ip_trial;
            break;
        }
    }

    MFEM_VERIFY(elem_center >= 0,
                "TraceElectronFieldLines: centre point not found in any element");

    mfem::ElementTransformation *Tcheck =
        global_mesh.GetElementTransformation(elem_center);
    Tcheck->SetIntPoint(&ip_center);

    mfem::Vector x_check(dim);
    Tcheck->Transform(ip_center, x_check);

    std::cout << "[DEBUG:TRACE] Centre from config : r=" << r_c
              << ", z=" << z_c << "\n";
    std::cout << "[DEBUG:TRACE] Mapped back centre : r=" << x_check[0]
              << ", z=" << x_check[1] << "\n";

    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:TRACE] Single-seed debug mode active.\n"
                  << "               Centre point: r=" << r_c
                  << ", z=" << z_c << "\n"
                  << "               Element id : " << elem_center << "\n";
    }
}

static void DebugSingleSeedTrace(mfem::ParMesh                       &global_mesh,
                                 int                                 dim,
                                 mfem::GridFunction                 &global_phi,
                                 const mfem::FiniteElementCollection *fec_phi,
                                 int                                 ordering_phi,
                                 const ElementAdjacency             &adj,
                                 const TpcGeometry                  &geom,
                                 const ElectronTraceParams          &params,
                                 const Config                       &cfg,
                                 std::vector<ElectronTraceResult>   &out_results,
                                 bool                                axisymmetric,
                                 bool                                save_paths)
{
    using namespace mfem;

    int elem_center = -1;
    mfem::IntegrationPoint ip_center;

    DebugFindCentreSeed(global_mesh,
                        dim,
                        params,
                        cfg,
                        elem_center,
                        ip_center);

    // 3) Local mesh and potential
    mfem::ParMesh local_mesh(global_mesh);

    mfem::ParFiniteElementSpace local_fes_phi(
        &local_mesh, fec_phi, /*vdim=*/1, ordering_phi);

    mfem::GridFunction local_phi(&local_fes_phi);
    local_phi = global_phi; // copy potential DOFs

    // Local E coefficient: E = -grad(phi)
    ElectricFieldCoeff local_E_coeff(local_phi, 1.0);

    const ElementAdjacency &local_adj = adj;
    for (std::size_t i = 0; i < local_adj.h.size(); ++i)
    {
        if (local_adj.h[i] <= 0.0)
        {
            throw std::runtime_error("Determind Mesh size non positive");
        }
    }

    // 4) Local params
    ElectronTraceParams local_params = params;

    // 5) Resize outputs and trace exactly one electron
    out_results.clear();
    out_results.resize(1);

    out_results[0] = TraceSingleElectronLine(local_mesh,
                                             local_E_coeff,
                                             local_adj,
                                             geom,
                                             local_params,
                                             elem_center,
                                             ip_center,
                                             axisymmetric,
                                             save_paths);
}

static void PrintExitConditionSummary(const Config                     &cfg,
                                      const CivSeeds                   &seeds,
                                      const std::vector<ElectronTraceResult> &out_results)
{
    if (!cfg.debug.debug) { return; }

    std::size_t count_hit_lgi      = 0;
    std::size_t count_hit_cathode  = 0;
    std::size_t count_hit_wall     = 0;
    std::size_t count_left_volume  = 0;
    std::size_t count_max_steps    = 0;
    std::size_t count_deg_dt       = 0;
    std::size_t count_none         = 0;
    std::size_t count_hit_axis     = 0;

    for (const auto &res : out_results)
    {
        switch (res.exit_code)
        {
            case ElectronExitCode::HitLiquidGas:       count_hit_lgi++;     break;
            case ElectronExitCode::HitCathode:         count_hit_cathode++; break;
            case ElectronExitCode::HitWall:            count_hit_wall++;    break;
            case ElectronExitCode::LeftVolume:         count_left_volume++; break;
            case ElectronExitCode::MaxSteps:           count_max_steps++;   break;
            case ElectronExitCode::DegenerateTimeStep: count_deg_dt++;      break;
            case ElectronExitCode::HitAxis:            count_hit_axis++;    break;
            case ElectronExitCode::None:               count_none++;        break;
        }
    }

    const double N = static_cast<double>(out_results.size());
    auto frac = [&](std::size_t c) -> double {
        return (N > 0.0 ? static_cast<double>(c) / N : 0.0);
    };

    double V_total        = 0.0;
    double V_hit_lgi      = 0.0;
    double V_hit_cathode  = 0.0;
    double V_hit_wall     = 0.0;
    double V_left_volume  = 0.0;
    double V_max_steps    = 0.0;
    double V_deg_dt       = 0.0;
    double V_none         = 0.0;
    double V_hit_axis     = 0.0;

    for (std::size_t i = 0; i < seeds.positions.size(); ++i)
    {
        const double dV = seeds.volumes[i];
        const auto  &res = out_results[i];

        V_total += dV;

        switch (res.exit_code)
        {
            case ElectronExitCode::HitLiquidGas:       V_hit_lgi     += dV; break;
            case ElectronExitCode::HitCathode:         V_hit_cathode += dV; break;
            case ElectronExitCode::HitWall:            V_hit_wall    += dV; break;
            case ElectronExitCode::LeftVolume:         V_left_volume += dV; break;
            case ElectronExitCode::MaxSteps:           V_max_steps   += dV; break;
            case ElectronExitCode::DegenerateTimeStep: V_deg_dt      += dV; break;
            case ElectronExitCode::HitAxis:            V_hit_axis    += dV; break;
            case ElectronExitCode::None:               V_none        += dV; break;
        }
    }

    auto vfrac = [&](double V) -> double {
        return (V_total > 0.0 ? V / V_total : 0.0);
    };

    std::cout << "\n---------------- Electron Tracing Summary ----------------\n";
    std::cout << "Total seeds traced: " << out_results.size() << "\n\n";

    std::cout << std::left
              << std::setw(22) << "Exit Type"
              << std::setw(18) << "Seed Fraction"
              << std::setw(18) << "Volume Fraction\n";

    std::cout << "-----------------------------------------------------------------\n";

    auto row = [&](const char *label, double f_seed, double f_vol) {
        std::cout << std::left
                  << std::setw(22) << label
                  << std::setw(18) << f_seed
                  << std::setw(18) << f_vol
                  << "\n";
    };

    row("Hit Liquid-Gas",   frac(count_hit_lgi),     vfrac(V_hit_lgi));
    row("Hit Cathode",      frac(count_hit_cathode), vfrac(V_hit_cathode));
    row("Hit Wall",         frac(count_hit_wall),    vfrac(V_hit_wall));
    row("Hit Axis",         frac(count_hit_axis),    vfrac(V_hit_axis));
    row("Left Volume",      frac(count_left_volume), vfrac(V_left_volume));
    row("Max Steps",        frac(count_max_steps),   vfrac(V_max_steps));
    row("Degenerate dt",    frac(count_deg_dt),      vfrac(V_deg_dt));
    row("None (unexpected)",frac(count_none),        vfrac(V_none));

    std::cout << "-----------------------------------------------------------------\n";
}

static void DumpElectronPathsCSV(const Config                          &cfg,
                                 const std::vector<ElectronTraceResult> &out_results)
{
    if (!cfg.debug.dumpdata) { return; }

    std::cout << "[DEBUG] Dumping Electron Paths in CIV" << std::endl;
    namespace fs = std::filesystem;
    fs::path outdir(cfg.save_path);
    fs::path outfile = outdir / "electron_paths_debug.csv";
    std::ofstream ofs(outfile);
    ofs << "# id, step, r, z, exit_code\n";

    for (std::size_t i = 0; i < out_results.size(); ++i)
    {
        const auto &res = out_results[i];
        const auto &pts = res.points;

        for (std::size_t k = 0; k < pts.size(); ++k)
        {
            const mfem::Vector &x = pts[k];
            const double r = x(0);
            const double z = x(1);

            ofs << i << "," << k << "," << r << "," << z << ","
                << static_cast<int>(res.exit_code) << "\n";
        }
    }
}


// -----------------------------------------------------------------------------
// INNER: does EVERYTHING except object creation needed to satisfy its signature
// -----------------------------------------------------------------------------
void TraceElectronFieldLinesInner(mfem::ParMesh                   &global_mesh,
                                  mfem::GridFunction              &global_phi,
                                  const ElementAdjacency          &adj,
                                  const mfem::FiniteElementCollection *fec_phi,
                                  int                              ordering_phi,
                                  const ElectronTraceParams        &params,
                                  const Config                     &cfg,
                                  const CivSeeds                   &seeds,
                                  std::vector<ElectronTraceResult> &out_results,
                                  bool                             axisymmetric,
                                  bool                             save_paths,
                                  const double                     *z_max_overrides) // size = n_seeds or nullptr
{
    using namespace mfem;

    const int dim = global_mesh.SpaceDimension();
    const std::size_t n_seeds = seeds.positions.size();

    if (z_max_overrides != nullptr)
    {
        MFEM_VERIFY(out_results.size() == n_seeds,
                    "TraceElectronFieldLinesInner: out_results must be sized to n_seeds");
    }

    // -------------------------------------------------------------------------
    // Single-seed debug mode
    // -------------------------------------------------------------------------
    if (cfg.debug.debug_single_seed)
    {
        // For debug_single_seed there isn't a natural seed index; use [0] if provided.
        ElectronTraceParams local_params = params;
        if (z_max_overrides != nullptr && n_seeds > 0)
        {
            local_params.z_max = z_max_overrides[0];
        }

        TpcGeometry geom(local_params);

        DebugSingleSeedTrace(global_mesh,
                             dim,
                             global_phi,
                             fec_phi,
                             ordering_phi,
                             adj,
                             geom,
                             local_params,
                             cfg,
                             out_results,
                             axisymmetric,
                             save_paths);
    }
    else
    {
        // ---------------------------------------------------------------------
        // Normal mode: trace all seeds in parallel
        // ---------------------------------------------------------------------
        #pragma omp parallel
        {
            ParMesh local_mesh(global_mesh);

            ParFiniteElementSpace local_fes_phi(&local_mesh, fec_phi, /*vdim=*/1, ordering_phi);
            GridFunction local_phi(&local_fes_phi);
            local_phi = global_phi;

            ElectricFieldCoeff local_E_coeff(local_phi, 1.0);

            const ElementAdjacency &local_adj = adj;

            #pragma omp for schedule(static)
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(n_seeds); ++i)
            {
                const std::size_t ui = static_cast<std::size_t>(i);
                const int start_elem = seeds.elements[ui];
                const IntegrationPoint &start_ip = seeds.ips[ui];

                ElectronTraceParams local_params = params;
                if (z_max_overrides != nullptr)
                {
                    local_params.z_max = z_max_overrides[ui];
                }

                // Per-seed geometry so ClassifyBoundary sees the per-seed z_max
                TpcGeometry local_geom(local_params);

                ElectronTraceResult res = TraceSingleElectronLine(local_mesh,
                                                                  local_E_coeff,
                                                                  local_adj,
                                                                  local_geom,
                                                                  local_params,
                                                                  start_elem,
                                                                  start_ip,
                                                                  axisymmetric,
                                                                  save_paths);

                out_results[ui] = std::move(res);
            }
        } // end parallel
    }

    PrintExitConditionSummary(cfg, seeds, out_results);
    DumpElectronPathsCSV(cfg, out_results);
}

// -----------------------------------------------------------------------------
// OUTER: only constructs objects needed for the inner signature, then calls it
// (nonstatic free function)
// -----------------------------------------------------------------------------
void TraceElectronFieldLines(const SimulationResult           &sim,
                             const Config                     &cfg,
                             const CivSeeds                   &seeds,
                             std::vector<ElectronTraceResult> &out_results)
{
    using namespace mfem;

    ElectronTraceParams params = cfg.tracing_params;
    params.z_max = cfg.tracing_params.tracing_z_max;

    ParMesh &global_mesh = *sim.mesh;
    GridFunction &global_phi = *sim.V;

    const std::size_t n_seeds = seeds.positions.size();
    MFEM_VERIFY(seeds.elements.size() == n_seeds &&
                seeds.ips.size()      == n_seeds,
                "TraceElectronFieldLines: CivSeeds size mismatch");

    out_results.clear();
    out_results.resize(n_seeds);

    ElementAdjacency adj = BuildAdjacency(global_mesh);

    TpcGeometry geom(params);

    const FiniteElementSpace      *fes_phi      = global_phi.FESpace();
    const FiniteElementCollection *fec_phi      = fes_phi->FEColl();
    const int                      ordering_phi = fes_phi->GetOrdering();

    const bool axisymmetric = cfg.solver.axisymmetric;
    const bool save_paths   = cfg.debug.dumpdata;

    TraceElectronFieldLinesInner(global_mesh,
                                 global_phi,
                                 adj,
                                 fec_phi,
                                 ordering_phi,
                                 params,
                                 cfg,
                                 seeds,
                                 out_results,
                                 axisymmetric,
                                 save_paths,
                                /*z_max_overrides=*/nullptr);
}
