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
#include <stdexcept>

#include "CIV.h"


// Return true for any form of charge insensitivity (cathode or wall)
static bool IsChargeInsensitive(const ElectronTraceResult &res)
{
    return (res.exit_code == ElectronExitCode::HitWall);
}

static Seeds ExtractCivSeeds(const Config &cfg, const SimulationResult &result)
{
    using namespace mfem;
    const int sdim = result.mesh->SpaceDimension();
    const auto &tp = cfg.tracing_params;
    const double r0 = tp.r_min + cfg.tracing_params.geom_tol, r1 = tp.r_max - cfg.tracing_params.geom_tol;
    const double z0 = tp.z_min + cfg.tracing_params.geom_tol, z1 = tp.z_max - cfg.tracing_params.geom_tol;

    int n = cfg.civ_params.num_seed_elements;

    std::mt19937 rng;
    if (cfg.civ_params.rng_seed != 0)
    {
        rng.seed(static_cast<std::mt19937::result_type>(cfg.civ_params.rng_seed));
    }
    else
    {
        std::random_device rd;
        rng.seed(rd());
    }

    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const bool axisymmetric = cfg.solver.axisymmetric;
    // Total ROI "volume" used for weights.
    const double V_total = (!axisymmetric)
        ? (r1 - r0) * (z1 - z0)
        : (M_PI * (r1 * r1 - r0 * r0) * (z1 - z0)); // ∫∫ 2π r dr dz

    const double w = (n > 0) ? (V_total / static_cast<double>(n)) : 0.0;

    Seeds seeds;
    seeds.positions.reserve(static_cast<std::size_t>(n));
    seeds.volumes.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        const double uz = U01(rng);
        const double z  = z0 + uz * (z1 - z0);
        double r = 0.0;
        if (!axisymmetric)
        {
            const double ur = U01(rng);
            r = r0 + ur * (r1 - r0);
        }
        else
        {
            // Sample r with pdf(r) ∝ r on [r0, r1] for uniform cylindrical volume.
            const double u = U01(rng);
            const double r2 = (r0 * r0) + u * ((r1 * r1) - (r0 * r0));
            r = std::sqrt(r2);
        }
        Vector p(sdim);
        p = 0.0;
        p[0] = r;
        p[1] = z;
        seeds.positions.push_back(std::move(p));
        seeds.volumes.push_back(w);
    }
    return seeds;
}

static Seeds ExtractCivFixedGridSeeds(const Config            &cfg,
                                  const SimulationResult &result)
{
    using namespace mfem;

    if (!result.mesh)
    {
        throw std::runtime_error("ExtractCivFixedGridSeeds: result.mesh is null");
    }

    const int sdim = result.mesh->SpaceDimension();
    if (sdim < 2)
    {
        throw std::runtime_error("ExtractCivFixedGridSeeds: mesh.SpaceDimension() < 2");
    }

    const auto &tp = cfg.tracing_params;

    const double r0 = tp.r_min;
    const double r1 = tp.r_max;
    const double z0 = tp.z_min;
    const double z1 = tp.z_max;

    const int nr = cfg.civ_params.nr;
    const int nz = cfg.civ_params.nz;

    if (nr <= 0 || nz <= 0)
    {
        throw std::runtime_error("ExtractCivFixedGridSeeds: cfg.civ_params.nr/nz must be > 0");
    }

    const double dr = (r1 - r0) / static_cast<double>(nr);
    const double dz = (z1 - z0) / static_cast<double>(nz);

    Seeds seeds;
    seeds.positions.resize(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));
    seeds.volumes.resize(seeds.positions.size());

    const bool axisymmetric = cfg.solver.axisymmetric;

    // Cell-center sampling, row-major indexing: k = iz*nr + ir
    for (int iz = 0; iz < nz; ++iz)
    {
        const double z = z0 + (static_cast<double>(iz) + 0.5) * dz;

        for (int ir = 0; ir < nr; ++ir)
        {
            const double r = r0 + (static_cast<double>(ir) + 0.5) * dr;

            const std::size_t k = static_cast<std::size_t>(iz) * static_cast<std::size_t>(nr)
                                + static_cast<std::size_t>(ir);

            Vector p(sdim);
            p = 0.0;
            p[0] = r;
            p[1] = z;
            // If sdim == 3, p[2] stays 0.0; fixed-grid CIV is intended for (r,z).
            seeds.positions[k] = std::move(p);

            double dV = dr * dz;
            if (axisymmetric)
            {
                dV *= (2.0 * M_PI * r);
            }
            seeds.volumes[k] = dV;
        }
    }

    return seeds;
}

static Seeds ExtractCivSeeds_Row(const Config &cfg,
                            const SimulationResult &result,
                            int iz_row,
                            int ir_begin,
                            int ir_end_inclusive)
{
    using namespace mfem;

    if (!result.mesh)
    {
        throw std::runtime_error("ExtractCivSeeds_Row: result.mesh is null");
    }

    const int sdim = result.mesh->SpaceDimension();
    if (sdim < 2)
    {
        throw std::runtime_error("ExtractCivSeeds_Row: mesh.SpaceDimension() < 2");
    }

    const int nr = cfg.civ_params.nr;
    const int nz = cfg.civ_params.nz;

    if (nr <= 0 || nz <= 0)
    {
        throw std::runtime_error("ExtractCivSeeds_Row: cfg.civ_params.nr/nz must be > 0");
    }
    if (iz_row < 0 || iz_row >= nz)
    {
        throw std::runtime_error("ExtractCivSeeds_Row: iz_row out of range");
    }

    const auto &tp = cfg.tracing_params;
    const double r0 = tp.r_min, r1 = tp.r_max;
    const double z0 = tp.z_min, z1 = tp.z_max;

    if (!(r1 > r0) || !(z1 > z0))
    {
        throw std::runtime_error("ExtractCivSeeds_Row: invalid bounds (r/z)");
    }

    const double dr = (r1 - r0) / static_cast<double>(nr);
    const double dz = (z1 - z0) / static_cast<double>(nz);

    // Normalize segment bounds to a list of ir values.
    const int step = (ir_begin <= ir_end_inclusive) ? 1 : -1;
    const int nseg = std::abs(ir_end_inclusive - ir_begin) + 1;

    if (ir_begin < 0 || ir_begin >= nr || ir_end_inclusive < 0 || ir_end_inclusive >= nr)
    {
        throw std::runtime_error("ExtractCivSeeds_Row: ir range out of bounds");
    }

    Seeds seeds;
    seeds.positions.resize(static_cast<std::size_t>(nseg));
    seeds.volumes.resize(seeds.positions.size());

    const bool axisymmetric = cfg.solver.axisymmetric;

    const double z = z0 + (static_cast<double>(iz_row) + 0.5) * dz;

    int ir = ir_begin;
    for (int k = 0; k < nseg; ++k, ir += step)
    {
        const double r = r0 + (static_cast<double>(ir) + 0.5) * dr;

        Vector p(sdim);
        p = 0.0;
        p[0] = r;
        p[1] = z;

        seeds.positions[static_cast<std::size_t>(k)] = std::move(p);

        double dV = dr * dz;
        if (axisymmetric)
        {
            dV *= (2.0 * M_PI * r);
        }
        seeds.volumes[static_cast<std::size_t>(k)] = dV;
    }

    return seeds;
}

// ------------------- Adaptive Sampling Helpers --------------
Seeds ExtractCivSeeds_Column(const Config &cfg,
                               const SimulationResult &result,
                               int ir_col)
{
    using namespace mfem;

    if (!result.mesh)
    {
        throw std::runtime_error("ExtractCivSeeds_Column: result.mesh is null");
    }

    const int sdim = result.mesh->SpaceDimension();
    if (sdim < 2)
    {
        throw std::runtime_error("ExtractCivSeeds_Column: mesh.SpaceDimension() < 2");
    }

    const int nr = cfg.civ_params.nr;
    const int nz = cfg.civ_params.nz;

    if (nr <= 0 || nz <= 0)
    {
        throw std::runtime_error("ExtractCivSeeds_Column: cfg.civ_params.nr/nz must be > 0");
    }
    if (ir_col < 0 || ir_col >= nr)
    {
        throw std::runtime_error("ExtractCivSeeds_Column: ir_col out of range");
    }

    const auto &tp = cfg.tracing_params;
    const double r0 = tp.r_min, r1 = tp.r_max;
    const double z0 = tp.z_min, z1 = tp.z_max;

    if (!(r1 > r0) || !(z1 > z0))
    {
        throw std::runtime_error("ExtractCivSeeds_Column: invalid bounds (r/z)");
    }

    const double dr = (r1 - r0) / static_cast<double>(nr);
    const double dz = (z1 - z0) / static_cast<double>(nz);

    const double r = r0 + (static_cast<double>(ir_col) + 0.5) * dr;

    Seeds seeds;
    seeds.positions.resize(static_cast<std::size_t>(nz));
    seeds.volumes.resize(seeds.positions.size());

    const bool axisymmetric = cfg.solver.axisymmetric;

    for (int iz = 0; iz < nz; ++iz)
    {
        const double z = z0 + (static_cast<double>(iz) + 0.5) * dz;

        Vector p(sdim);
        p = 0.0;
        p[0] = r;
        p[1] = z;
        // If sdim==3, p[2]=0 remains.

        seeds.positions[static_cast<std::size_t>(iz)] = std::move(p);

        double dV = dr * dz;
        if (axisymmetric)
        {
            dV *= (2.0 * M_PI * r);
        }
        seeds.volumes[static_cast<std::size_t>(iz)] = dV;
    }

    return seeds;
}

void UpdateLeftmostCivFromColumn(const Config &cfg,
                                 const SimulationResult & /*result*/,
                                 int ir_col,
                                 const Seeds &col_seeds,
                                 const std::vector<ElectronTraceResult> &col_results,
                                 std::vector<int> &leftmost_civ_per_row,
                                 bool &column_has_any_civ)
{
    const int nz = cfg.civ_params.nz;

    if (nz <= 0)
    {
        throw std::runtime_error("UpdateLeftmostCivFromColumn: cfg.civ_params.nz must be > 0");
    }
    if (static_cast<int>(leftmost_civ_per_row.size()) != nz)
    {
        throw std::runtime_error("UpdateLeftmostCivFromColumn: leftmost_civ_per_row has wrong size");
    }
    if (col_seeds.positions.size() != static_cast<std::size_t>(nz) ||
        col_results.size() != static_cast<std::size_t>(nz))
    {
        throw std::runtime_error("UpdateLeftmostCivFromColumn: column seeds/results size mismatch");
    }

    column_has_any_civ = false;

    for (int iz = 0; iz < nz; ++iz)
    {
        const auto &res = col_results[static_cast<std::size_t>(iz)];

        if (IsChargeInsensitive(res))
        {
            column_has_any_civ = true;

            int &lm = leftmost_civ_per_row[static_cast<std::size_t>(iz)];
            if (lm < 0 || ir_col < lm)
            {
                lm = ir_col;
            }
        }
    }
}


Seeds ExtractCivSeeds_BoundaryBand(const Config &cfg,
                                     const SimulationResult &result,
                                     int level,
                                     const std::vector<int> &leftmost_civ_prev_level)
{
    using namespace mfem;

    if (level <= 0)
    {
        throw std::runtime_error("ExtractCivSeeds_BoundaryBand: level must be >= 1");
    }
    if (!result.mesh)
    {
        throw std::runtime_error("ExtractCivSeeds_BoundaryBand: result.mesh is null");
    }

    const int sdim = result.mesh->SpaceDimension();
    if (sdim < 2)
    {
        throw std::runtime_error("ExtractCivSeeds_BoundaryBand: mesh.SpaceDimension() < 2");
    }

    const int nr0 = cfg.civ_params.nr;
    const int nz0 = cfg.civ_params.nz;
    if (nr0 <= 0 || nz0 <= 0)
    {
        throw std::runtime_error("ExtractCivSeeds_BoundaryBand: cfg.civ_params.nr/nz must be > 0");
    }

    // Refined dimensions at this level.
    const int scale = 1 << level;       // 2^level
    const int NrL   = nr0 * scale;
    const int NzL   = nz0 * scale;

    // Previous level dimensions must match leftmost_civ_prev_level.
    const int prev_scale = 1 << (level - 1);
    const int NzPrev     = nz0 * prev_scale;

    if (static_cast<int>(leftmost_civ_prev_level.size()) != NzPrev)
    {
        throw std::runtime_error("ExtractCivSeeds_BoundaryBand: leftmost_civ_prev_level wrong size");
    }

    const auto &tp = cfg.tracing_params;
    const double r0 = tp.r_min, r1 = tp.r_max;
    const double z0 = tp.z_min, z1 = tp.z_max;

    if (!(r1 > r0) || !(z1 > z0))
    {
        throw std::runtime_error("ExtractCivSeeds_BoundaryBand: invalid bounds (r/z)");
    }

    const double dr = (r1 - r0) / static_cast<double>(NrL);
    const double dz = (z1 - z0) / static_cast<double>(NzL);

    const bool axisymmetric = cfg.solver.axisymmetric;

    // Collect unique (ir,iz) in this level's band.
    // Key encoding: high 32 bits = iz, low 32 bits = ir.
    std::vector<std::uint64_t> keys;
    keys.reserve(static_cast<std::size_t>(NzL) * 18); // rough: 2 r-bases * 3x3 halo

    bool any_parent_has_civ = false;

    for (int iz = 0; iz < NzL; ++iz)
    {
        const int iz_parent = iz / 2; // parent row at previous level
        const int b_prev    = leftmost_civ_prev_level[static_cast<std::size_t>(iz_parent)];

        if (b_prev < 0)
        {
            continue; // no CIV known in this parent row -> no band for its children
        }

        any_parent_has_civ = true;

        // Map previous boundary index to this level:
        // The parent cell splits into two in r, so boundary candidates are {2*b, 2*b+1}.
        const int ir_base0 = 2 * b_prev;
        const int ir_base1 = 2 * b_prev + 1;

        const int ir_bases[2] = { ir_base0, ir_base1 };

        for (int ib = 0; ib < 2; ++ib)
        {
            const int ir_base = ir_bases[ib];

            // Halo: ±1 in both r and z at current level.
            for (int diz = -1; diz <= 1; ++diz)
            {
                const int iz2 = iz + diz;
                if (iz2 < 0 || iz2 >= NzL) { continue; }

                for (int dir = -1; dir <= 1; ++dir)
                {
                    const int ir2 = ir_base + dir;
                    if (ir2 < 0 || ir2 >= NrL) { continue; }

                    const std::uint64_t key =
                        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(iz2)) << 32) |
                        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(ir2)));

                    keys.push_back(key);
                }
            }
        }
    }

    if (!any_parent_has_civ || keys.empty())
    {
        return Seeds{}; // empty => nothing to refine
    }

    // Unique keys to avoid duplicate seeds.
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    Seeds seeds;
    seeds.positions.resize(keys.size());
    seeds.volumes.resize(keys.size());

    for (std::size_t k = 0; k < keys.size(); ++k)
    {
        const std::uint32_t iz = static_cast<std::uint32_t>(keys[k] >> 32);
        const std::uint32_t ir = static_cast<std::uint32_t>(keys[k] & 0xFFFFFFFFu);

        const double r = r0 + (static_cast<double>(ir) + 0.5) * dr;
        const double z = z0 + (static_cast<double>(iz) + 0.5) * dz;

        Vector p(sdim);
        p = 0.0;
        p[0] = r;
        p[1] = z;

        seeds.positions[k] = std::move(p);

        double dV = dr * dz;
        if (axisymmetric)
        {
            dV *= (2.0 * M_PI * r);
        }
        seeds.volumes[k] = dV;
    }

    return seeds;
}


void UpdateLeftmostCivFromBand(const Config &cfg,
                               const SimulationResult & /*result*/,
                               int level,
                               const Seeds &band_seeds,
                               const std::vector<ElectronTraceResult> &band_results,
                               std::vector<int> &leftmost_civ_this_level)
{
    if (level <= 0)
    {
        throw std::runtime_error("UpdateLeftmostCivFromBand: level must be >= 1");
    }

    const int nr0 = cfg.civ_params.nr;
    const int nz0 = cfg.civ_params.nz;
    if (nr0 <= 0 || nz0 <= 0)
    {
        throw std::runtime_error("UpdateLeftmostCivFromBand: cfg.civ_params.nr/nz must be > 0");
    }

    const int scale = 1 << level;
    const int NrL   = nr0 * scale;
    const int NzL   = nz0 * scale;

    if (band_seeds.positions.size() != band_results.size())
    {
        throw std::runtime_error("UpdateLeftmostCivFromBand: seeds/results size mismatch");
    }

    leftmost_civ_this_level.assign(static_cast<std::size_t>(NzL), -1);

    const auto &tp = cfg.tracing_params;
    const double r0 = tp.r_min, r1 = tp.r_max;
    const double z0 = tp.z_min, z1 = tp.z_max;

    const double dr = (r1 - r0) / static_cast<double>(NrL);
    const double dz = (z1 - z0) / static_cast<double>(NzL);

    // Map each seed position back to its (ir,iz) at this level.
    // Since we generated exact centers, this should be stable.
    for (std::size_t k = 0; k < band_results.size(); ++k)
    {
        if (!IsChargeInsensitive(band_results[k]))
        {
            continue;
        }

        const auto &p = band_seeds.positions[k];
        if (p.Size() < 2)
        {
            throw std::runtime_error("UpdateLeftmostCivFromBand: seed position has Size() < 2");
        }

        const double r = p[0];
        const double z = p[1];

        // Compute indices by inverting center formula.
        // ir = floor((r - r0)/dr), iz = floor((z - z0)/dz)
        int ir = static_cast<int>(std::floor((r - r0) / dr));
        int iz = static_cast<int>(std::floor((z - z0) / dz));

        // Clamp indices defensively to grid bounds (not tolerance clamping).
        ir = std::min(std::max(ir, 0), NrL - 1);
        iz = std::min(std::max(iz, 0), NzL - 1);

        int &lm = leftmost_civ_this_level[static_cast<std::size_t>(iz)];
        if (lm < 0 || ir < lm)
        {
            lm = ir;
        }
    }
}



double IntegrateCIVFromBoundary(const Config &cfg,
                                const SimulationResult & /*result*/,
                                int level,
                                const std::vector<int> &leftmost_civ_per_row)
{
    const int nr0 = cfg.civ_params.nr;
    const int nz0 = cfg.civ_params.nz;

    if (nr0 <= 0 || nz0 <= 0)
    {
        throw std::runtime_error("IntegrateCIVFromBoundary: cfg.civ_params.nr/nz must be > 0");
    }
    if (level < 0)
    {
        throw std::runtime_error("IntegrateCIVFromBoundary: level must be >= 0");
    }

    const int scale = 1 << level;
    const int NrL   = nr0 * scale;
    const int NzL   = nz0 * scale;

    if (static_cast<int>(leftmost_civ_per_row.size()) != NzL)
    {
        throw std::runtime_error("IntegrateCIVFromBoundary: leftmost_civ_per_row has wrong size for level");
    }

    const auto &tp = cfg.tracing_params;
    const double r0 = tp.r_min, r1 = tp.r_max;
    const double z0 = tp.z_min, z1 = tp.z_max;

    if (!(r1 > r0) || !(z1 > z0))
    {
        throw std::runtime_error("IntegrateCIVFromBoundary: invalid bounds (r/z)");
    }

    const double dr = (r1 - r0) / static_cast<double>(NrL);
    const double dz = (z1 - z0) / static_cast<double>(NzL);

    const bool axisymmetric = cfg.solver.axisymmetric;

    // Total volume/area of ROI under chosen weighting.
    double V_total = 0.0;
    if (!axisymmetric)
    {
        V_total = (r1 - r0) * (z1 - z0);
    }
    else
    {
        // ∫∫ 2π r dr dz = π (r1^2 - r0^2) (z1 - z0)
        V_total = M_PI * (r1 * r1 - r0 * r0) * (z1 - z0);
    }

    double V_civ = 0.0;

    for (int iz = 0; iz < NzL; ++iz)
    {
        int ir_b = leftmost_civ_per_row[static_cast<std::size_t>(iz)];

        if (ir_b < 0)
        {
            continue; // no CIV in this row
        }

        // Defensive bounds (not tolerance clamping).
        ir_b = std::min(std::max(ir_b, 0), NrL);

        // Boundary at the left edge of leftmost CIV cell.
        const double r_boundary = r0 + static_cast<double>(ir_b) * dr;

        if (!axisymmetric)
        {
            V_civ += (r1 - r_boundary) * dz;
        }
        else
        {
            // ∫_{r_boundary}^{r1} 2π r dr * dz = π (r1^2 - r_boundary^2) dz
            V_civ += M_PI * (r1 * r1 - r_boundary * r_boundary) * dz;
        }
    }

    return (V_total > 0.0) ? (V_civ / V_total) : 0.0;
}


// -------------------- Random Sample -----------------------
double ComputeCIV_RandomSample(const Config            &cfg,
                               const SimulationResult &result)
{
    using namespace mfem;

    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (RandomSample)\n";
    }

    Seeds seeds = ExtractCivSeeds(cfg, result);

    // Construct Tracing Object
    ElectronFieldLineTracer tracer;
    tracer.Setup(result, cfg.tracing_params, cfg, cfg.solver.axisymmetric);

    if (seeds.positions.empty())
    {
        throw std::runtime_error("No seeds to compute CIV for");
    }

    std::vector<ElectronTraceResult> trace_results;
    tracer.Trace(seeds, trace_results,
                 cfg.solver.axisymmetric, 
                 cfg.debug.dumpdata,
                 nullptr);

    double V_total = 0.0;
    double V_civ   = 0.0;

    for (std::size_t i = 0; i < seeds.positions.size(); ++i)
    {
        const double               dV  = seeds.volumes[i];
        const ElectronTraceResult &res = trace_results[i];

        V_total += dV;

        if (IsChargeInsensitive(res))
        {
            V_civ += dV;
        }
    }
    if (cfg.debug.debug) { std::cout << "Traced " << seeds.positions.size() << " seed points" << std::endl;}
    return (V_total > 0.0) ? (V_civ / V_total) : 0.0;
}
// --------------------------  Fixed Grid ---------------------------
double ComputeCIV_FixedGrid(const Config            &cfg,
                            const SimulationResult &result)
{
    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (FixedGrid) "
                  << "nr=" << cfg.civ_params.nr << " nz=" << cfg.civ_params.nz << "\n";
    }

    Seeds seeds = ExtractCivFixedGridSeeds(cfg, result);

    // Construct Tracing Object
    ElectronFieldLineTracer tracer;
    tracer.Setup(result, cfg.tracing_params, cfg, cfg.solver.axisymmetric);


    if (seeds.positions.empty())
    {
        throw std::runtime_error("ComputeCIV_FixedGrid: no seeds generated");
    }
    if (seeds.volumes.size() != seeds.positions.size())
    {
        throw std::runtime_error("ComputeCIV_FixedGrid: volumes/positions size mismatch");
    }

    std::vector<ElectronTraceResult> trace_results;
    tracer.Trace(seeds, trace_results,
                 cfg.solver.axisymmetric, 
                 cfg.debug.dumpdata,
                 nullptr);

    if (trace_results.size() != seeds.positions.size())
    {
        throw std::runtime_error("ComputeCIV_FixedGrid: trace_results size mismatch");
    }

    double V_total = 0.0;
    double V_civ   = 0.0;

    for (std::size_t i = 0; i < seeds.positions.size(); ++i)
    {
        const double dV = seeds.volumes[i];
        V_total += dV;

        if (IsChargeInsensitive(trace_results[i]))
        {
            V_civ += dV;
        }
    }

    if (cfg.debug.debug) { std::cout << "Traced " << seeds.positions.size() << " seed points" << std::endl;}
    return (V_total > 0.0) ? (V_civ / V_total) : 0.0;
}

// ----------------------- Compute CIV Adaptive Grid -------------------

double ComputeCIV_AdaptiveGrid(const Config            &cfg,
                               const SimulationResult &result)
{
    if (cfg.civ_params.nr <= 0 || cfg.civ_params.nz <= 0)
    {
        throw std::runtime_error("ComputeCIV_AdaptiveGrid: cfg.civ_params.nr/nz must be > 0");
    }
    int seed_count = 0;

    const int nr0 = cfg.civ_params.nr;
    const int nz0 = cfg.civ_params.nz;

    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (AdaptiveGrid)\n"
                  << "  coarse nr=" << nr0 << " nz=" << nz0 << "\n";
    }

    // Construct Tracing Object
    ElectronFieldLineTracer tracer;
    tracer.Setup(result, cfg.tracing_params, cfg, cfg.solver.axisymmetric);

    // -------------------------------------------------------------------------
    // Level 0: right-to-left coarse column sweep
    // -------------------------------------------------------------------------
    std::vector<int> leftmost_civ_lvl0(static_cast<std::size_t>(nz0), -1);

    bool prev_col_had_any_civ = true; // allow at least one column to be processed
    int  last_traced_col      = nr0 - 1;

    for (int ir = nr0 - 1; ir >= 0; --ir)
    {
        // Termination rule: stop once the *previous* column had no CIV.
        if (!prev_col_had_any_civ)
        {
            break;
        }

        Seeds col_seeds = ExtractCivSeeds_Column(cfg, result, ir);
        seed_count += col_seeds.positions.size();
        if (col_seeds.positions.empty())
        {
            throw std::runtime_error("ComputeCIV_AdaptiveGrid: ExtractCivSeeds_Column produced no seeds");
        }

        std::vector<ElectronTraceResult> col_results;
        tracer.Trace(col_seeds, col_results,
                 cfg.solver.axisymmetric, 
                 cfg.debug.dumpdata,
                 nullptr);

        if (col_results.size() != col_seeds.positions.size())
        {
            throw std::runtime_error("ComputeCIV_AdaptiveGrid: column trace size mismatch");
        }

        bool col_has_any_civ = false;
        UpdateLeftmostCivFromColumn(cfg, result, ir,
                                    col_seeds, col_results,
                                    leftmost_civ_lvl0,
                                    col_has_any_civ);

        prev_col_had_any_civ = col_has_any_civ;
        last_traced_col      = ir;
    }

    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:OPTIMIZATION] AdaptiveGrid: level0 sweep complete. "
                  << "last_traced_col=" << last_traced_col << "\n";
    }

    // -------------------------------------------------------------------------
    // Refinement levels: bisection around boundary band (one trace call per level)
    // -------------------------------------------------------------------------
    const int max_levels = cfg.civ_params.max_levels;
    const double tol     = cfg.tracing_params.geom_tol;       

    // If you prefer: max_levels can be derived from tol and coarse dr/dz.
    // For now we treat it as config-driven.

    std::vector<int> leftmost_prev = leftmost_civ_lvl0;
    int final_level = 0;

    for (int level = 1; level <= max_levels; ++level)
    {
        // Stop condition by tolerance (conceptual here; exact check will be implemented
        // once we have consistent dr/dz computation for each level).
        // We keep the skeleton: "if reached tolerance -> break".
        // The actual criterion will live in a helper later.
        (void)tol;

        Seeds band_seeds = ExtractCivSeeds_BoundaryBand(cfg, result, level, leftmost_prev);

        seed_count += band_seeds.positions.size();
        if (band_seeds.positions.empty())
        {
            // No interface to refine; accept previous boundary as final.
            final_level = level - 1;
            break;
        }

        std::vector<ElectronTraceResult> band_results;
        tracer.Trace(band_seeds, band_results,
                 cfg.solver.axisymmetric, 
                 cfg.debug.dumpdata,
                 nullptr);

        if (band_results.size() != band_seeds.positions.size())
        {
            throw std::runtime_error("ComputeCIV_AdaptiveGrid: band trace size mismatch");
        }

        // Build next-level boundary representation. Size should correspond to nz(level).
        std::vector<int> leftmost_this; // will be sized by updater
        UpdateLeftmostCivFromBand(cfg, result, level,
                                  band_seeds, band_results,
                                  leftmost_this);

        // If updater decides boundary didn't change or nothing to refine further,
        // it can return identical boundary and/or mark a condition. For now, we just proceed.
        leftmost_prev = std::move(leftmost_this);
        final_level = level;

        if (cfg.debug.debug)
        {
            std::cout << "[DEBUG:OPTIMIZATION] AdaptiveGrid: completed refinement level "
                      << level << "\n";
        }
    }

    // If refinement loop never ran, boundary is level0.
    const std::vector<int> &leftmost_final = (final_level == 0) ? leftmost_civ_lvl0 : leftmost_prev;

    // -------------------------------------------------------------------------
    // Final integration using boundary representation (Option A)
    // -------------------------------------------------------------------------
    const double civ_fraction = IntegrateCIVFromBoundary(cfg, result, final_level, leftmost_final);

    if (cfg.debug.debug) { std::cout << "Traced " << seed_count << " seed points" << std::endl;}
    return civ_fraction;
}


// Row by Row Sweep fixed grid 
double ComputeCIV_RowSweep(const Config            &cfg,
                           const SimulationResult &result)
{
    const int nr = cfg.civ_params.nr;
    const int nz = cfg.civ_params.nz;

    int seed_count = 0;

    // Construct Tracing Object
    ElectronFieldLineTracer tracer;
    tracer.Setup(result, cfg.tracing_params, cfg, cfg.solver.axisymmetric);


    if (nr <= 0 || nz <= 0)
    {
        throw std::runtime_error("ComputeCIV_RowSweep: cfg.civ_params.nr/nz must be > 0");
    }

    if (cfg.debug.debug)
    {
        std::cout << "[DEBUG:OPTIMIZATION] Computing CIV (RowSweep) "
                  << "nr=" << nr << " nz=" << nz << "\n";
    }

    std::vector<int> leftmost_civ_per_row(static_cast<std::size_t>(nz), -1);

    // Optional monotonic carry from previous (lower) row:
    // boundary cannot move right when going up? (not specified)
    // We do not enforce cross-row monotonicity beyond the stated "below is CIV".
    // However, we can use the lower-row boundary as a starting point to reduce work:
    // if below-row boundary is at ir_below, then in the current row CIV (if exists)
    // cannot start to the right of ir_below (i.e., boundary index <= ir_below) is NOT guaranteed.
    // So we do NOT apply this unless you confirm. For now, always start at right edge.
    (void)leftmost_civ_per_row;

    // Bottom -> top
    for (int iz = 0; iz < nz; ++iz)
    {
        int boundary = -1;

        // Scan right -> left in blocks to keep tracing batches sizable.
        // Block size can be tuned; keep it simple and deterministic.
        const int block = std::max(1, cfg.civ_params.block_size); // assumed to exist; otherwise set e.g. 64

        int ir_right = nr - 1;
        bool seen_any_civ = false;
        bool done = false;

        while (!done && ir_right >= 0)
        {
            const int ir_left = std::max(0, ir_right - (block - 1));

            Seeds row_seeds = ExtractCivSeeds_Row(cfg, result, iz, ir_right, ir_left);
            seed_count += row_seeds.positions.size();
            std::vector<ElectronTraceResult> row_results;
            tracer.Trace(row_seeds, row_results,
                 cfg.solver.axisymmetric, 
                 cfg.debug.dumpdata,
                 nullptr);

            if (row_results.size() != row_seeds.positions.size())
            {
                throw std::runtime_error("ComputeCIV_RowSweep: row trace size mismatch");
            }

            // row_seeds are ordered right->left by construction in this call.
            for (std::size_t k = 0; k < row_results.size(); ++k)
            {
                const int ir = ir_right - static_cast<int>(k);

                if (IsChargeInsensitive(row_results[k]))
                {
                    seen_any_civ = true;
                    boundary = ir; // keep updating; we want leftmost CIV
                }
                else
                {
                    // If we've already seen CIV in this row and now we see non-CIV while moving left,
                    // then the boundary is between this cell and the previous one. We can stop.
                    if (seen_any_civ)
                    {
                        done = true;
                        break;
                    }
                }
            }

            ir_right = ir_left - 1;
        }

        leftmost_civ_per_row[static_cast<std::size_t>(iz)] = boundary;
    }

    // Integrate using the boundary representation at level 0.
    if (cfg.debug.debug) { std::cout << "Traced " << seed_count << " seed points" << std::endl;}
    return IntegrateCIVFromBoundary(cfg, result, /*level=*/0, leftmost_civ_per_row);
}