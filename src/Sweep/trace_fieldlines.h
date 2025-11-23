#pragma once

#include <mfem.hpp>
#include <memory>

#include "solver_api.h"

// Tracer interface types
enum class ElectronExitCode
{
    None = 0,
    HitCathode,
    HitLiquidGas,
    HitAxis,   
    HitWall,   
    LeftVolume,
    WeakField,
    MaxSteps,
    InvalidSeed
};

struct ElectronTraceParams
{
    double c_step         = 1;        // in Δt ≤ c_step * h_K / ||v||
    double max_time       = 1e3;      // max integration "time"
    int    max_steps      = 200;      // safety cap
    double min_field_norm = 0.01;      // stop if ||E|| < threshold
    double geom_tol       = 1e-6;     // tolerance for r,z boundary checks
};

struct ElectronTraceResult
{
    std::vector<mfem::Vector> points; // (r,z) 
    ElectronExitCode          exit_code    = ElectronExitCode::None;
    int                       exit_element = -1;
};

// ------------------------ Internal Helpers ------------------------
// For this we make one of each cell to speed up the lookup of next mesh elem
struct ElementAdjacency
{
    // list of face-neighbor
    std::vector<std::vector<int>> neighbors;
    // characteristic size of element
    std::vector<double>           h;
    // active[e]: whether element e lies inside the TPC volume of interest
    std::vector<bool>             active;
};
// Mesh element seed container for things we care about in CIV 
struct CivSeeds
{
    std::vector<mfem::Vector>        positions; // (r,z) seeds
    std::vector<double>              volumes;   // dV per seed
    std::vector<int>                 elements;  // element index for each seed
    std::vector<mfem::IntegrationPoint> ips;    // reference IP for each seed
};



// Electric field coefficient wrapper - just to have less to pass around
class ElectricFieldCoeff : public mfem::VectorCoefficient
{
public:
    explicit ElectricFieldCoeff(const mfem::GridFunction &E);

    void Eval(mfem::Vector &V,
              mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override;

private:
    const mfem::GridFunction &E_;
};

// Single electron tracing 
static ElectronTraceResult TraceSingleElectronLine(const mfem::ParMesh &mesh,const ElectricFieldCoeff &E_coeff,const ElementAdjacency &adj,const Config &cfg,int start_elem,const mfem::IntegrationPoint &start_ip,const ElectronTraceParams &params);

// Wrapper to trace all supplied electric field lines from starting points
void TraceElectronFieldLines(const SimulationResult &sim,const Config &cfg,const std::vector<mfem::Vector> &seed_points,const std::vector<int> &seed_elements,const std::vector<mfem::IntegrationPoint> &seed_ips,const ElectronTraceParams &params,std::vector<ElectronTraceResult> &out_results);

// Extract all the seeds for CIV
CivSeeds ExtractCivSeeds(const Config &cfg, const SimulationResult &result, int ir_order = 1);