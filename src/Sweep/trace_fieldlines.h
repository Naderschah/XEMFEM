/*


FIXME  Lots of dead code here

*/
#pragma once

#include <mfem.hpp>
#include <vector>
#include "solver_api.h"
#include "Config.h"
#include "./common_tracing/tracing_objects.h"
#include "trace_fieldlines_MPI.h"

#include <memory>
#include <optional>
#include <functional>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkCellData.h>
#include <vtkStreamTracer.h>
#include <vtkAbstractInterpolatedVelocityField.h>

#include <vtkStaticCellLocator.h>
#include <vtkStreamTracer.h>
#include <vtkRungeKutta4.h>
#include <vtkRungeKutta2.h>

#include <vtkExtractGeometry.h>
#include <vtkBox.h>

#include "mfem/mesh/vtk.hpp"

// Trace  A single electron through the TPC
ElectronTraceResult TraceSingleElectronLine(
    mfem::ParMesh                    &mesh,
    mfem::FindPointsGSLIB            &finder,
    const mfem::ParGridFunction      &E_gf,
    const TpcGeometry                &geom,
    const ElectronTraceParams        &params,
    const mfem::Vector               &x0_in,          // (r,z,*) or (x,y,z)
    bool                              axisymmetric,
    bool                              save_pathlines,
    int                               max_traversals,
    mfem::Vector                     &pos_scratch,
    mfem::Vector                     &E_scratch,
    const double                      h_ref);

// struct dispatcher

#include <mfem.hpp>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <limits>

#include "trace_fieldlines.h" // for CivSeeds, ElectronTraceResult, Config, ElectronTraceParams, etc.

// Forward declarations for VTK types (to avoid pulling VTK headers here).
class vtkUnstructuredGrid;
class vtkStaticCellLocator;
template <typename T> class vtkSmartPointer;

// ----------------------------
// Cached BOOST context
// ----------------------------
struct BoostTraceContext
{
    std::unique_ptr<mfem::FindPointsGSLIB> finder;
    mfem::Vector pos;   // scratch
    mfem::Vector Eout;  // scratch
    double h_ref = 1.0;

    void Build(mfem::ParMesh &mesh, bool debug);
};

// ----------------------------
// Cached VTK context
// ----------------------------
struct VTKTraceContext
{
    vtkSmartPointer<vtkUnstructuredGrid> grid;
    vtkSmartPointer<vtkStaticCellLocator> locator;
    double h_ref = 1.0;

    void Build(mfem::ParMesh &pmesh,
               const mfem::ParGridFunction &E_gf,
               bool axisymmetric,
               bool debug);

    bool Ready() const;
};

// ----------------------------
// Unified tracer wrapper
// ----------------------------
struct ElectronFieldLineTracer
{
    using TraceFn = std::function<void(
        const Seeds &seeds,
        std::vector<ElectronTraceResult> &out_results,
        bool axisymmetric,
        bool save_paths,
        const double *z_max_overrides)>;

    // Non-owning mesh pointer (owned by SimulationResult)
    mfem::ParMesh *pmesh = nullptr;

    // Owned electric field objects (must persist after Setup)
    std::unique_ptr<mfem::H1_FECollection>       fec_vec;
    std::unique_ptr<mfem::ParFiniteElementSpace> fes_vec;
    std::unique_ptr<mfem::ParGridFunction>       E_gf_owned;

    // Non-owning view used everywhere else
    const mfem::ParGridFunction *E_gf = nullptr;


    // Params/config are stored (unchanged across calls):
    ElectronTraceParams params{};
    Config cfg{};

    // Cached contexts:
    std::optional<BoostTraceContext> boost;
    std::optional<VTKTraceContext>   vtk;
    std::optional<MPITraceContext>   mpitracer;

    // Hot-path callable:
    TraceFn trace_fn;

    // Lifecycle
    void Reset();

    // One-call setup: attach + build + select dispatcher
    void Setup(const SimulationResult &result,
               const ElectronTraceParams &tp,
               const Config &cfg_,
               bool axisymmetric);

    // Trace many times after Setup
    void Trace(const Seeds &seeds,
               std::vector<ElectronTraceResult> &out_results,
               bool axisymmetric,
               bool save_paths,
               const double *z_max_overrides = nullptr) const;

private:
    void BuildBOOST_(bool debug);
    void BuildVTK_(bool axisymmetric, bool debug);
    void BuildMPITracer_(bool debug);
    void SelectProvider_(const std::string &provider, bool axisymmetric);
};
