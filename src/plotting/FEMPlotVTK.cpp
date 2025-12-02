// FEMPlotVTK.cpp
//
// Implementation of the FEMPlot API plus an example main().

#include "FEMPlotVTK.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>  
#include <cstring>

// VTK includes
#include <vtkPlane.h>
#include <vtkPlanes.h>
#include <vtkTextProperty.h>
#include <vtkPlanes.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCubeAxesActor.h>
#include <vtkCamera.h>
#include <vtkActorCollection.h>
#include <vtkActor2DCollection.h>
#include <vtkScalarBarActor.h>
#include <vtkProperty.h>
#include <vtkDoubleArray.h>
#include <vtkPointSource.h>
#include <vtkRungeKutta4.h>
#include <vtkStreamTracer.h>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkGradientFilter.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkScalarBarActor.h>
#include <vtkCubeAxesActor.h>
#include <vtkContourFilter.h>
#include <vtkColorSeries.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkGraphicsFactory.h>
#include <vtkActor.h>
#include <vtkGraphicsFactory.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkClipDataSet.h>
#include <vtkBox.h>
#include <vtkExtractGeometry.h>
// Mesh IO
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
// Scalars
#include <vtkPointData.h>
#include <vtkDataArray.h>
// Seeding for Streamlines
#include <vtkBoundedPointSource.h>

namespace FEMPlot
{


// -----------------------------------------------------------------------------
// Core IO
// -----------------------------------------------------------------------------

vtkSmartPointer<vtkUnstructuredGrid>
LoadGridFromVTU(const std::string& filename)
{
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid* ug = reader->GetOutput();
    if (!ug)
    {
        std::cerr << "[LoadGridFromVTU] ERROR: reader->GetOutput() is null for '"
                  << filename << "'.\n";
        return nullptr;
    }

    // Wrap in smart pointer explicitly (reader owns it internally, but this
    // pattern is what we already use elsewhere).
    auto result = vtkSmartPointer<vtkUnstructuredGrid>::New();
    result->ShallowCopy(ug);
    return result;
}

// -----------------------------------------------------------------------------
// Geometry clipping (ROI support)
// -----------------------------------------------------------------------------

static vtkSmartPointer<vtkUnstructuredGrid>
ClipToBox(vtkUnstructuredGrid* grid, const double bounds[6])
{
    if (!grid)
    {
        std::cerr << "[ClipToBox] ERROR: grid is null.\n";
        return nullptr;
    }

    double b[6] = {
        bounds[0], bounds[1],
        bounds[2], bounds[3],
        bounds[4], bounds[5]
    };

    // If this is effectively a 2D slice in z, give it a tiny thickness so
    // VTK filters don't get confused by zero extent.
    if (b[4] == b[5])
    {
        const double dz = 1e-6 * std::max({b[1] - b[0], b[3] - b[2], 1.0});
        b[4] -= dz * 0.5;
        b[5] += dz * 0.5;
    }

    vtkNew<vtkBox> box;
    box->SetBounds(b); // xmin,xmax, ymin,ymax, zmin,zmax

    vtkNew<vtkExtractGeometry> extract;
    extract->SetInputData(grid);
    extract->SetImplicitFunction(box);
    extract->ExtractInsideOn();         // keep cells inside the box
    extract->ExtractBoundaryCellsOn();
    extract->Update();

    auto clipped = vtkSmartPointer<vtkUnstructuredGrid>::New();
    clipped->ShallowCopy(extract->GetOutput());
    return clipped;
}

// -----------------------------------------------------------------------------
// Streamlines overlay
// -----------------------------------------------------------------------------

void AddStreamlines(
    vtkDataSet* dataset,               // MUST contain the geometry
    vtkRenderer* renderer,
    const Result& result,              // E-field at VTK points
    const StreamlineConfig& cfg)
{
    if (!dataset || !renderer)
    {
        std::cerr << "[AddStreamlines] ERROR: dataset or renderer is null.\n";
        return;
    }

    vtkIdType npts = dataset->GetNumberOfPoints();
    if (npts <= 0)
    {
        std::cerr << "[AddStreamlines] WARNING: dataset has no points.\n";
        return;
    }

    if (result.E_at_points.size() != static_cast<std::size_t>(npts))
    {
        std::cerr << "[AddStreamlines] ERROR: E_at_points size ("
                  << result.E_at_points.size()
                  << ") does not match number of VTK points ("
                  << npts << ").\n";
        return;
    }

    // ---------------------------------------------------------------------
    // 1) Attach / update vector field E(x) on the dataset
    // ---------------------------------------------------------------------
    auto pd = dataset->GetPointData();
    if (!pd)
    {
        std::cerr << "[AddStreamlines] ERROR: dataset has no point data.\n";
        return;
    }

    vtkDataArray* existingE = pd->GetArray("E");
    vtkSmartPointer<vtkDoubleArray> evec;

    if (existingE)
    {
        // Reuse existing array if it is 3-component and same size.
        if (existingE->GetNumberOfComponents() == 3 &&
            existingE->GetNumberOfTuples() == npts)
        {
            evec = vtkDoubleArray::SafeDownCast(existingE);
        }
        else
        {
            // Replace with a fresh, correctly-sized array.
            evec = vtkSmartPointer<vtkDoubleArray>::New();
            evec->SetName("E");
            evec->SetNumberOfComponents(3);
            evec->SetNumberOfTuples(npts);
            pd->RemoveArray("E");
            pd->AddArray(evec);
        }
    }
    else
    {
        evec = vtkSmartPointer<vtkDoubleArray>::New();
        evec->SetName("E");
        evec->SetNumberOfComponents(3);
        evec->SetNumberOfTuples(npts);
        pd->AddArray(evec);
    }

    for (vtkIdType i = 0; i < npts; ++i)
    {
        const auto& e = result.E_at_points[static_cast<std::size_t>(i)];
        double tuple[3] = { e[0], e[1], e[2] };
        evec->SetTuple(i, tuple);
    }

    pd->SetActiveVectors("E");

    // ---------------------------------------------------------------------
    // 2) Compute bounding box → characteristic radius for tube radius
    // ---------------------------------------------------------------------
    double bounds[6];
    dataset->GetBounds(bounds);

    const double dx = bounds[1] - bounds[0];
    const double dy = bounds[3] - bounds[2];
    const double dz = bounds[5] - bounds[4];
    const double radius = 0.5 * std::max({dx, dy, dz, 1e-6});

    // ---------------------------------------------------------------------
    // 3) Seed generator (within dataset bounds or, after clipping, within ROI)
    // ---------------------------------------------------------------------
    auto seeds = vtkSmartPointer<vtkBoundedPointSource>::New();
    seeds->SetNumberOfPoints(cfg.n_seeds);
    seeds->SetBounds(bounds);         // current dataset bounding box

    // ---------------------------------------------------------------------
    // 4) Stream tracer
    // ---------------------------------------------------------------------
    auto rk4 = vtkSmartPointer<vtkRungeKutta4>::New();

    auto tracer = vtkSmartPointer<vtkStreamTracer>::New();
    tracer->SetInputData(dataset);            // dataset WITH vectors
    tracer->SetSourceConnection(seeds->GetOutputPort());
    tracer->SetIntegrator(rk4);
    tracer->SetIntegrationDirectionToBoth();
    tracer->SetMaximumPropagation(cfg.max_propagation);
    tracer->SetInitialIntegrationStep(cfg.initial_step);
    tracer->SetMinimumIntegrationStep(cfg.min_step);
    tracer->SetMaximumIntegrationStep(cfg.max_step);
    tracer->SetMaximumNumberOfSteps(cfg.max_steps);
    tracer->SetComputeVorticity(false);       // avoid extra computation

    // ---------------------------------------------------------------------
    // 5) Tubes
    // ---------------------------------------------------------------------
    auto tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputConnection(tracer->GetOutputPort());
    tubes->SetRadius(radius * cfg.tube_radius_rel);
    tubes->SetNumberOfSides(10);
    tubes->CappingOn();

    // ---------------------------------------------------------------------
    // 6) Mapper + Actor
    // ---------------------------------------------------------------------
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tubes->GetOutputPort());
    mapper->ScalarVisibilityOff();            // plain-colored tubes

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    actor->GetProperty()->SetOpacity(0.75);
    actor->GetProperty()->SetInterpolationToPhong();

    renderer->AddActor(actor);
}

// -----------------------------------------------------------------------------
// Offscreen rendering
// -----------------------------------------------------------------------------

// New, more general overload for a pre-configured render window.
// This will be used later for multi-viewport layouts (field + colorbar).
static void RenderOffscreenPNG(vtkRenderWindow* renderWindow,
                               const std::string& filename)
{
    if (!renderWindow)
    {
        std::cerr << "RenderOffscreenPNG(window) ERROR: renderWindow is null\n";
        return;
    }

    // Ensure offscreen mode and render
    renderWindow->SetOffScreenRendering(1);
    renderWindow->Render();

    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    windowToImageFilter->SetInput(renderWindow);
    //windowToImageFilter->SetInputBufferTypeToRGBA();
    windowToImageFilter->ReadFrontBufferOff();
    windowToImageFilter->Update();

    vtkImageData* image = windowToImageFilter->GetOutput();
    if (!image || !image->GetScalarPointer())
    {
        std::cerr << "RenderOffscreenPNG(window) ERROR: "
                     "WindowToImageFilter produced null image\n";
        return;
    }

    vtkNew<vtkPNGWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();
}

// Backwards-compatible wrapper for the existing single-renderer usage.
// Internally forwards to the window-based overload above.
void RenderOffscreenPNG(vtkRenderer* renderer,
                        const std::string& filename,
                        int width,
                        int height)
{
    if (!renderer)
    {
        std::cerr << "RenderOffscreenPNG(renderer) ERROR: renderer is null\n";
        return;
    }

    // If renderer is already attached to some window, detach it first.
    if (renderer->GetRenderWindow())
    {
        std::cerr << "RenderOffscreenPNG(renderer): renderer already "
                     "has a RenderWindow; detaching it before capture.\n";
        renderer->GetRenderWindow()->RemoveRenderer(renderer);
    }

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetOffScreenRendering(1);
    renderWindow->SetSize(width, height);
    renderWindow->AddRenderer(renderer);
    // Ensure the window background is white as well

    // Use the new, window-based helper
    RenderOffscreenPNG(renderWindow, filename);

    // Detach the renderer so it is in a clean state if you ever reuse it.
    renderWindow->RemoveRenderer(renderer);
}


vtkSmartPointer<vtkActor>
CreateMeshActor(vtkSmartPointer<vtkDataSetMapper> mapper,
                bool show_edges,
                double line_width)
{
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    auto prop = actor->GetProperty();
    prop->SetEdgeVisibility(show_edges ? 1 : 0);
    prop->SetLineWidth(line_width);

    // Optionally tweak surface appearance here if needed, e.g.:
    // prop->SetInterpolationToPhong();

    return actor;
}

void AddXYAxes(vtkRenderer* renderer,
               const double bounds[6],
               const std::string& x_label,
               const std::string& y_label,
               const std::string& z_label,
               int axis_title_font_size,
               int axis_label_font_size)
{
    if (!renderer)
        return;

    auto axes = vtkSmartPointer<vtkCubeAxesActor>::New();
    axes->SetBounds(bounds);
    axes->SetCamera(renderer->GetActiveCamera());

    axes->SetScreenSize(15.0);

    axes->SetXTitle(x_label.c_str());
    axes->SetYTitle(y_label.c_str());
    axes->SetZTitle(z_label.c_str());

    axes->XAxisLabelVisibilityOn();
    axes->YAxisLabelVisibilityOn();
    axes->ZAxisLabelVisibilityOn();

    axes->SetFlyModeToClosestTriad();

    // Symmetric tick marks (both sides)
    axes->SetTickLocation(vtkCubeAxesActor::VTK_TICKS_BOTH);

    // Text styling for all axes with per-view font sizes
    for (int i = 0; i < 3; ++i)
    {
        if (auto tprop = axes->GetTitleTextProperty(i))
        {
            tprop->SetFontFamilyToArial();
            tprop->SetColor(0.0, 0.0, 0.0);
            tprop->SetFontSize(axis_title_font_size);
        }
        if (auto lprop = axes->GetLabelTextProperty(i))
        {
            lprop->SetFontFamilyToArial();
            lprop->SetColor(0.0, 0.0, 0.0);
            lprop->SetFontSize(axis_label_font_size);
        }
    }

    axes->GetXAxesLinesProperty()->SetColor(0.0, 0.0, 0.0);
    axes->GetYAxesLinesProperty()->SetColor(0.0, 0.0, 0.0);
    axes->GetZAxesLinesProperty()->SetColor(0.0, 0.0, 0.0);

    renderer->AddActor(axes);
}

// Adjust image width/height so that the FIELD viewport has roughly the same
// aspect ratio as the physical region, while reserving some room for title
// and colorbar. The longest side is capped to 3840 pixels.
static PlottingOptions
FitPlottingOptionsToRegion(const PlottingOptions& base,
                           const double[6],
                           bool)
{
    return base;  // fixed resolution for all plots
}

// Compute scalar statistics on the given dataset for the given array name.
// Currently: uses all points of the dataset; the "ROI" is whatever dataset
// you pass (so if you cropped first, this is already local).
static ScalarStats ComputeScalarStats(vtkDataSet* ds,
                                      const std::string& scalar_name)
{
    ScalarStats stats;

    if (!ds)
    {
        std::cerr << "[ComputeScalarStats] ERROR: dataset is null.\n";
        return stats;
    }

    vtkPointData* pd = ds->GetPointData();
    if (!pd)
    {
        std::cerr << "[ComputeScalarStats] ERROR: no point data.\n";
        return stats;
    }

    vtkDataArray* arr = pd->GetArray(scalar_name.c_str());
    if (!arr)
    {
        std::cerr << "[ComputeScalarStats] ERROR: scalar array '"
                  << scalar_name << "' not found.\n";
        return stats;
    }

    vtkIdType n = arr->GetNumberOfTuples();
    if (n <= 0)
    {
        std::cerr << "[ComputeScalarStats] WARNING: scalar array '"
                  << scalar_name << "' has no tuples.\n";
        return stats;
    }

    std::vector<double> vals;
    vals.reserve(static_cast<std::size_t>(n));

    double minv = std::numeric_limits<double>::infinity();
    double maxv = -std::numeric_limits<double>::infinity();

    for (vtkIdType i = 0; i < n; ++i)
    {
        double v = arr->GetComponent(i, 0);
        vals.push_back(v);
        if (v < minv) minv = v;
        if (v > maxv) maxv = v;
    }

    stats.min_roi = minv;
    stats.max_roi = maxv;

    // Percentile-based clipping to avoid tiny outliers dominating the colorbar.
    // Use 2%–98% as a first reasonable default.
    std::sort(vals.begin(), vals.end());

    auto get_percentile = [&](double p) -> double {
        if (vals.empty()) return 0.0;
        if (p <= 0.0) return vals.front();
        if (p >= 100.0) return vals.back();
        double pos = p * (vals.size() - 1) / 100.0;
        auto idx = static_cast<std::size_t>(pos);
        if (idx >= vals.size() - 1) return vals.back();
        double alpha = pos - idx;
        return (1.0 - alpha) * vals[idx] + alpha * vals[idx + 1];
    };

    double p_low  = get_percentile(5.0);
    double p_high = get_percentile(95.0);

    // If the spread is too small, fall back to full min/max.
    const double eps = 1e-12;
    if (std::fabs(p_high - p_low) < eps)
    {
        stats.min_used = minv;
        stats.max_used = maxv;
    }
    else
    {
        stats.min_used = p_low;
        stats.max_used = p_high;
    }

    return stats;
}
static ScalarStats ComputeScalarStatsInBounds(vtkDataSet* ds,
                                              const std::string& scalar_name,
                                              const double bounds[6])
{
    ScalarStats stats;

    if (!ds)
    {
        std::cerr << "[ComputeScalarStatsInBounds] ERROR: dataset is null.\n";
        return stats;
    }

    vtkPointData* pd = ds->GetPointData();
    if (!pd)
    {
        std::cerr << "[ComputeScalarStatsInBounds] ERROR: no point data.\n";
        return stats;
    }

    vtkDataArray* arr = pd->GetArray(scalar_name.c_str());
    if (!arr)
    {
        std::cerr << "[ComputeScalarStatsInBounds] ERROR: scalar array '"
                  << scalar_name << "' not found.\n";
        return stats;
    }

    vtkIdType npts = ds->GetNumberOfPoints();
    if (npts <= 0)
    {
        std::cerr << "[ComputeScalarStatsInBounds] WARNING: no points.\n";
        return stats;
    }

    std::vector<double> vals;
    vals.reserve(static_cast<std::size_t>(npts));

    double minv = std::numeric_limits<double>::infinity();
    double maxv = -std::numeric_limits<double>::infinity();

    for (vtkIdType i = 0; i < npts; ++i)
    {
        double p[3];
        ds->GetPoint(i, p);

        if (p[0] < bounds[0] || p[0] > bounds[1] ||
            p[1] < bounds[2] || p[1] > bounds[3] ||
            p[2] < bounds[4] || p[2] > bounds[5])
        {
            continue;
        }

        double v = arr->GetComponent(i, 0);
        vals.push_back(v);
        if (v < minv) minv = v;
        if (v > maxv) maxv = v;
    }

    if (vals.empty())
    {
        std::cerr << "[ComputeScalarStatsInBounds] WARNING: no points in ROI, "
                     "falling back to global stats.\n";
        return ComputeScalarStats(ds, scalar_name);
    }

    stats.min_roi = minv;
    stats.max_roi = maxv;

    std::sort(vals.begin(), vals.end());

    auto get_percentile = [&](double p) -> double {
        if (vals.empty()) return 0.0;
        if (p <= 0.0) return vals.front();
        if (p >= 100.0) return vals.back();
        double pos = p * (vals.size() - 1) / 100.0;
        auto idx = static_cast<std::size_t>(pos);
        if (idx >= vals.size() - 1) return vals.back();
        double alpha = pos - idx;
        return (1.0 - alpha) * vals[idx] + alpha * vals[idx + 1];
    };

    double p_low  = get_percentile(5.0);
    double p_high = get_percentile(95.0);

    const double eps = 1e-12;
    if (std::fabs(p_high - p_low) < eps)
    {
        stats.min_used = minv;
        stats.max_used = maxv;
    }
    else
    {
        stats.min_used = p_low;
        stats.max_used = p_high;
    }

    return stats;
}

// Basic LUT builder. Right now we support a few simple "names";
// otherwise we fall back to a default blue-red like map.
static vtkSmartPointer<vtkLookupTable>
BuildLookupTable(const std::string& color_map_name,
                 double range_min,
                 double range_max)
{
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues(256);
    lut->SetRange(range_min, range_max);

    // Simple hard-coded styles for now.
    // You can refine this later with vtkColorSeries, etc.
    if (color_map_name == "gray" || color_map_name == "grey")
    {
        lut->SetHueRange(0.0, 0.0);      // no hue
        lut->SetSaturationRange(0.0, 0.0);
        lut->SetValueRange(0.0, 1.0);    // black -> white
    }
    else if (color_map_name == "blue-red" || color_map_name == "coolwarm")
    {
        lut->SetHueRange(0.6667, 0.0);   // blue -> red
        lut->SetSaturationRange(1.0, 1.0);
        lut->SetValueRange(1.0, 1.0);
    }
    else if (color_map_name == "viridis")
    {
        // Approximate "viridis"-like (greenish-blue to yellowish).
        lut->SetHueRange(0.7, 0.1);
        lut->SetSaturationRange(1.0, 1.0);
        lut->SetValueRange(0.3, 1.0);
    }
    else
    {
        // Default: VTK classic (blue->red)
        lut->SetHueRange(0.6667, 0.0);
        lut->SetSaturationRange(1.0, 1.0);
        lut->SetValueRange(1.0, 1.0);
    }

    lut->Build();
    return lut;
}

// Helper to create a scalar bar actor for a given lookup table and title.
static vtkSmartPointer<vtkScalarBarActor>
CreateScalarBar(vtkLookupTable* lut,
                const std::string& title,
                int title_font_size,
                int label_font_size,
                int num_labels)
{
    auto scalar_bar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalar_bar->SetLookupTable(lut);
    scalar_bar->SetTitle(title.c_str());
    scalar_bar->SetNumberOfLabels(num_labels > 0 ? num_labels : 5);

    auto tprop = scalar_bar->GetTitleTextProperty();
    tprop->SetFontFamilyToArial();
    tprop->SetFontSize(title_font_size);
    tprop->SetColor(0.0, 0.0, 0.0);

    auto lprop = scalar_bar->GetLabelTextProperty();
    lprop->SetFontFamilyToArial();
    lprop->SetFontSize(label_font_size);
    lprop->SetColor(0.0, 0.0, 0.0);

    // Default geometry (overridden when using separate cbar viewport)
    scalar_bar->SetWidth(0.1);
    scalar_bar->SetHeight(0.7);
    scalar_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalar_bar->GetPositionCoordinate()->SetValue(0.90, 0.15);

    return scalar_bar;
}
bool AttachGradientEAndMagnitude(vtkUnstructuredGrid* grid,
                                 const std::string& potential_array_name,
                                 const std::string& e_array_name,
                                 const std::string& emag_array_name)
{
    if (!grid)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: grid is null.\n";
        return false;
    }

    vtkPointData* pd = grid->GetPointData();
    if (!pd)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: no point data.\n";
        return false;
    }

    vtkDataArray* pot = pd->GetArray(potential_array_name.c_str());
    if (!pot)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: potential array '"
                  << potential_array_name << "' not found.\n";
        return false;
    }

    vtkIdType npts_grid = grid->GetNumberOfPoints();
    if (pot->GetNumberOfTuples() != npts_grid)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: potential array '"
                  << potential_array_name << "' size mismatch.\n";
        return false;
    }

    // ---------------------------------------------------------------------
    // 1) Compute gradient of potential with vtkGradientFilter
    // ---------------------------------------------------------------------
    auto grad = vtkSmartPointer<vtkGradientFilter>::New();
    grad->SetInputData(grid);
    grad->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS,
                          potential_array_name.c_str());
    grad->SetResultArrayName("GradTemp");
    grad->Update();

    vtkDataSet* ds_grad = grad->GetOutput();
    if (!ds_grad)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: "
                     "gradient filter returned null output.\n";
        return false;
    }

    vtkPointData* pd_grad = ds_grad->GetPointData();
    if (!pd_grad)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: "
                     "no point data on gradient output.\n";
        return false;
    }

    vtkDataArray* gradArray = pd_grad->GetArray("GradTemp");
    if (!gradArray)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: "
                     "GradTemp array not found.\n";
        return false;
    }

    int comps = gradArray->GetNumberOfComponents();
    if (comps < 1)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: GradTemp has "
                  << comps << " components.\n";
        return false;
    }

    vtkIdType npts_grad = ds_grad->GetNumberOfPoints();
    if (npts_grad != npts_grid)
    {
        std::cerr << "[AttachGradientEAndMagnitude] ERROR: point count mismatch, "
                  << "grid=" << npts_grid << " grad=" << npts_grad << "\n";
        return false;
    }

    // ---------------------------------------------------------------------
    // 2) Build E = -∇V as a vector array on the original grid
    // ---------------------------------------------------------------------
    auto evec = vtkSmartPointer<vtkDoubleArray>::New();
    evec->SetName(e_array_name.c_str());
    evec->SetNumberOfComponents(3);
    evec->SetNumberOfTuples(npts_grid);

    for (vtkIdType i = 0; i < npts_grid; ++i)
    {
        double g[3] = {0.0, 0.0, 0.0};
        int ncopy = std::min(comps, 3);
        for (int c = 0; c < ncopy; ++c)
            g[c] = gradArray->GetComponent(i, c);

        // E = -∇V
        double e[3] = { -g[0], -g[1], -g[2] };
        evec->SetTuple(i, e);
    }

    pd->AddArray(evec);
    pd->SetActiveVectors(e_array_name.c_str());

    // ---------------------------------------------------------------------
    // 3) Build |E| as scalar
    // ---------------------------------------------------------------------
    auto emag = vtkSmartPointer<vtkDoubleArray>::New();
    emag->SetName(emag_array_name.c_str());
    emag->SetNumberOfComponents(1);
    emag->SetNumberOfTuples(npts_grid);

    for (vtkIdType i = 0; i < npts_grid; ++i)
    {
        double e[3];
        evec->GetTuple(i, e);
        double mag = std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
        emag->SetValue(i, mag);
    }

    pd->AddArray(emag);

    return true;
}
void PlotScalarFieldView(vtkUnstructuredGrid* grid,
                         const ScalarViewRequest& request,
                         const PlottingOptions& plotting,
                         const StreamlineConfig* /*stream_cfg*/)
{
    if (!grid)
    {
        std::cerr << "[PlotScalarFieldView] ERROR: grid is null.\n";
        return;
    }

    vtkDataSet* dataset_to_plot = grid;

    vtkPointData* pd = dataset_to_plot->GetPointData();
    if (!pd)
    {
        std::cerr << "[PlotScalarFieldView] ERROR: dataset has no point data.\n";
        return;
    }

    vtkDataArray* scalars = pd->GetArray(request.scalar_name.c_str());
    if (!scalars)
    {
        std::cerr << "[PlotScalarFieldView] ERROR: scalar array '"
                  << request.scalar_name << "' not found on dataset.\n";
        return;
    }

    // ---------------------------------------------------------------------
    // 1) Scalar statistics and LUT
    // ---------------------------------------------------------------------
    ScalarStats stats;
    if (request.crop_to_region)
    {
        stats = ComputeScalarStatsInBounds(dataset_to_plot,
                                           request.scalar_name,
                                           request.region_bounds);
    }
    else
    {
        stats = ComputeScalarStats(dataset_to_plot, request.scalar_name);
    }

    auto lut = BuildLookupTable(request.color_map_name,
                                stats.min_used,
                                stats.max_used);

    // ---------------------------------------------------------------------
    // 2) Mapper and mesh actor (with optional clipping planes)
    // ---------------------------------------------------------------------
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(dataset_to_plot);
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(request.scalar_name.c_str());
    mapper->SetScalarRange(stats.min_used, stats.max_used);
    mapper->SetLookupTable(lut);
    mapper->ScalarVisibilityOn();

    if (request.crop_to_region)
    {
        double xmin = request.region_bounds[0];
        double xmax = request.region_bounds[1];
        double ymin = request.region_bounds[2];
        double ymax = request.region_bounds[3];

        mapper->RemoveAllClippingPlanes();

        vtkNew<vtkPlane> pxMin;
        pxMin->SetOrigin(xmin, 0.0, 0.0);
        pxMin->SetNormal(1.0, 0.0, 0.0);
        mapper->AddClippingPlane(pxMin);

        vtkNew<vtkPlane> pxMax;
        pxMax->SetOrigin(xmax, 0.0, 0.0);
        pxMax->SetNormal(-1.0, 0.0, 0.0);
        mapper->AddClippingPlane(pxMax);

        vtkNew<vtkPlane> pyMin;
        pyMin->SetOrigin(0.0, ymin, 0.0);
        pyMin->SetNormal(0.0, 1.0, 0.0);
        mapper->AddClippingPlane(pyMin);

        vtkNew<vtkPlane> pyMax;
        pyMax->SetOrigin(0.0, ymax, 0.0);
        pyMax->SetNormal(0.0, -1.0, 0.0);
        mapper->AddClippingPlane(pyMax);
    }
    else
    {
        mapper->RemoveAllClippingPlanes();
    }

    auto mesh_actor = CreateMeshActor(mapper,
                                      plotting.show_edges,
                                      plotting.edge_width);

    // ---------------------------------------------------------------------
    // 3) Render window + layout + renderers
    // ---------------------------------------------------------------------
    int img_w = (request.image_width  > 0) ? request.image_width  : plotting.image_width;
    int img_h = (request.image_height > 0) ? request.image_height : plotting.image_height;

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetOffScreenRendering(1);
    renderWindow->SetAlphaBitPlanes(0);
    renderWindow->SetMultiSamples(0);
    renderWindow->SetSize(img_w, img_h);

    // Layout for this plot
    PlotViewLayout layout = MakeDefaultPlotViewLayout(
        request.cbar_horizontal,
        request.separate_cbar_viewport);

    auto makeRenderer = [&](const ViewportRect& r) {
        auto ren = vtkSmartPointer<vtkRenderer>::New();
        ren->SetBackground(1.0, 1.0, 1.0);
        ren->SetBackground2(1.0, 1.0, 1.0);
        ren->GradientBackgroundOff();
        ren->SetViewport(r.x0, r.y0, r.x1, r.y1);
        renderWindow->AddRenderer(ren);
        return ren;
    };

    // Main field renderer
    vtkSmartPointer<vtkRenderer> field_renderer = makeRenderer(layout.field);
    field_renderer->AddActor(mesh_actor);

    // Title / cbar / label renderers
    vtkSmartPointer<vtkRenderer> title_renderer           = makeRenderer(layout.title);
    vtkSmartPointer<vtkRenderer> cbar_renderer            = makeRenderer(layout.cbar);
    vtkSmartPointer<vtkRenderer> cbar_title_renderer      = makeRenderer(layout.cbarTitle);
    vtkSmartPointer<vtkRenderer> cbar_label_top_renderer  = makeRenderer(layout.cbarLabelTopOrLeft);
    vtkSmartPointer<vtkRenderer> cbar_label_bot_renderer  = makeRenderer(layout.cbarLabelBottomOrRight);

    // ---------------------------------------------------------------------
    // 4) Axes (remain inside field viewport)
    // ---------------------------------------------------------------------
    AddXYAxes(field_renderer, request.region_bounds,
              request.x_label, request.y_label, request.z_label,
              request.axis_title_font_size,
              request.axis_label_font_size);

    // ---------------------------------------------------------------------
    // 5) Plot title in its own viewport
    // ---------------------------------------------------------------------
    if (!request.title.empty())
    {
        auto titleActor = vtkSmartPointer<vtkTextActor>::New();
        titleActor->SetInput(request.title.c_str());

        auto tprop = titleActor->GetTextProperty();
        tprop->SetFontFamilyToArial();
        tprop->SetFontSize(request.title_font_size);
        tprop->SetColor(0.0, 0.0, 0.0);
        tprop->SetJustificationToCentered();
        tprop->SetVerticalJustificationToCentered();

        titleActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        titleActor->SetPosition(0.5, 0.5);  // centered in title viewport

        title_renderer->AddActor2D(titleActor);
    }

    // ---------------------------------------------------------------------
    // 6) Scalar bar (colorbar viewport)
    // ---------------------------------------------------------------------
    auto scalar_bar = CreateScalarBar(
        lut,
        "",  // we show the title in a separate viewport below
        request.cbar_title_font_size,
        request.cbar_label_font_size,
        request.cbar_num_labels);

    if (request.cbar_horizontal)
    {
        scalar_bar->SetOrientationToHorizontal();
        scalar_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        scalar_bar->GetPositionCoordinate()->SetValue(0.10, 0.25);
        scalar_bar->SetWidth(0.80);
        scalar_bar->SetHeight(0.5);
    }
    else
    {
        scalar_bar->SetOrientationToVertical();
        scalar_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        scalar_bar->GetPositionCoordinate()->SetValue(0.35, 0.05);
        scalar_bar->SetWidth(0.30);
        scalar_bar->SetHeight(0.90);
    }

    if (request.show_contours)
    {
        // Placeholder: contour-aligned labels could be set here later.
        scalar_bar->UseCustomLabelsOff();
    }
    else
    {
        scalar_bar->UseCustomLabelsOff();
    }

    cbar_renderer->AddActor2D(scalar_bar);

    // ---------------------------------------------------------------------
    // 7) Colorbar title in its own viewport
    // ---------------------------------------------------------------------
    if (!request.cbar_title.empty())
    {
        auto cbarTitleActor = vtkSmartPointer<vtkTextActor>::New();
        cbarTitleActor->SetInput(request.cbar_title.c_str());

        auto tprop = cbarTitleActor->GetTextProperty();
        tprop->SetFontFamilyToArial();
        tprop->SetFontSize(request.cbar_title_font_size);
        tprop->SetColor(0.0, 0.0, 0.0);
        tprop->SetJustificationToCentered();
        tprop->SetVerticalJustificationToCentered();

        cbarTitleActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        cbarTitleActor->SetPosition(0.5, 0.5);  // centered in cbarTitle viewport

        cbar_title_renderer->AddActor2D(cbarTitleActor);
    }

    // ---------------------------------------------------------------------
    // 8) True min/max annotations in dedicated label viewports
    // ---------------------------------------------------------------------
    if (stats.min_used > stats.min_roi)
    {
        auto txt = vtkSmartPointer<vtkTextActor>::New();
        std::ostringstream ss;
        ss << "min = " << stats.min_roi;
        txt->SetInput(ss.str().c_str());

        auto tp = txt->GetTextProperty();
        tp->SetFontFamilyToArial();
        tp->SetFontSize(request.cbar_minmax_font_size);
        tp->SetColor(0.0, 0.0, 0.0);
        tp->SetJustificationToCentered();
        tp->SetVerticalJustificationToCentered();

        txt->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        txt->SetPosition(0.5, 0.5);

        cbar_label_bot_renderer->AddActor2D(txt);
    }

    if (stats.max_used < stats.max_roi)
    {
        auto txt = vtkSmartPointer<vtkTextActor>::New();
        std::ostringstream ss;
        ss << "max = " << stats.max_roi;
        txt->SetInput(ss.str().c_str());

        auto tp = txt->GetTextProperty();
        tp->SetFontFamilyToArial();
        tp->SetFontSize(request.cbar_minmax_font_size);
        tp->SetColor(0.0, 0.0, 0.0);
        tp->SetJustificationToCentered();
        tp->SetVerticalJustificationToCentered();

        txt->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        txt->SetPosition(0.5, 0.5);

        cbar_label_top_renderer->AddActor2D(txt);
    }

    // ---------------------------------------------------------------------
    // 9) Streamlines (not yet wired in this refactor)
    // ---------------------------------------------------------------------
    if (request.show_streamlines)
    {
        std::cerr << "[PlotScalarFieldView] NOTE: show_streamlines requested "
                     "but streamlines are not yet wired into the new view "
                     "pipeline.\n";
    }

    // ---------------------------------------------------------------------
    // 10) Camera setup and zoom for field renderer
    // ---------------------------------------------------------------------
    double b[6];
    std::copy(std::begin(request.region_bounds),
              std::end(request.region_bounds), b);

    if (b[4] == b[5])
    {
        const double dz = 1e-3;
        b[4] -= dz * 0.8;
        b[5] += dz * 0.8;
    }

    field_renderer->ResetCamera(b);
    if (auto cam = field_renderer->GetActiveCamera())
    {
        cam->SetParallelProjection(true);
        field_renderer->ResetCameraClippingRange();

        double zf = request.zoom_factor;
        if (zf <= 0.0) zf = 1.0;
        double ps = cam->GetParallelScale();
        cam->SetParallelScale(ps * zf);
    }

    // ---------------------------------------------------------------------
    // 11) Render and write PNG
    // ---------------------------------------------------------------------
    RenderOffscreenPNG(renderWindow, request.output_path);
}

void PlotStandardViewsForPotentialAndEnorm(vtkUnstructuredGrid* grid,
                                           const PlotInput& input,
                                           const PlottingOptions& base_options)
{
    if (!grid)
    {
        std::cerr << "[PlotStandardViewsForPotentialAndEnorm] ERROR: grid is null.\n";
        return;
    }

    double full_bounds[6];
    grid->GetBounds(full_bounds);

    if (!AttachGradientEAndMagnitude(grid, "V", "E", "Enorm"))
    {
        std::cerr << "[PlotStandardViewsForPotentialAndEnorm] WARNING: "
                     "failed to attach E/Enorm from V. |E| views may be skipped.\n";
    }

    const std::string xlab = "x [m]";
    const std::string ylab = "y [m]";
    const std::string zlab = "z [m]";

    double full[6];
    std::copy(std::begin(full_bounds), std::end(full_bounds), full);

    double ymin = full_bounds[2];
    double ymax = full_bounds[3];
    double xmin = full_bounds[0];
    double xmax = full_bounds[1];
    double zmin = full_bounds[4];
    double zmax = full_bounds[5];

    double top[6]       = { xmin, xmax, -0.13,  0.07,  zmin, zmax };
    double bottom[6]    = { xmin, xmax, -1.6,  -1.4,  zmin, zmax };
    double right_bar[6] = { xmax - 0.1, xmax, -1.45, -0.05, zmin, zmax };

    const int img_w = base_options.image_width;
    const int img_h = base_options.image_height;

    const int title_fs       = 32;
    const int axis_title_fs  = 28;
    const int axis_label_fs  = 22;
    const int cbar_title_fs  = 28;
    const int cbar_label_fs  = 22;
    const int cbar_minmax_fs = 20;
    const int cbar_labels    = 5;

    // ---------------- V views ----------------
        // --- Full-domain Potential V view ---------------------------------------
    ScalarViewRequest r_full;
    r_full.scalar_name           = "V";                     // dataset array to visualize
    r_full.cbar_title            = "V";                     // title text shown above the colorbar
    r_full.cbar_horizontal       = false;                   // colorbar orientation (vertical here)
    r_full.title                 = "Potential V";           // plot title shown at top-left
    r_full.x_label               = xlab;                    // axis label for X
    r_full.y_label               = ylab;                    // axis label for Y
    r_full.z_label               = zlab;                    // axis label for Z (unused in 2D)
    r_full.color_map_name        = "viridis";               // LUT preset name
    std::copy(full, full + 6, r_full.region_bounds);        // physical region plotted (full domain)
    r_full.crop_to_region        = false;                   // do NOT crop mesh; full geometry visible
    r_full.separate_cbar_viewport= true;                    // use dedicated viewport for colorbar
    r_full.show_contours         = true;                    // overlay contour lines
    r_full.show_streamlines      = false;                   // no streamlines for potential
    r_full.n_contours            = 10;                      // number of contour isolines
    r_full.output_path           = input.out_dir + "/V_full.png";  // output PNG path
    r_full.image_width           = 1920;                    // final rendered image width  (px)
    r_full.image_height          = 1080;                    // final rendered image height (px)
    r_full.zoom_factor           = 0.5;                     // orthographic camera zoom (<1 = zoom in)
    r_full.title_font_size       = 70;                      // title text size in pixels
    r_full.axis_title_font_size  = 700;                      // axis title font size (X,Y,Z titles)
    r_full.axis_label_font_size  = 700;                      // tick-label font size
    r_full.cbar_title_font_size  = 20;                      // colorbar title font size
    r_full.cbar_label_font_size  = 400;           // colorbar tick-label font size
    r_full.cbar_minmax_font_size = cbar_minmax_fs;          // annotation font size for full min/max
    r_full.cbar_num_labels       = cbar_labels;             // number of labels shown on the colorbar


    ScalarViewRequest r_top = r_full;
    r_top.cbar_horizontal       = true;
    r_top.title                 = "Potential V (top stack)";
    std::copy(top, top + 6, r_top.region_bounds);
    r_top.crop_to_region        = true;
    r_top.output_path           = input.out_dir + "/V_top.png";

    ScalarViewRequest r_bottom = r_full;
    r_bottom.cbar_horizontal    = true;
    r_bottom.title              = "Potential V (bottom stack)";
    std::copy(bottom, bottom + 6, r_bottom.region_bounds);
    r_bottom.crop_to_region     = true;
    r_bottom.output_path        = input.out_dir + "/V_bottom.png";

    ScalarViewRequest r_bar = r_full;
    r_bar.cbar_horizontal      = false;
    r_bar.title                = "Potential V (Field Cage)";
    std::copy(right_bar, right_bar + 6, r_bar.region_bounds);
    r_bar.crop_to_region       = true;
    r_bar.output_path          = input.out_dir + "/V_bar.png";

    PlotScalarFieldView(grid, r_full,   base_options);
    PlotScalarFieldView(grid, r_top,    base_options);
    PlotScalarFieldView(grid, r_bottom, base_options);
    PlotScalarFieldView(grid, r_bar,    base_options);

    // ---------------- |E| views ----------------
    StreamlineConfig stream_cfg; // placeholder, kept for signature

    // Start from the V-views and override only what's different.

    // Full |E|
    ScalarViewRequest r_fullE = r_full;
    r_fullE.scalar_name           = "Enorm";
    r_fullE.cbar_title            = "E [V/m]";
    r_fullE.title                 = "|E|";
    r_fullE.show_contours         = false;
    r_fullE.show_streamlines      = true;
    r_fullE.n_contours            = 0;
    r_fullE.output_path           = input.out_dir + "/Enorm_full.png";

    // Top |E|
    ScalarViewRequest r_topE = r_top;
    r_topE.scalar_name           = "Enorm";
    r_topE.cbar_title            = "E [V/m]";
    r_topE.title                 = "|E| (top stack)";
    r_topE.show_contours         = false;
    r_topE.show_streamlines      = true;
    r_topE.n_contours            = 0;
    r_topE.output_path           = input.out_dir + "/Enorm_top.png";

    // Bottom |E|
    ScalarViewRequest r_bottomE = r_bottom;
    r_bottomE.scalar_name           = "Enorm";
    r_bottomE.cbar_title            = "E [V/m]";
    r_bottomE.title                 = "|E| (bottom stack)";
    r_bottomE.show_contours         = false;
    r_bottomE.show_streamlines      = true;
    r_bottomE.n_contours            = 0;
    r_bottomE.output_path           = input.out_dir + "/Enorm_bottom.png";

    // Field cage |E|
    ScalarViewRequest r_barE = r_bar;
    r_barE.scalar_name           = "Enorm";
    r_barE.cbar_title            = "E [V/m]";
    r_barE.title                 = "|E| (Field Cage)";
    r_barE.show_contours         = false;
    r_barE.show_streamlines      = true;
    r_barE.n_contours            = 0;
    r_barE.output_path           = input.out_dir + "/Enorm_bar.png";

    PlotScalarFieldView(grid, r_fullE,   base_options, &stream_cfg);
    PlotScalarFieldView(grid, r_topE,    base_options, &stream_cfg);
    PlotScalarFieldView(grid, r_bottomE, base_options, &stream_cfg);
    PlotScalarFieldView(grid, r_barE,    base_options, &stream_cfg);

    std::cout << "Standard V / |E| views written to: " << input.out_dir << "\n";
}

static PlotInput PreparePlotInput(const char* raw_path)
{
    PlotInput R;

    // --------------------------
    // 1. Normalize input path
    // --------------------------
    std::string path = raw_path ? std::string(raw_path) : std::string();

    // Expand '~'
    if (!path.empty() && path[0] == '~')
    {
        const char* home = std::getenv("HOME");
        if (!home)
        {
            std::cerr << "[PreparePlotInput] ERROR: HOME environment variable not set.\n";
            std::exit(1);
        }
        path = std::string(home) + path.substr(1);
    }

    // Strip trailing slashes
    while (!path.empty() && path.back() == '/')
        path.pop_back();

    auto ends_with = [](const std::string& s, const std::string& suf) {
        return s.size() >= suf.size() &&
               s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
    };

    auto dirname_of = [](const std::string& s) -> std::string {
        auto pos = s.find_last_of('/');
        if (pos == std::string::npos) return ".";
        if (pos == 0) return "/";
        return s.substr(0, pos);
    };

    // --------------------------
    // 2. Determine VTU + base dir
    // --------------------------
    if (ends_with(path, ".pvd"))
    {
        R.base_dir = dirname_of(path);
        // This follows your previous convention; adjust if your directory
        // layout changed.
        R.vtu_file = R.base_dir + "/Cycle000000/proc000000.vtu";
    }
    else
    {
        R.vtu_file = path;
        R.base_dir = dirname_of(path);
    }

    std::cout << "[PreparePlotInput] Reading VTU: " << R.vtu_file << "\n";

    // --------------------------
    // 3. Load the grid
    // --------------------------
    R.grid = LoadGridFromVTU(R.vtu_file);
    if (!R.grid)
    {
        std::cerr << "[PreparePlotInput] ERROR: Failed to load grid from "
                  << R.vtu_file << "\n";
        std::exit(1);
    }

    // --------------------------
    // 4. Prepare output directory
    // --------------------------
    R.out_dir = R.base_dir + "/plots";

    {
        std::string cmd = "mkdir -p \"" + R.out_dir + "\"";
        int ret = std::system(cmd.c_str());
        if (ret != 0)
        {
            std::cerr << "[PreparePlotInput] ERROR: Failed to create output directory: "
                      << R.out_dir << "\n";
            std::exit(1);
        }
    }

    return R;
}


//------------------------- New Things ------------------------------ 


PlotViewLayout MakeDefaultPlotViewLayout(bool cbarHorizontal,
                                         bool separateCbarViewport)
{
    PlotViewLayout layout;
    layout.cbarHorizontal = cbarHorizontal;

    // For now each plot uses the full window.
    layout.plot = {0.0, 0.0, 1.0, 1.0};

    if (cbarHorizontal)
    {
        // Horizontal cbar at the bottom.
        const double titleHeight = 0.10;                         // top band
        const double cbarHeight  = separateCbarViewport ? 0.18 : 0.12;
        const double gap         = 0.02;

        double px0 = layout.plot.x0;
        double py0 = layout.plot.y0;
        double px1 = layout.plot.x1;
        double py1 = layout.plot.y1;

        // Title band across the top
        layout.title = { px0, py1 - titleHeight, px1, py1 };

        // Colorbar band at the bottom (short, centered)
        layout.cbar = {
            px0 + 0.10,
            py0 + 0.02,
            px1 - 0.10,
            py0 + 0.02 + cbarHeight
        };

        // Colorbar title band above the bar
        layout.cbarTitle = {
            layout.cbar.x0,
            layout.cbar.y1 + gap,
            layout.cbar.x1,
            layout.cbar.y1 + gap + 0.06
        };

        // Labels just below and above the bar (for min/max etc.)
        layout.cbarLabelBottomOrRight = {
            layout.cbar.x0,
            layout.cbar.y0 - 0.06,
            layout.cbar.x1,
            layout.cbar.y0
        };

        layout.cbarLabelTopOrLeft = {
            layout.cbar.x0,
            layout.cbar.y1,
            layout.cbar.x1,
            layout.cbar.y1 + 0.06
        };

        // Field occupies the remaining area between cbar and title
        layout.field = {
            px0,
            layout.cbar.y1 + gap,
            px1,
            layout.title.y0 - gap
        };
    }
    else
    {
        // Vertical cbar on the right.
        const double titleHeight = 0.10;
        const double cbarWidth   = separateCbarViewport ? 0.18 : 0.12;
        const double sideMargin  = 0.04;
        const double gap         = 0.02;

        double px0 = layout.plot.x0;
        double py0 = layout.plot.y0;
        double px1 = layout.plot.x1;
        double py1 = layout.plot.y1;

        // Title band across the top
        layout.title = { px0, py1 - titleHeight, px1, py1 };

        // Colorbar band on the right, leaving some margin at top/bottom
        layout.cbar = {
            px1 - cbarWidth,
            py0 + sideMargin,
            px1 - sideMargin,
            py1 - titleHeight - sideMargin
        };

        // Colorbar title band above the cbar
        layout.cbarTitle = {
            layout.cbar.x0,
            layout.cbar.y1 + gap,
            layout.cbar.x1,
            layout.cbar.y1 + gap + 0.06
        };

        // Labels to left/right of cbar
        layout.cbarLabelTopOrLeft = {
            layout.cbar.x0 - 0.06,
            layout.cbar.y0,
            layout.cbar.x0,
            layout.cbar.y1
        };

        layout.cbarLabelBottomOrRight = {
            layout.cbar.x1,
            layout.cbar.y0,
            layout.cbar.x1 + 0.06,
            layout.cbar.y1
        };

        // Field: everything left of cbar and below title
        layout.field = {
            px0,
            py0,
            layout.cbar.x0 - gap,
            layout.title.y0 - gap
        };
    }

    // Clamp to [0,1] just in case
    auto clampRect = [](ViewportRect& r) {
        r.x0 = std::max(0.0, std::min(1.0, r.x0));
        r.y0 = std::max(0.0, std::min(1.0, r.y0));
        r.x1 = std::max(0.0, std::min(1.0, r.x1));
        r.y1 = std::max(0.0, std::min(1.0, r.y1));
    };

    clampRect(layout.plot);
    clampRect(layout.field);
    clampRect(layout.cbar);
    clampRect(layout.title);
    clampRect(layout.cbarTitle);
    clampRect(layout.cbarLabelTopOrLeft);
    clampRect(layout.cbarLabelBottomOrRight);

    return layout;
}



}

int make_plots(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <simulation .vtu or Simulation.pvd>\n";
        return 1;
    }
    using namespace FEMPlot;

    // ---------------------------------------------------------------------
    // 1) Offscreen rendering backend (Mesa)
    // ---------------------------------------------------------------------
    vtkNew<vtkGraphicsFactory> graphics_factory;
    graphics_factory->SetOffScreenOnlyMode(1);
    graphics_factory->SetUseMesaClasses(1);

    // ---------------------------------------------------------------------
    // 2) Load grid and path handling
    // ---------------------------------------------------------------------
    // PreparePlotInput is already defined in FEMPlotVTK.cpp and returns
    // a PlotInput with { base_dir, vtu_file, grid, out_dir }.
    FEMPlot::PlotInput input = PreparePlotInput(argv[1]);

    if (!input.grid)
    {
        std::cerr << "ERROR: grid is null after PreparePlotInput.\n";
        return 1;
    }

    // ---------------------------------------------------------------------
    // 3) Base plotting options (image size, edges, etc.)
    // ---------------------------------------------------------------------
    FEMPlot::PlottingOptions opt;
    // You can tweak defaults here if desired, e.g.:
    opt.image_width  = 3840;
    opt.image_height = 2160;
    // opt.show_edges   = false;

    // ---------------------------------------------------------------------
    // 4) Generate standard views for V and |E|
    // ---------------------------------------------------------------------
    PlotStandardViewsForPotentialAndEnorm(input.grid, input, opt);

    std::cout << "Saved standard V / |E| plots to: " << input.out_dir << "\n";
    return 0;
}