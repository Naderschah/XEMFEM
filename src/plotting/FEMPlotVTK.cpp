// FEMPlotVTK.cpp
//
// Implementation of the FEMPlot API plus an example main().

#include "FEMPlotVTK.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// VTK includes
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
#include <cstring> 

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>


namespace FEMPlot
{
// ---------------------------------------------------------------------------
// Core helpers
// ---------------------------------------------------------------------------

vtkSmartPointer<vtkUnstructuredGrid>
LoadGridFromVTU(const std::string& filename)
{
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    return reader->GetOutput();
}

vtkSmartPointer<vtkDataSetMapper>
CreateScalarMapper(vtkUnstructuredGrid* grid,
                   const std::string& scalar_name,
                   double range_min,
                   double range_max)
{
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(grid);

    // Use named point-data array instead of relying on "active" scalars.
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(scalar_name.c_str());
    mapper->ScalarVisibilityOn();

    // If explicit range given, use it; otherwise use global range of the array.
    if (!std::isnan(range_min) && !std::isnan(range_max))
    {
        mapper->SetScalarRange(range_min, range_max);
    }
    else
    {
        auto pd = grid->GetPointData();
        if (auto arr = pd->GetArray(scalar_name.c_str()))
        {
            double r[2];
            arr->GetRange(r);
            mapper->SetScalarRange(r);
        }
        else
        {
            std::cerr << "Warning: scalar array '" << scalar_name
                      << "' not found on grid; mapper will have no scalars.\n";
            mapper->ScalarVisibilityOff();
        }
    }

    return mapper;
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

    return actor;
}

vtkSmartPointer<vtkRenderer>
CreateRendererWithActor(vtkSmartPointer<vtkActor> actor,
                        const double* background)
{
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);

    double bg[3] = {1.0, 1.0, 1.0};
    if (background)
    {
        bg[0] = background[0];
        bg[1] = background[1];
        bg[2] = background[2];
    }
    renderer->SetBackground(bg[0], bg[1], bg[2]);

    return renderer;
}

std::pair<double, double>
ComputeLocalScalarRange(vtkUnstructuredGrid* grid,
                        const std::string& scalar_name,
                        const double region_bounds[6])
{
    auto pd = grid->GetPointData();
    auto scalars = pd->GetArray(scalar_name.c_str());
    if (!scalars)
    {
        std::cerr << "Warning: scalar array '" << scalar_name
                  << "' not found; using global range.\n";
        double global_range[2] = {0.0, 0.0};
        pd->GetScalars() ? pd->GetScalars()->GetRange(global_range)
                         : void();
        return {global_range[0], global_range[1]};
    }

    double local_min = std::numeric_limits<double>::infinity();
    double local_max = -std::numeric_limits<double>::infinity();

    vtkIdType npts = grid->GetNumberOfPoints();
    for (vtkIdType i = 0; i < npts; ++i)
    {
        double p[3];
        grid->GetPoint(i, p);

        if (p[0] < region_bounds[0] || p[0] > region_bounds[1] ||
            p[1] < region_bounds[2] || p[1] > region_bounds[3] ||
            p[2] < region_bounds[4] || p[2] > region_bounds[5])
        {
            continue;
        }

        double val = scalars->GetComponent(i, 0);
        if (val < local_min) local_min = val;
        if (val > local_max) local_max = val;
    }

    if (!std::isfinite(local_min) || !std::isfinite(local_max))
    {
        double r[2];
        scalars->GetRange(r);
        return {r[0], r[1]};
    }

    return {local_min, local_max};
}

// ---------------------------------------------------------------------------
// High-level plot builders
// ---------------------------------------------------------------------------

vtkSmartPointer<vtkRenderer>
CreatePotentialRenderer(vtkUnstructuredGrid* grid,
                        const std::string& scalar_name)
{
    auto mapper = CreateScalarMapper(grid, scalar_name);
    auto actor  = CreateMeshActor(mapper, true, 1.0);
    auto renderer = CreateRendererWithActor(actor);

    double bounds[6];
    grid->GetBounds(bounds);
    renderer->ResetCamera(bounds);

    EnsureScalarBar(renderer, scalar_name);

    return renderer;
}

vtkSmartPointer<vtkRenderer>
CreateZoomedRenderer(vtkUnstructuredGrid* grid,
                     const std::string& scalar_name,
                     const double region_bounds[6])
{
    auto range = ComputeLocalScalarRange(grid, scalar_name, region_bounds);
    auto mapper = CreateScalarMapper(grid, scalar_name, range.first, range.second);
    auto actor  = CreateMeshActor(mapper, true, 1.0);
    auto renderer = CreateRendererWithActor(actor);

    // Zoom camera to region bounds
    double b[6];
    std::copy(region_bounds, region_bounds + 6, b);
    renderer->ResetCamera(b);

    EnsureScalarBar(renderer, scalar_name);

    return renderer;
}

// ---------------------------------------------------------------------------
// Decorations and overlays
// ---------------------------------------------------------------------------

void AddXYAxes(vtkRenderer* renderer,
               const double bounds[6],
               const std::string& x_label,
               const std::string& y_label,
               const std::string& z_label)
{
    if (!renderer)
        return;

    auto axes = vtkSmartPointer<vtkCubeAxesActor>::New();
    axes->SetBounds(bounds);
    axes->SetCamera(renderer->GetActiveCamera());

    axes->SetXTitle(x_label.c_str());
    axes->SetYTitle(y_label.c_str());
    axes->SetZTitle(z_label.c_str());

    axes->XAxisLabelVisibilityOn();
    axes->YAxisLabelVisibilityOn();
    axes->ZAxisLabelVisibilityOn();

    axes->SetFlyModeToClosestTriad();

    renderer->AddActor(axes);
}

void EnsureScalarBar(vtkRenderer* renderer,
                     const std::string& title)
{
    if (!renderer)
        return;

    // Check if scalar bar already exists
    auto actors2D = renderer->GetActors2D();
    actors2D->InitTraversal();
    while (auto a2d = actors2D->GetNextActor2D())
    {
        if (vtkScalarBarActor::SafeDownCast(a2d))
            return;
    }

    // Find a DataSetMapper with LUT
    vtkDataSetMapper* mapper_with_lut = nullptr;
    auto actors3D = renderer->GetActors();
    actors3D->InitTraversal();
    while (auto a = actors3D->GetNextActor())
    {
        if (auto m = vtkDataSetMapper::SafeDownCast(a->GetMapper()))
        {
            if (m->GetLookupTable())
            {
                mapper_with_lut = m;
                break;
            }
        }
    }

    if (!mapper_with_lut)
        return;

    auto scalar_bar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalar_bar->SetLookupTable(mapper_with_lut->GetLookupTable());
    if (!title.empty())
        scalar_bar->SetTitle(title.c_str());
    scalar_bar->SetNumberOfLabels(4);

    renderer->AddActor2D(scalar_bar);
}

void AddStreamlines(vtkUnstructuredGrid* grid,
                    vtkRenderer* renderer,
                    const Result& result,
                    const StreamlineConfig& cfg)
{
    if (!grid || !renderer)
        return;

    vtkIdType npts = grid->GetNumberOfPoints();
    if (static_cast<size_t>(npts) != result.E_at_points.size())
    {
        std::cerr << "Error: AddStreamlines - E_at_points size mismatch.\n";
        return;
    }

    // Attach vector field as point-data array
    auto vec_array = vtkSmartPointer<vtkDoubleArray>::New();
    vec_array->SetName("E");
    vec_array->SetNumberOfComponents(3);
    vec_array->SetNumberOfTuples(npts);

    for (vtkIdType i = 0; i < npts; ++i)
    {
        const auto& e = result.E_at_points[static_cast<size_t>(i)];
        double v[3] = { e[0], e[1], e[2] };
        vec_array->SetTuple(i, v);
    }

    auto pd = grid->GetPointData();
    pd->AddArray(vec_array);
    pd->SetActiveVectors("E");

    // Domain center and radius for seeding and tube scaling
    double bounds[6];
    grid->GetBounds(bounds);
    double center[3] = {
        0.5 * (bounds[0] + bounds[1]),
        0.5 * (bounds[2] + bounds[3]),
        0.5 * (bounds[4] + bounds[5])
    };
    double dx = bounds[1] - bounds[0];
    double dy = bounds[3] - bounds[2];
    double dz = bounds[5] - bounds[4];
    double radius = 0.5 * std::max({dx, dy, dz});

    auto seeds = vtkSmartPointer<vtkPointSource>::New();
    seeds->SetCenter(center);
    seeds->SetRadius(radius);
    seeds->SetNumberOfPoints(cfg.n_seeds);
    seeds->SetDistributionToUniform();

    auto integrator = vtkSmartPointer<vtkRungeKutta4>::New();

    auto stream = vtkSmartPointer<vtkStreamTracer>::New();
    stream->SetInputData(grid);
    stream->SetSourceConnection(seeds->GetOutputPort());
    stream->SetIntegrator(integrator);
    stream->SetIntegrationDirectionToBoth();
    stream->SetMaximumPropagation(cfg.max_propagation);
    stream->SetInitialIntegrationStep(cfg.initial_step);
    stream->SetMinimumIntegrationStep(cfg.min_step);
    stream->SetMaximumIntegrationStep(cfg.max_step);
    stream->SetMaximumNumberOfSteps(cfg.max_steps);

    auto tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputConnection(stream->GetOutputPort());
    tubes->SetRadius(cfg.tube_radius_rel * radius);
    tubes->SetNumberOfSides(8);
    tubes->CappingOn();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tubes->GetOutputPort());
    mapper->ScalarVisibilityOff();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.0, 0.0, 0.0);

    renderer->AddActor(actor);
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

void RenderOffscreenPNG(vtkRenderer* /*unused*/,
                        const std::string& /*filename*/,
                        int width,
                        int height)
{
    std::cerr << "=== Interactive VTK test: drawing y = x ===\n";

    // 1. Build a simple line from (0,0,0) to (1,1,0)
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    points->InsertNextPoint(1.0, 1.0, 0.0);

    auto lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(2);
    lines->InsertCellPoint(0);
    lines->InsertCellPoint(1);

    auto polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, 0.0, 0.0); // red line
    actor->GetProperty()->SetLineWidth(2.0);

    // 2. Renderer + window
    auto ren = vtkSmartPointer<vtkRenderer>::New();
    ren->AddActor(actor);
    ren->SetBackground(1.0, 1.0, 1.0);
    ren->ResetCamera();

    auto window = vtkSmartPointer<vtkRenderWindow>::New();
    window->AddRenderer(ren);
    window->SetSize(width, height);
    window->SetWindowName("VTK y = x minimal test");

    // 3. Interactor to show the window
    auto iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(window);

    std::cerr << "Calling Render()...\n";
    window->Render();
    std::cerr << "Render() returned, starting interactor...\n";

    iren->Initialize();
    iren->Start();  // <-- This should pop up a window and block until you close it.

    std::cerr << "Interactor finished, returning from RenderOffscreenPNG.\n";
}





void RenderThreeSideBySide(vtkRenderer* left,
                           vtkRenderer* middle,
                           vtkRenderer* right,
                           const std::string& filename,
                           int width,
                           int height)
{
    if (!left || !middle || !right)
        return;

    auto window = vtkSmartPointer<vtkRenderWindow>::New();
    window->SetSize(width, height);

    left->SetViewport(0.0,         0.0, 1.0 / 3.0, 1.0);
    middle->SetViewport(1.0 / 3.0,  0.0, 2.0 / 3.0, 1.0);
    right->SetViewport(2.0 / 3.0,   0.0, 1.0,       1.0);

    window->AddRenderer(left);
    window->AddRenderer(middle);
    window->AddRenderer(right);

    window->Render();

    int* sz = window->GetSize();
    int w = sz[0];
    int h = sz[1];
    if (w <= 0 || h <= 0)
    {
        std::cerr << "RenderThreeSideBySide: invalid window size (" << w
                  << "x" << h << ")\n";
        return;
    }

    unsigned char* rgb = window->GetPixelData(0, 0, w - 1, h - 1, 1);
    if (!rgb)
    {
        std::cerr << "RenderThreeSideBySide: GetPixelData returned null.\n";
        return;
    }

    auto image = vtkSmartPointer<vtkImageData>::New();
    image->SetDimensions(w, h, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3); // RGB

    unsigned char* dst =
        static_cast<unsigned char*>(image->GetScalarPointer());
    const int row_bytes = w * 3;

    for (int y = 0; y < h; ++y)
    {
        unsigned char* src_row = rgb + (h - 1 - y) * row_bytes;
        unsigned char* dst_row = dst + y * row_bytes;
        std::memcpy(dst_row, src_row, row_bytes);
    }

    delete [] rgb;

    auto writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(image);
    writer->Write();
}

} // namespace FEMPlot
int make_plots(int argc, char** argv)
{
    using namespace FEMPlot;

    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <simulation .vtu or Simulation.pvd>\n";
        return 1;
    }

    // ---------------------------------------------------------------------
    // 1. Normalize input path (expand ~, strip trailing /)
    // ---------------------------------------------------------------------
    std::string path = argv[1];

    // Expand "~" manually
    if (!path.empty() && path[0] == '~')
    {
        const char* home = std::getenv("HOME");
        if (!home)
        {
            std::cerr << "Error: HOME environment variable not set.\n";
            return 1;
        }
        path = std::string(home) + path.substr(1);
    }

    while (!path.empty() && path.back() == '/')
        path.pop_back();

    auto ends_with = [](const std::string& s, const std::string& suffix) {
        return s.size() >= suffix.size() &&
               s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    auto dirname_of = [](const std::string& s) -> std::string {
        auto pos = s.find_last_of('/');
        if (pos == std::string::npos) return ".";
        if (pos == 0) return "/";
        return s.substr(0, pos);
    };

    // ---------------------------------------------------------------------
    // 2. Determine which VTU file to load and where to put plots
    // ---------------------------------------------------------------------
    std::string vtu_file;
    std::string base_dir;  // where plots will go

    if (ends_with(path, ".pvd"))
    {
        // MFEM ParaViewDataCollection layout:
        //   <dir>/Simulation.pvd
        //   <dir>/Cycle000000/proc000000.vtu
        base_dir = dirname_of(path);  // e.g. .../Simulation
        vtu_file = base_dir + "/Cycle000000/proc000000.vtu";
    }
    else
    {
        // Treat the argument as the actual VTU path
        vtu_file = path;
        base_dir = dirname_of(path);
    }

    std::cout << "Reading VTU: " << vtu_file << "\n";

    auto grid = LoadGridFromVTU(vtu_file);
    if (!grid)
    {
        std::cerr << "Failed to load grid from " << vtu_file << "\n";
        return 1;
    }

    // Output directory "plots" next to the data file
    std::string out_dir = base_dir + "/plots";

    {
        std::string cmd = "mkdir -p \"" + out_dir + "\"";
        int ret = std::system(cmd.c_str());
        if (ret != 0)
        {
            std::cerr << "Failed to create output directory: "
                      << out_dir << "\n";
            return 1;
        }
    }

    // ---------------------------------------------------------------------
    // 3. Full-domain potential renderer (scalar "V")
    // ---------------------------------------------------------------------
    const std::string potential_name = "V";

    auto r_full = CreatePotentialRenderer(grid, potential_name);
    double full_bounds[6];
    grid->GetBounds(full_bounds);
    AddXYAxes(r_full, full_bounds, "x", "y", "z");

    // ---------------------------------------------------------------------
    // 4. Zoomed renderer (center 50%) with local color range
    // ---------------------------------------------------------------------
    double xmid = 0.5 * (full_bounds[0] + full_bounds[1]);
    double ymid = 0.5 * (full_bounds[2] + full_bounds[3]);
    double zmid = 0.5 * (full_bounds[4] + full_bounds[5]);

    double xhalf = 0.25 * (full_bounds[1] - full_bounds[0]);
    double yhalf = 0.25 * (full_bounds[3] - full_bounds[2]);
    double zhalf = 0.25 * (full_bounds[5] - full_bounds[4]);

    double zoom_bounds[6] = {
        xmid - xhalf, xmid + xhalf,
        ymid - yhalf, ymid + yhalf,
        zmid - zhalf, zmid + zhalf
    };

    auto r_zoom = CreateZoomedRenderer(grid, potential_name, zoom_bounds);
    AddXYAxes(r_zoom, zoom_bounds, "x", "y", "z");

    // ---------------------------------------------------------------------
    // 5. Full-domain + streamlines computed from potential V
    // ---------------------------------------------------------------------
    auto r_stream = CreatePotentialRenderer(grid, potential_name);
    AddXYAxes(r_stream, full_bounds, "x", "y", "z");

    // Use VTK to compute E = -∇V on the points, then feed that to AddStreamlines
    FEMPlot::Result res;
    res.E_at_points.resize(static_cast<size_t>(grid->GetNumberOfPoints()));

    // Gradient filter: input scalar "V" on points, output array "GradV"
    auto grad = vtkSmartPointer<vtkGradientFilter>::New();
    grad->SetInputData(grid);
    grad->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS,
                          potential_name.c_str());
    grad->SetResultArrayName("GradV");
    grad->Update();

    vtkDataSet* ds = grad->GetOutput();
    if (!ds)
    {
        std::cerr << "Error: vtkGradientFilter produced null output.\n";
        return 1;
    }

    auto pd = ds->GetPointData();
    auto gradArray = pd->GetArray("GradV");
    if (!gradArray || gradArray->GetNumberOfComponents() < 3)
    {
        std::cerr << "Error: gradient array 'GradV' not found or invalid.\n";
        return 1;
    }

    vtkIdType npts = ds->GetNumberOfPoints();
    if (npts != static_cast<vtkIdType>(res.E_at_points.size()))
    {
        std::cerr << "Error: gradient output point count mismatch.\n";
        return 1;
    }

    // Fill E(x) = -∇V into res.E_at_points
    for (vtkIdType i = 0; i < npts; ++i)
    {
        double g[3];
        gradArray->GetTuple(i, g); // g = ∇V
        res.E_at_points[static_cast<size_t>(i)] = {
            -g[0], -g[1], -g[2]
        };
    }

    StreamlineConfig cfg; // defaults
    AddStreamlines(grid, r_stream, res, cfg);

    // ---------------------------------------------------------------------
    // 6. Save plots
    // ---------------------------------------------------------------------
    RenderOffscreenPNG(r_full,
                       out_dir + "/potential_full.png",
                       1920, 1080);

    RenderOffscreenPNG(r_zoom,
                       out_dir + "/potential_zoom.png",
                       1920, 1080);

    RenderOffscreenPNG(r_stream,
                       out_dir + "/potential_streamlines.png",
                       1920, 1080);

    RenderThreeSideBySide(r_full, r_zoom, r_stream,
                          out_dir + "/three_plots.png",
                          2400, 800);

    std::cout << "Saved plots to: " << out_dir << "\n";
    return 0;
}
