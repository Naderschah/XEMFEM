// FEMPlotVTK.cpp
//
// Implementation of the FEMPlot API plus an example main().

#include "FEMPlotVTK.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <fstream>
#include <sstream>
#include <cctype>

// Mappers and Actors
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkAppendPolyData.h>

#include <vtkGridAxesActor3D.h>
#include <vtkGridAxesHelper.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCellData.h>


#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkMath.h>
#include <vtkStripper.h>

#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkActor2D.h>
#include <vtkLineSource.h>
#include <vtkCoordinate.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkActor2DCollection.h>
#include <vtkActorCollection.h>
#include <vtkBoundedPointSource.h>
#include <vtkBox.h>
#include <vtkProperty2D.h>
#include <vtkCamera.h>
#include <vtkCubeAxesActor.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkExtractGeometry.h>
#include <vtkGradientFilter.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRungeKutta4.h>
#include <vtkScalarBarActor.h>
#include <vtkStreamTracer.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTubeFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXMLUnstructuredGridReader.h>

namespace FEMPlot
{

// -----------------------------------------------------------------------------
// Plot config loading (simple INI-style)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Plot config loading (INI-style)
// -----------------------------------------------------------------------------

namespace
{

    static LayoutConfig g_layoutCfg;

    std::string GetSourceDir()
    {
        std::string f = __FILE__;
        std::size_t pos = f.find_last_of("/\\");
        if (pos == std::string::npos)
            return ".";
        return f.substr(0, pos);
    }
    inline void Trim(std::string& s)
    {
        auto is_space = [](unsigned char c){ return std::isspace(c); };
        auto it = std::find_if_not(s.begin(), s.end(), is_space);
        s.erase(s.begin(), it);
        auto rit = std::find_if_not(s.rbegin(), s.rend(), is_space);
        s.erase(rit.base(), s.end());
    }

    inline void StripComment(std::string& s)
    {
        std::size_t pos = s.find_first_of("#;");
        if (pos != std::string::npos)
            s.erase(pos);
    }

    inline bool ParseBool(const std::string& v, bool& out)
    {
        if (v.empty()) return false;
        if (v == "0") { out = false; return true; }
        if (v == "1") { out = true;  return true; }
        std::string low = v;
        std::transform(low.begin(), low.end(), low.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        if (low == "true" || low == "yes" || low == "on")
        {
            out = true; return true;
        }
        if (low == "false" || low == "no" || low == "off")
        {
            out = false; return true;
        }
        return false;
    }

    enum class SectionKind
    {
        None,
        FrameFull,
        FrameStack,
        FrameBar,
        ContentV,
        ContentE,
        Axes,
        Text,
        Regions,
        LayoutHorizontal,
        LayoutVertical
    };

    SectionKind SectionFromName(const std::string& name)
    {
        if (name == "frame.full")   return SectionKind::FrameFull;
        if (name == "frame.stack")  return SectionKind::FrameStack;
        if (name == "frame.bar")    return SectionKind::FrameBar;
        if (name == "content.V")    return SectionKind::ContentV;
        if (name == "content.E")    return SectionKind::ContentE;
        if (name == "axes")         return SectionKind::Axes;
        if (name == "text")         return SectionKind::Text;
        if (name == "regions")      return SectionKind::Regions;
        if (name == "layout.horizontal")return SectionKind::LayoutHorizontal;
        if (name == "layout.vertical")  return SectionKind::LayoutVertical;
        return SectionKind::None;
    }

} // anonymous namespace

bool LoadPlotConfig(const std::string& path, PlotConfig& cfg)
{
    std::ifstream in(path);
    if (!in)
    {
        std::cerr << "[LoadPlotConfig] INFO: could not open '" << path
                  << "', using compiled defaults.\n";
        return false;
    }

    SectionKind section = SectionKind::None;
    std::string line;

    while (std::getline(in, line))
    {
        StripComment(line);
        Trim(line);
        if (line.empty())
            continue;

        // Section header: [name]
        if (line.front() == '[' && line.back() == ']')
        {
            std::string sec = line.substr(1, line.size() - 2);
            Trim(sec);
            section = SectionFromName(sec);
            continue;
        }

        std::size_t eq = line.find('=');
        if (eq == std::string::npos)
            continue;

        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        Trim(key);
        Trim(val);
        if (key.empty() || val.empty())
            continue;

        std::istringstream iss(val);

        switch (section)
        {
        case SectionKind::FrameFull:
        case SectionKind::FrameStack:
        case SectionKind::FrameBar:
        {
            FrameConfig* fc = nullptr;
            if (section == SectionKind::FrameFull)  fc = &cfg.frame_full;
            if (section == SectionKind::FrameStack) fc = &cfg.frame_stack;
            if (section == SectionKind::FrameBar)   fc = &cfg.frame_bar;
            if (!fc) break;

            if (key == "image_width")
            {
                iss >> fc->image_width;
            }
            else if (key == "image_height")
            {
                iss >> fc->image_height;
            }
            else if (key == "zoom_factor")
            {
                iss >> fc->zoom_factor;
            }
            else if (key == "cbar_horizontal")
            {
                bool b; if (ParseBool(val, b)) fc->cbar_horizontal = b;
            }
            else if (key == "crop_to_region")
            {
                bool b; if (ParseBool(val, b)) fc->crop_to_region = b;
            }
            else if (key == "separate_cbar_viewport")
            {
                bool b; if (ParseBool(val, b)) fc->separate_cbar_viewport = b;
            }
            else if (key == "title_font_size")
            {
                iss >> fc->text.title_font_size;
            }
            else if (key == "axis_title_font_size")
            {
                iss >> fc->text.axis_title_font_size;
            }
            else if (key == "axis_label_font_size")
            {
                iss >> fc->text.axis_label_font_size;
            }
            else if (key == "cbar_title_font_size")
            {
                iss >> fc->text.cbar_title_font_size;
            }
            else if (key == "cbar_label_font_size")
            {
                iss >> fc->text.cbar_label_font_size;
            }
            else if (key == "cbar_minmax_font_size")
            {
                iss >> fc->text.cbar_minmax_font_size;
            }
            else if (key == "cbar_num_labels")
            {
                int iv = 0;
                iss >> iv;
                if (iv > 0) fc->text.cbar_num_labels = iv;
            }
            else if (key == "x_num_labels")
            {
                int iv = 0;
                iss >> iv;
                if (iv > 0) fc->x_num_labels = iv;
            }
            break;
        }

        case SectionKind::ContentV:
        case SectionKind::ContentE:
        {
            ContentConfig* cc =
                (section == SectionKind::ContentV) ? &cfg.content_V : &cfg.content_E;

            if (key == "scalar_name")
            {
                cc->scalar_name = val;
            }
            else if (key == "cbar_title")
            {
                cc->cbar_title = val;
            }
            else if (key == "color_map")
            {
                cc->color_map_name = val;
            }
            else if (key == "show_contours")
            {
                bool b; if (ParseBool(val, b)) cc->show_contours = b;
            }
            else if (key == "show_streamlines")
            {
                bool b; if (ParseBool(val, b)) cc->show_streamlines = b;
            }
            else if (key == "n_contours")
            {
                iss >> cc->n_contours;
            }
            else if (key == "title_full")
            {
                cc->title_full = val;
            }
            else if (key == "title_stack_top")
            {
                cc->title_stack_top = val;
            }
            else if (key == "title_stack_bottom")
            {
                cc->title_stack_bottom = val;
            }
            else if (key == "title_bar")
            {
                cc->title_bar = val;
            }
            break;
        }

        case SectionKind::Axes:
        {
            if (key == "x_label")      cfg.axes.x_label = val;
            else if (key == "y_label") cfg.axes.y_label = val;
            else if (key == "z_label") cfg.axes.z_label = val;
            break;
        }

        case SectionKind::Text:
        {
            int iv = 0;
            if (key == "title_font_size")
            {
                iss >> cfg.text.title_font_size;
            }
            else if (key == "axis_title_font_size")
            {
                iss >> cfg.text.axis_title_font_size;
            }
            else if (key == "axis_label_font_size")
            {
                iss >> cfg.text.axis_label_font_size;
            }
            else if (key == "cbar_title_font_size")
            {
                iss >> cfg.text.cbar_title_font_size;
            }
            else if (key == "cbar_label_font_size")
            {
                iss >> cfg.text.cbar_label_font_size;
            }
            else if (key == "cbar_minmax_font_size")
            {
                iss >> cfg.text.cbar_minmax_font_size;
            }
            else if (key == "cbar_num_labels")
            {
                iss >> iv;
                if (iv > 0) cfg.text.cbar_num_labels = iv;
            }
            break;
        }

        case SectionKind::Regions:
        {
            if (key == "top_ymin")      iss >> cfg.regions.top_ymin;
            else if (key == "top_ymax") iss >> cfg.regions.top_ymax;
            else if (key == "bottom_ymin") iss >> cfg.regions.bottom_ymin;
            else if (key == "bottom_ymax") iss >> cfg.regions.bottom_ymax;
            else if (key == "bar_dx")      iss >> cfg.regions.bar_dx;
            else if (key == "bar_ymin")    iss >> cfg.regions.bar_ymin;
            else if (key == "bar_ymax")    iss >> cfg.regions.bar_ymax;
            break;
        }
        case SectionKind::LayoutHorizontal:
        {
            // Horizontal layout: bottom→top bands + cbar horizontal extent
            if (key == "label_y0")      iss >> cfg.layout.h_label_y0;
            else if (key == "label_y1") iss >> cfg.layout.h_label_y1;
            else if (key == "cbar_y0")  iss >> cfg.layout.h_cbar_y0;
            else if (key == "cbar_y1")  iss >> cfg.layout.h_cbar_y1;
            else if (key == "ctitle_y0")iss >> cfg.layout.h_ctitle_y0;
            else if (key == "ctitle_y1")iss >> cfg.layout.h_ctitle_y1;
            else if (key == "field_y0") iss >> cfg.layout.h_field_y0;
            else if (key == "field_y1") iss >> cfg.layout.h_field_y1;
            else if (key == "title_y0") iss >> cfg.layout.h_title_y0;
            else if (key == "title_y1") iss >> cfg.layout.h_title_y1;
            else if (key == "cbar_x0")  iss >> cfg.layout.h_cbar_x0;
            else if (key == "cbar_x1")  iss >> cfg.layout.h_cbar_x1;
            break;
        }
        case SectionKind::LayoutVertical:
        {
            // Vertical layout: left/right + bottom→top bands
            if (key == "field_x0")         iss >> cfg.layout.v_field_x0;
            else if (key == "field_x1")    iss >> cfg.layout.v_field_x1;
            else if (key == "label_bot_y0")iss >> cfg.layout.v_label_bot_y0;
            else if (key == "label_bot_y1")iss >> cfg.layout.v_label_bot_y1;
            else if (key == "cbar_y0")     iss >> cfg.layout.v_cbar_y0;
            else if (key == "cbar_y1")     iss >> cfg.layout.v_cbar_y1;
            else if (key == "label_top_y0")iss >> cfg.layout.v_label_top_y0;
            else if (key == "label_top_y1")iss >> cfg.layout.v_label_top_y1;
            else if (key == "title_y0")    iss >> cfg.layout.v_title_y0;
            else if (key == "title_y1")    iss >> cfg.layout.v_title_y1;
            else if (key == "cbar_x0")     iss >> cfg.layout.v_cbar_x0;
            else if (key == "cbar_x1")     iss >> cfg.layout.v_cbar_x1;
            break;
        }
        default:
            break;
        }
    }

    std::cerr << "[LoadPlotConfig] Loaded plot config from '" << path << "'.\n";
    return true;
} // namespace


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

    auto result = vtkSmartPointer<vtkUnstructuredGrid>::New();
    result->ShallowCopy(ug);
    return result;
}

// -----------------------------------------------------------------------------
// Viewport layout
// -----------------------------------------------------------------------------

PlotViewLayout MakeDefaultPlotViewLayout(bool cbarHorizontal,
                                         bool /*separateCbarViewport*/)
{
    PlotViewLayout layout;
    layout.cbarHorizontal = cbarHorizontal;

    layout.plot = {0.0, 0.0, 1.0, 1.0};

    if (cbarHorizontal)
    {
        // Horizontal colorbar at the bottom; all bands from INI / defaults
        const double y0_label = g_layoutCfg.h_label_y0;
        const double y1_label = g_layoutCfg.h_label_y1;
        const double y0_cbar  = g_layoutCfg.h_cbar_y0;
        const double y1_cbar  = g_layoutCfg.h_cbar_y1;
        const double y0_ctit  = g_layoutCfg.h_ctitle_y0;
        const double y1_ctit  = g_layoutCfg.h_ctitle_y1;
        const double y0_field = g_layoutCfg.h_field_y0;
        const double y1_field = g_layoutCfg.h_field_y1;
        const double y0_title = g_layoutCfg.h_title_y0;
        const double y1_title = g_layoutCfg.h_title_y1;

        const double x_cbar0  = g_layoutCfg.h_cbar_x0;
        const double x_cbar1  = g_layoutCfg.h_cbar_x1;

        // bottom band: left half=min, right half=max
        layout.cbarLabelBottomOrRight = {0.0,  y0_label, 0.5,  y1_label}; // min
        layout.cbarLabelTopOrLeft     = {0.5,  y0_label, 1.0,  y1_label}; // max

        // colorbar itself
        layout.cbar = {x_cbar0, y0_cbar,
                       x_cbar1, y1_cbar};

        // colorbar title above the bar
        layout.cbarTitle = {x_cbar0, y0_ctit,
                            x_cbar1, y1_ctit};

        // field and plot title
        layout.field = {0.0,  y0_field,
                        1.0,  y1_field};

        layout.title = {0.0,  y0_title,
                        1.0,  y1_title};
    }
    else
    {
        // Vertical colorbar on the right; all bands from INI / defaults
        const double x_field0 = g_layoutCfg.v_field_x0;
        const double x_field1 = g_layoutCfg.v_field_x1;

        const double x_cbar0  = g_layoutCfg.v_cbar_x0;
        const double x_cbar1  = g_layoutCfg.v_cbar_x1;

        const double y0_label_bot = g_layoutCfg.v_label_bot_y0;
        const double y1_label_bot = g_layoutCfg.v_label_bot_y1;
        const double y0_cbar      = g_layoutCfg.v_cbar_y0;
        const double y1_cbar      = g_layoutCfg.v_cbar_y1;
        const double y0_label_top = g_layoutCfg.v_label_top_y0;
        const double y1_label_top = g_layoutCfg.v_label_top_y1;
        const double y0_titleBand = g_layoutCfg.v_title_y0;
        const double y1_titleBand = g_layoutCfg.v_title_y1;

        // field on the left
        layout.field = {x_field0, 0.0,
                        x_field1, 1.0};

        // colorbar
        layout.cbar = {x_cbar0, y0_cbar,
                       x_cbar1, y1_cbar};

        // labels below / above bar
        layout.cbarLabelBottomOrRight = {x_cbar0, y0_label_bot,
                                         x_cbar1, y1_label_bot}; // min
        layout.cbarLabelTopOrLeft     = {x_cbar0, y0_label_top,
                                         x_cbar1, y1_label_top}; // max

        // plot title across full width at top
        layout.title = {0.0,        y0_titleBand,
                        1.0,        y1_titleBand};

        // cbar title, aligned with cbar, in same top band
        // TODO Add config entry 
        layout.cbarTitle = {x_cbar0, y0_titleBand -0.13,
                            x_cbar1, y1_titleBand -0.13};
    }

    return layout;
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

    if (b[4] == b[5])
    {
        const double dz = 1e-6 * std::max({b[1] - b[0], b[3] - b[2], 1.0});
        b[4] -= dz * 0.5;
        b[5] += dz * 0.5;
    }

    vtkNew<vtkBox> box;
    box->SetBounds(b);

    vtkNew<vtkExtractGeometry> extract;
    extract->SetInputData(grid);
    extract->SetImplicitFunction(box);
    extract->ExtractInsideOn();
    extract->ExtractBoundaryCellsOn();
    extract->Update();

    auto clipped = vtkSmartPointer<vtkUnstructuredGrid>::New();
    clipped->ShallowCopy(extract->GetOutput());
    return clipped;
}

// -----------------------------------------------------------------------------
// Streamlines overlay
// -----------------------------------------------------------------------------

void AddStreamlines(vtkDataSet* dataset,
                    vtkRenderer* renderer,
                    const Result& result,
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

    // 1) Attach / update vector field E(x) on the dataset
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
        if (existingE->GetNumberOfComponents() == 3 &&
            existingE->GetNumberOfTuples() == npts)
        {
            evec = vtkDoubleArray::SafeDownCast(existingE);
        }
        else
        {
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

    // 2) Bounding box → characteristic radius for tube radius
    double bounds[6];
    dataset->GetBounds(bounds);
    const double dx = bounds[1] - bounds[0];
    const double dy = bounds[3] - bounds[2];
    const double dz = bounds[5] - bounds[4];
    const double radius = 0.5 * std::max({dx, dy, dz, 1e-6});

    // 3) Seed generator within dataset bounds
    auto seeds = vtkSmartPointer<vtkBoundedPointSource>::New();
    seeds->SetNumberOfPoints(cfg.n_seeds);
    seeds->SetBounds(bounds);

    // 4) Stream tracer
    auto rk4 = vtkSmartPointer<vtkRungeKutta4>::New();

    auto tracer = vtkSmartPointer<vtkStreamTracer>::New();
    tracer->SetInputData(dataset);
    tracer->SetSourceConnection(seeds->GetOutputPort());
    tracer->SetIntegrator(rk4);
    tracer->SetIntegrationDirectionToBoth();
    tracer->SetMaximumPropagation(cfg.max_propagation);
    tracer->SetInitialIntegrationStep(cfg.initial_step);
    tracer->SetMinimumIntegrationStep(cfg.min_step);
    tracer->SetMaximumIntegrationStep(cfg.max_step);
    tracer->SetMaximumNumberOfSteps(cfg.max_steps);
    tracer->SetComputeVorticity(false);

    // 5) Tubes
    auto tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputConnection(tracer->GetOutputPort());
    tubes->SetRadius(radius * cfg.tube_radius_rel);
    tubes->SetNumberOfSides(10);
    tubes->CappingOn();

    // 6) Mapper + Actor
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tubes->GetOutputPort());
    mapper->ScalarVisibilityOff();

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

// Render an already-configured render window.
static void RenderOffscreenPNG(vtkRenderWindow* renderWindow,
                               const std::string& filename)
{
    if (!renderWindow)
    {
        std::cerr << "RenderOffscreenPNG(window) ERROR: renderWindow is null\n";
        return;
    }

    renderWindow->SetOffScreenRendering(1);
    renderWindow->Render();

    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    windowToImageFilter->SetInput(renderWindow);
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

// Wrapper for single-renderer usage.
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

    RenderOffscreenPNG(renderWindow, filename);

    renderWindow->RemoveRenderer(renderer);
}

// -----------------------------------------------------------------------------
// Mesh/axes helpers
// -----------------------------------------------------------------------------

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

inline double NiceStep(double x)
{
    if (x <= 0.0)
        return 0.0;

    const double expv = std::floor(std::log10(x));
    const double f    = x / std::pow(10.0, expv);

    double nf;
    if (f < 1.5)
        nf = 1.0;
    else if (f < 3.0)
        nf = 2.0;
    else if (f < 7.0)
        nf = 5.0;
    else
        nf = 10.0;

    return nf * std::pow(10.0, expv);
}
inline void ConfigureAxisTicks(vtkGridAxesActor3D* axes,
                               const double bounds[6],
                               int axis,
                               int requestedTicks)
{
    // Do not trust the docs its all wrong just look at the source code - also function names are nonsensical
    if (!axes || requestedTicks <= 0)
        return;
    if (axis < 0 || axis > 2)
        return;

    const int idx = 2 * axis;
    double minVal = bounds[idx];
    double maxVal = bounds[idx + 1];

    if (minVal > maxVal)
        std::swap(minVal, maxVal);

    const double span = maxVal - minVal;
    if (span <= 0.0)
    {
        // Degenerate: single label at minVal
        axes->SetUseCustomLabels(axis, true);
        axes->SetNumberOfLabels(axis, 1);
        axes->SetLabel(axis, 0, minVal);
        return;
    }

    // 1) Rough step from requested tick count
    if (requestedTicks < 2)
        requestedTicks = 2;

    const double rawStep = span / static_cast<double>(requestedTicks - 1);
    double step = NiceStep(rawStep);
    if (step <= 0.0)
        step = rawStep; // fallback

    // 2) Snap to an inner range [start, end] ⊆ [minVal, maxVal]
    const double start = std::ceil(minVal / step) * step;
    const double end   = std::floor(maxVal / step) * step;

    if (end < start)
    {
        // Step is too large; fall back to a single label at center.
        const double mid = 0.5 * (minVal + maxVal);
        axes->SetUseCustomLabels(axis, true);
        axes->SetNumberOfLabels(axis, 1);
        axes->SetLabel(axis, 0, mid);
        return;
    }

    // 3) Compute actual number of labels and fill them.
    const double nStepsReal = (end - start) / step;
    const vtkIdType nLabels = static_cast<vtkIdType>(std::floor(nStepsReal + 0.5)) + 1;

    axes->SetUseCustomLabels(axis, true);
    axes->SetNumberOfLabels(axis, nLabels);

    for (vtkIdType i = 0; i < nLabels; ++i)
    {
        const double val = start + static_cast<double>(i) * step;
        axes->SetLabel(axis, i, val);
    }
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

    auto axes = vtkSmartPointer<vtkGridAxesActor3D>::New();

    // Use your physical region as grid bounds
    double b[6] = {
        bounds[0], bounds[1],
        bounds[2], bounds[3],
        bounds[4], bounds[5]
    };

    // If this is effectively 2D in Z, give it a tiny thickness
    if (b[4] == b[5])
    {
        const double dz = 1e-6 * std::max({b[1] - b[0], b[3] - b[2], 1.0});
        b[4] -= 0.5 * dz;
        b[5] += 0.5 * dz;
    }

    axes->SetGridBounds(b);

    // Only a single XY face (2D view) – no extra axes floating around.
    axes->SetFaceMask(vtkGridAxesHelper::MIN_XY);

    // We want ONLY axes/edges/ticks, NO interior grid planes.
    axes->SetGenerateGrid(false);   // disable interior grid
    axes->SetGenerateEdges(true);   // keep the box edges
    axes->SetGenerateTicks(true);   // keep tick markers

    // Make axes lines clearly visible (black on white background)
    axes->GetProperty()->SetColor(0.0, 0.0, 0.0);
    axes->GetProperty()->SetLineWidth(1.0);

    // Titles: only X and Y; hide Z to avoid a third axis in a 2D plot.
    axes->SetXTitle(x_label);
    axes->SetYTitle(y_label + "  ");
    axes->SetZTitle("");  // hide Z-axis title completely

    // Reduce Tick frequency - 
    ConfigureAxisTicks(axes, b, /*axis=X*/ 0, /*requestedTicks*/ 5);
    ConfigureAxisTicks(axes, b, /*axis=Y*/ 1, /*requestedTicks*/ 10);
    

    // Text properties
    for (int i = 0; i < 3; ++i)
    {
        if (auto tprop = axes->GetTitleTextProperty(i))
        {
            tprop->SetFontFamilyToArial();
            tprop->SetColor(0.0, 0.0, 0.0);
            if (axis_title_font_size > 0)
                tprop->SetFontSize(axis_title_font_size);
        }
        if (auto lprop = axes->GetLabelTextProperty(i))
        {
            lprop->SetFontFamilyToArial();
            lprop->SetColor(0.0, 0.0, 0.0);
            if (axis_label_font_size > 0)
                lprop->SetFontSize(axis_label_font_size);
        }
    }

    // Labels only on bottom X and left Y edges:
    // use vtkGridAxesHelper::LabelMasks bits: MIN_X and MIN_Y.
    unsigned int labelMask = 0;
    labelMask |= vtkGridAxesHelper::MIN_X;  // left vertical edge
    labelMask |= vtkGridAxesHelper::MIN_Y;  // bottom horizontal edge
    axes->SetLabelMask(labelMask);

    // Ensure we don’t suppress those edges’ labels
    axes->SetLabelUniqueEdgesOnly(true);

    renderer->AddActor(axes);
}









// -----------------------------------------------------------------------------
// Scalar statistics + LUT
// -----------------------------------------------------------------------------

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

// Remove closed loops: lines where first and last point index are equal.
inline vtkSmartPointer<vtkPolyData>
ExtractOpenContours(vtkPolyData* input)
{
    if (!input)
        return nullptr;

    // Make sure we have polylines instead of lots of tiny segments.
    auto stripper = vtkSmartPointer<vtkStripper>::New();
    stripper->SetInputData(input);
    stripper->Update();

    vtkPolyData* stripped = stripper->GetOutput();
    if (!stripped)
        return nullptr;

    auto output = vtkSmartPointer<vtkPolyData>::New();
    output->ShallowCopy(stripped);  // copy points, arrays, etc.

    vtkCellArray* lines = stripped->GetLines();
    if (!lines)
    {
        output->SetLines(nullptr);
        return output;
    }

    auto newLines = vtkSmartPointer<vtkCellArray>::New();

    lines->InitTraversal();
    vtkIdType npts;
    const vtkIdType* pts = nullptr;
    while (lines->GetNextCell(npts, pts))
    {
        if (npts < 2)
            continue;

        // Closed polyline: first and last vertex are the same
        if (pts[0] == pts[npts - 1])
            continue;  // drop closed loop

        newLines->InsertNextCell(npts, pts);
    }

    output->SetLines(newLines);
    return output;
}

// Basic LUT builder.
static vtkSmartPointer<vtkLookupTable>
BuildLookupTable(const std::string& color_map_name,
                 double range_min,
                 double range_max)
{
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues(256);
    lut->SetRange(range_min, range_max);

    if (color_map_name == "gray" || color_map_name == "grey")
    {
        lut->SetHueRange(0.0, 0.0);
        lut->SetSaturationRange(0.0, 0.0);
        lut->SetValueRange(0.0, 1.0);
    }
    else if (color_map_name == "blue-red" || color_map_name == "coolwarm")
    {
        lut->SetHueRange(0.6667, 0.0);
        lut->SetSaturationRange(1.0, 1.0);
        lut->SetValueRange(1.0, 1.0);
    }
    else if (color_map_name == "viridis")
    {
        lut->SetHueRange(0.7, 0.1);
        lut->SetSaturationRange(1.0, 1.0);
        lut->SetValueRange(0.3, 1.0);
    }
    else
    {
        lut->SetHueRange(0.6667, 0.0);
        lut->SetSaturationRange(1.0, 1.0);
        lut->SetValueRange(1.0, 1.0);
    }

    lut->Build();
    return lut;
}

static vtkSmartPointer<vtkDoubleArray>
MakeNiceScalarBarLabels(double vmin, double vmax, int n_desired)
{
    auto labels = vtkSmartPointer<vtkDoubleArray>::New();

    if (!std::isfinite(vmin) || !std::isfinite(vmax) || vmin == vmax)
    {
        labels->InsertNextValue(vmin);
        return labels;
    }

    if (n_desired < 2)
        n_desired = 2;

    if (vmax < vmin)
        std::swap(vmin, vmax);

    double range = vmax - vmin;

    auto niceNum = [](double x, bool round) -> double {
        double expv  = std::floor(std::log10(x));
        double frac  = x / std::pow(10.0, expv);
        double niceFrac;

        if (round)
        {
            if (frac < 1.5)      niceFrac = 1.0;
            else if (frac < 3.0) niceFrac = 2.0;
            else if (frac < 7.0) niceFrac = 5.0;
            else                 niceFrac = 10.0;
        }
        else
        {
            if (frac <= 1.0)      niceFrac = 1.0;
            else if (frac <= 2.0) niceFrac = 2.0;
            else if (frac <= 5.0) niceFrac = 5.0;
            else                  niceFrac = 10.0;
        }

        return niceFrac * std::pow(10.0, expv);
    };

    // "Nice" step size
    double step = niceNum(range / (n_desired - 1), true);

    // Extend to nice min/max that still cover [vmin, vmax]
    double nice_min = std::floor(vmin / step) * step;
    double nice_max = std::ceil (vmax / step) * step;

    // Build labels
    for (double v = nice_min; v <= nice_max + 0.5 * step; v += step)
    {
        double value = v;

        // If value is exactly zero, shift it to a small positive value
        if (std::abs(value) < 1e-12)
            value = 0.1 * step;

        labels->InsertNextValue(value);
    }

    return labels;
}

inline double ComputeCellLength(vtkCell* cell, vtkPoints* pts)
{
    const vtkIdType npts = cell->GetNumberOfPoints();
    if (npts < 2 || !pts)
        return 0.0;

    double p0[3], p1[3];
    double length = 0.0;

    for (vtkIdType i = 0; i < npts - 1; ++i)
    {
        pts->GetPoint(cell->GetPointId(i),     p0);
        pts->GetPoint(cell->GetPointId(i + 1), p1);
        length += std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
    }

    return length;
}

// Keep only line cells whose length >= minRelativeLength * maxLength.
// minRelativeLength in (0,1], e.g. 0.2 = keep >= 20% of longest.
inline vtkSmartPointer<vtkPolyData>
FilterContoursByRelativeLength(vtkPolyData* input, double minRelativeLength)
{
    if (!input || minRelativeLength <= 0.0)
        return nullptr;

    vtkPoints* pts = input->GetPoints();
    vtkCellArray* inLines = input->GetLines();

    if (!pts || !inLines)
    {
        auto empty = vtkSmartPointer<vtkPolyData>::New();
        empty->ShallowCopy(input);
        empty->SetLines(nullptr);
        return empty;
    }

    // First pass: compute lengths and maximum.
    std::vector<double> lengths;
    lengths.reserve(inLines->GetNumberOfCells());

    inLines->InitTraversal();
    vtkIdType npts;
    const vtkIdType* cellPts = nullptr;
    double maxLen = 0.0;

    while (inLines->GetNextCell(npts, cellPts))
    {
        if (npts < 2)
        {
            lengths.push_back(0.0);
            continue;
        }

        double length = 0.0;
        double p0[3], p1[3];

        for (vtkIdType i = 0; i < npts - 1; ++i)
        {
            pts->GetPoint(cellPts[i],     p0);
            pts->GetPoint(cellPts[i + 1], p1);
            length += std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
        }

        lengths.push_back(length);
        if (length > maxLen)
            maxLen = length;
    }

    auto output = vtkSmartPointer<vtkPolyData>::New();
    output->ShallowCopy(input);

    if (maxLen <= 0.0)
    {
        output->SetLines(nullptr);
        return output;
    }

    const double lengthThreshold = minRelativeLength * maxLen;

    // Second pass: copy only long-enough cells.
    auto outLines = vtkSmartPointer<vtkCellArray>::New();
    inLines->InitTraversal();

    size_t idx = 0;
    while (inLines->GetNextCell(npts, cellPts))
    {
        double len = (idx < lengths.size() ? lengths[idx] : 0.0);
        if (len >= lengthThreshold && npts > 1)
        {
            outLines->InsertNextCell(npts, cellPts);
        }
        ++idx;
    }

    output->SetLines(outLines);
    return output;
}



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
    scalar_bar->SetLabelFormat("%.0f");

    auto tprop = scalar_bar->GetTitleTextProperty();
    tprop->SetFontFamilyToArial();
    tprop->SetFontSize(title_font_size);
    tprop->SetColor(0.0, 0.0, 0.0);

    auto lprop = scalar_bar->GetLabelTextProperty();
    lprop->SetFontFamilyToArial();
    lprop->SetFontSize(label_font_size);
    lprop->SetColor(0.0, 0.0, 0.0);

    // Its bold for some reason
    lprop->BoldOff();
    lprop->ItalicOff();
    lprop->ShadowOff();
    // Generic white background behind the bar + labels
    scalar_bar->DrawBackgroundOn();
    scalar_bar->GetBackgroundProperty()->SetColor(1.0, 1.0, 1.0);
    scalar_bar->GetBackgroundProperty()->SetOpacity(1.0);
    // Optional: drop the frame if you don't want a box
    scalar_bar->DrawFrameOff();

    // Default geometry (overridden in PlotScalarFieldView)
    scalar_bar->SetWidth(0.1);
    scalar_bar->SetHeight(0.7);
    scalar_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalar_bar->GetPositionCoordinate()->SetValue(0.90, 0.15);

    return scalar_bar;
}

// -----------------------------------------------------------------------------
// Gradient-based E and |E| attachment
// -----------------------------------------------------------------------------

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

        double e[3] = { -g[0], -g[1], -g[2] };
        evec->SetTuple(i, e);
    }

    pd->AddArray(evec);
    pd->SetActiveVectors(e_array_name.c_str());

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

// -----------------------------------------------------------------------------
// Main scalar-field view (viewport-based)
// -----------------------------------------------------------------------------

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

    // 1) Scalar statistics and LUT
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

    // 2) Mapper and mesh actor (with optional clipping planes)
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(dataset_to_plot);
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(request.scalar_name.c_str());
    mapper->SetScalarRange(stats.min_used, stats.max_used);
    mapper->SetLookupTable(lut);
    mapper->ScalarVisibilityOn();

    // Clean way to not plot out of bounds of interest 
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

    // 3) Render window + layout
    int img_w = (request.image_width  > 0) ? request.image_width  : plotting.image_width;
    int img_h = (request.image_height > 0) ? request.image_height : plotting.image_height;

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetOffScreenRendering(1);
    renderWindow->SetAlphaBitPlanes(0);
    renderWindow->SetMultiSamples(0);
    renderWindow->SetSize(img_w, img_h);
    renderWindow->SetNumberOfLayers(2);

    // Background renderer: full window, pure white
    auto bgRenderer = vtkSmartPointer<vtkRenderer>::New();
    bgRenderer->SetLayer(0);  // background layer
    bgRenderer->SetViewport(0.0, 0.0, 1.0, 1.0);
    bgRenderer->SetBackground(1.0, 1.0, 1.0);
    bgRenderer->SetBackground2(1.0, 1.0, 1.0);
    bgRenderer->GradientBackgroundOff();
    renderWindow->AddRenderer(bgRenderer);


    PlotViewLayout layout = MakeDefaultPlotViewLayout(
        request.cbar_horizontal,
        request.separate_cbar_viewport);

    auto autoFont = [&](int requested,
                        const ViewportRect& r,
                        double fracOfHeight) -> int
    {
        if (requested > 0)
            return requested;

        double hnorm = (r.y1 - r.y0);
        int h = static_cast<int>(hnorm * img_h + 0.5);
        if (h < 1) h = 1;

        int fs = static_cast<int>(h * fracOfHeight + 0.5);
        if (fs < 6)  fs = 6;
        if (fs > 72) fs = 72;
        return fs;
    };

    int titleFS      = autoFont(request.title_font_size,       layout.title,               0.50);
    int axisTitleFS  = autoFont(request.axis_title_font_size,  layout.field,               0.08);
    int axisLabelFS  = autoFont(request.axis_label_font_size,  layout.field,               0.06);
    int cbarTitleFS  = autoFont(request.cbar_title_font_size,  layout.cbarTitle,           0.50);
    int cbarLabelFS  = autoFont(request.cbar_label_font_size,  layout.cbar,                0.25);
    int cbarMinMaxFS = autoFont(request.cbar_minmax_font_size, layout.cbarLabelBottomOrRight, 0.45);

    auto makeRenderer = [&](const ViewportRect& r) {
        auto ren = vtkSmartPointer<vtkRenderer>::New();
        ren->SetBackground(1.0, 1.0, 1.0);
        ren->SetBackground2(1.0, 1.0, 1.0);
        ren->GradientBackgroundOff();
        ren->SetViewport(r.x0, r.y0, r.x1, r.y1);
        renderWindow->AddRenderer(ren);
        return ren;
    };

    vtkSmartPointer<vtkRenderer> field_renderer          = makeRenderer(layout.field);
    vtkSmartPointer<vtkRenderer> title_renderer          = makeRenderer(layout.title);
    vtkSmartPointer<vtkRenderer> cbar_renderer           = makeRenderer(layout.cbar);
    vtkSmartPointer<vtkRenderer> cbar_title_renderer     = makeRenderer(layout.cbarTitle);
    vtkSmartPointer<vtkRenderer> cbar_label_top_renderer = makeRenderer(layout.cbarLabelTopOrLeft);
    vtkSmartPointer<vtkRenderer> cbar_label_bot_renderer = makeRenderer(layout.cbarLabelBottomOrRight);
    
    field_renderer->AddActor(mesh_actor);

    // 4) Axes
    AddXYAxes(field_renderer, request.region_bounds,
            request.x_label, request.y_label, request.z_label,
            request.axis_title_font_size,
            request.axis_label_font_size);




    // 5) Plot title
    if (!request.title.empty())
    {
        auto titleActor = vtkSmartPointer<vtkTextActor>::New();
        titleActor->SetInput(request.title.c_str());

        auto tprop = titleActor->GetTextProperty();
        tprop->SetFontFamilyToArial();
        tprop->SetFontSize(titleFS);
        tprop->SetColor(0.0, 0.0, 0.0);
        tprop->SetJustificationToCentered();
        tprop->SetVerticalJustificationToCentered();

        titleActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        titleActor->SetPosition(0.5, 0.5);

        title_renderer->AddViewProp(titleActor);
    }

    // 6) Scalar bar
    auto scalar_bar = CreateScalarBar(
        lut,
        "",
        cbarTitleFS,
        cbarLabelFS,
        request.cbar_num_labels);

    if (request.cbar_horizontal)
    {
        // Horizontal bar: fill the entire cbar viewport
        scalar_bar->SetOrientationToHorizontal();
        scalar_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        scalar_bar->GetPositionCoordinate()->SetValue(0.0, 0.0);
        scalar_bar->SetWidth(1.0);
        scalar_bar->SetHeight(1.0);
    }
    else
    {
        // Vertical bar: fill the entire cbar viewport
        scalar_bar->SetOrientationToVertical();
        scalar_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        scalar_bar->GetPositionCoordinate()->SetValue(0.0, 0.0);
        scalar_bar->SetWidth(1.0);
        scalar_bar->SetHeight(1.0);
    }

    // Decide how many contour/label levels we want
    int n_for_levels = (request.n_contours > 0)
                        ? request.n_contours
                        : (request.cbar_num_labels > 0 ? request.cbar_num_labels : 5);

    // Build "nice" rounded levels over the clipped range
    auto levelValues = MakeNiceScalarBarLabels(
        stats.min_used,
        stats.max_used,
        n_for_levels);

    // Drive scalar bar labels from these values
    if (levelValues && levelValues->GetNumberOfTuples() > 0)
    {
        scalar_bar->SetNumberOfLabels(levelValues->GetNumberOfTuples());
        scalar_bar->SetCustomLabels(levelValues);
        scalar_bar->UseCustomLabelsOn();
    }
    else
    {
        scalar_bar->UseCustomLabelsOff();
    }

    cbar_renderer->AddViewProp(scalar_bar);

    // 7) Contour lines in the field view (same levels as scalar bar)
    if (request.show_contours && levelValues && levelValues->GetNumberOfTuples() > 0)
    {
        auto surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surface->SetInputData(dataset_to_plot);

        // Append polydata to collect cleaned contours from all levels
        auto appendAll = vtkSmartPointer<vtkAppendPolyData>::New();

        const int nLevels = static_cast<int>(levelValues->GetNumberOfTuples());

        for (int i = 0; i < nLevels; ++i)
        {
            double v = levelValues->GetValue(i);

            // Contour this single level
            auto contour = vtkSmartPointer<vtkContourFilter>::New();
            contour->SetInputConnection(surface->GetOutputPort());
            contour->SetInputArrayToProcess(
                0, 0, 0,
                vtkDataObject::FIELD_ASSOCIATION_POINTS,
                request.scalar_name.c_str());
            contour->SetNumberOfContours(1);
            contour->SetValue(0, v);

            // Join small line segments into polylines
            auto stripper = vtkSmartPointer<vtkStripper>::New();
            stripper->SetInputConnection(contour->GetOutputPort());
            stripper->JoinContiguousSegmentsOn();   // important
            // no need to Update(); appendAll will drive it

            appendAll->AddInputConnection(stripper->GetOutputPort());
        }

        appendAll->Update();

        // Filter short contours - loops and what not
        vtkPolyData* allContours = appendAll->GetOutput();
        auto filteredContours =
            FilterContoursByRelativeLength(allContours, /*minRelativeLength*/ 0.001); // keep >= 20% of longest

        auto contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        contourMapper->SetInputData(filteredContours);
        contourMapper->ScalarVisibilityOff();   // uniform color for all lines

        if (request.crop_to_region)
        {
            double xmin = request.region_bounds[0];
            double xmax = request.region_bounds[1];
            double ymin = request.region_bounds[2];
            double ymax = request.region_bounds[3];

            contourMapper->RemoveAllClippingPlanes();

            vtkNew<vtkPlane> pxMin;
            pxMin->SetOrigin(xmin, 0.0, 0.0);
            pxMin->SetNormal(1.0, 0.0, 0.0);
            contourMapper->AddClippingPlane(pxMin);

            vtkNew<vtkPlane> pxMax;
            pxMax->SetOrigin(xmax, 0.0, 0.0);
            pxMax->SetNormal(-1.0, 0.0, 0.0);
            contourMapper->AddClippingPlane(pxMax);

            vtkNew<vtkPlane> pyMin;
            pyMin->SetOrigin(0.0, ymin, 0.0);
            pyMin->SetNormal(0.0, 1.0, 0.0);
            contourMapper->AddClippingPlane(pyMin);

            vtkNew<vtkPlane> pyMax;
            pyMax->SetOrigin(0.0, ymax, 0.0);
            pyMax->SetNormal(0.0, -1.0, 0.0);
            contourMapper->AddClippingPlane(pyMax);
        }

        auto contourActor = vtkSmartPointer<vtkActor>::New();
        contourActor->SetMapper(contourMapper);

        // Visually distinct from viridis
        contourActor->GetProperty()->SetColor(1.0, 0.0, 1.0); // magenta
        contourActor->GetProperty()->SetLineWidth(1.0);
        contourActor->GetProperty()->SetOpacity(1.0);

        field_renderer->AddActor(contourActor);
    }

    // 7) Colorbar title actor
    if (!request.cbar_title.empty())
    {
        auto cbarTitleActor = vtkSmartPointer<vtkTextActor>::New();
        cbarTitleActor->SetInput(request.cbar_title.c_str());

        auto tprop = cbarTitleActor->GetTextProperty();
        tprop->SetFontFamilyToArial();
        tprop->SetBold(1);
        tprop->SetFontSize(cbarTitleFS);
        tprop->SetColor(0.0, 0.0, 0.0);
        tprop->SetJustificationToCentered();
        tprop->SetVerticalJustificationToCentered();

        cbarTitleActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        cbarTitleActor->SetPosition(0.5, 0.5);

        cbar_title_renderer->AddViewProp(cbarTitleActor);
    }

    // 8) True min/max annotations
    if (stats.min_used > stats.min_roi)
    {
        auto txt = vtkSmartPointer<vtkTextActor>::New();
        std::ostringstream ss;
        ss << "min = " << stats.min_roi;
        txt->SetInput(ss.str().c_str());

        auto tp = txt->GetTextProperty();
        tp->SetFontFamilyToArial();
        tp->SetFontSize(cbarMinMaxFS);
        tp->SetColor(0.0, 0.0, 0.0);
        tp->SetJustificationToCentered();
        tp->SetVerticalJustificationToCentered();

        txt->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        txt->SetPosition(0.5, 0.5);

        cbar_label_bot_renderer->AddViewProp(txt);
    }

    if (stats.max_used < stats.max_roi)
    {
        auto txt = vtkSmartPointer<vtkTextActor>::New();
        std::ostringstream ss;
        ss << "max = " << stats.max_roi;
        txt->SetInput(ss.str().c_str());

        auto tp = txt->GetTextProperty();
        tp->SetFontFamilyToArial();
        tp->SetFontSize(cbarMinMaxFS);
        tp->SetColor(0.0, 0.0, 0.0);
        tp->SetJustificationToCentered();
        tp->SetVerticalJustificationToCentered();

        txt->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        txt->SetPosition(0.5, 0.5);

        cbar_label_top_renderer->AddViewProp(txt);
    }

    // 9) Streamlines — not yet wired into new pipeline
    if (request.show_streamlines)
    {
        std::cerr << "[PlotScalarFieldView] NOTE: show_streamlines requested "
                     "but streamlines are not yet wired into the new view pipeline.\n";
    }

    // 10) Camera setup and zoom
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

    // 11) Render and write PNG
    RenderOffscreenPNG(renderWindow, request.output_path);
}

// -----------------------------------------------------------------------------
// High-level standard views
// -----------------------------------------------------------------------------
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

    // ---------------------------------------------------------------------
    // Detect which scalar array to treat as the potential V
    // ---------------------------------------------------------------------
    vtkPointData* pd = grid->GetPointData();
    if (!pd)
    {
        std::cerr << "[PlotStandardViewsForPotentialAndEnorm] ERROR: no point data on grid.\n";
        return;
    }

    std::string potentialName;
    int scalarCandidateCount = 0;
    std::string lastScalarCandidateName;

    const int nArrays = pd->GetNumberOfArrays();
    for (int i = 0; i < nArrays; ++i)
    {
        vtkDataArray* arr = pd->GetArray(i);
        if (!arr)
            continue;

        if (arr->GetNumberOfComponents() != 1)
            continue; // only scalar arrays

        const char* cname = arr->GetName();
        if (!cname || !*cname)
            continue;

        std::string name(cname);
        ++scalarCandidateCount;
        lastScalarCandidateName = name;

        if (name == "V")
            potentialName = "V";
        else if (name == "Electric_potential")
            potentialName = "Electric_potential";
    }

    if (potentialName.empty() && scalarCandidateCount == 1)
    {
        // Only one scalar present → assume it's the potential
        potentialName = lastScalarCandidateName;
    }

    if (potentialName.empty())
    {
        std::cerr << "[PlotStandardViewsForPotentialAndEnorm] ERROR: "
                     "could not identify a potential scalar array "
                     "(no V, no Electric_potential, and not exactly one scalar array).\n";
        return;
    }

    std::cout << "[PlotStandardViewsForPotentialAndEnorm] Using scalar '"
              << potentialName << "' as potential.\n";

    // ---------------------------------------------------------------------
    // Compute E and Enorm from the chosen potential (if not already present)
    // ---------------------------------------------------------------------
    bool haveEnorm = (pd->GetArray("Enorm") != nullptr);
    if (!haveEnorm)
    {
        if (!AttachGradientEAndMagnitude(grid, potentialName, "E", "Enorm"))
        {
            std::cerr << "[PlotStandardViewsForPotentialAndEnorm] WARNING: "
                         "failed to attach E/Enorm from '" << potentialName << "'.\n";
        }
        else
        {
            haveEnorm = true;
        }
    }

    // ---------------------------------------------------------------------
    // Load configuration (optional)
    // ---------------------------------------------------------------------
    PlotConfig cfg;
    const std::string cfg_path = GetSourceDir() + "/plot_config.ini";
    bool have_cfg = LoadPlotConfig(cfg_path, cfg);

    // Layout config becomes global
    g_layoutCfg = cfg.layout;

    // Defaults if config does not specify values
    if (cfg.axes.x_label.empty()) cfg.axes.x_label = "x [m]";
    if (cfg.axes.y_label.empty()) cfg.axes.y_label = "y [m]";
    if (cfg.axes.z_label.empty()) cfg.axes.z_label = "z [m]";

    // Text defaults
    if (cfg.text.cbar_num_labels <= 0) cfg.text.cbar_num_labels = 5;

    // Frame defaults (if ini missing, keep legacy behavior)
    if (!have_cfg)
    {
        // full frame
        cfg.frame_full.image_width            = base_options.image_width;
        cfg.frame_full.image_height           = base_options.image_height;
        cfg.frame_full.zoom_factor            = 1.0;
        cfg.frame_full.cbar_horizontal        = false;
        cfg.frame_full.crop_to_region         = false;
        cfg.frame_full.separate_cbar_viewport = true;

        // stack frame (top/bottom)
        cfg.frame_stack = cfg.frame_full;
        cfg.frame_stack.cbar_horizontal = true;
        cfg.frame_stack.crop_to_region  = true;

        // bar frame
        cfg.frame_bar = cfg.frame_full;
        cfg.frame_bar.cbar_horizontal = false;
        cfg.frame_bar.crop_to_region  = true;

        // Regions (legacy constants)
        cfg.regions.top_ymin    = -0.13;
        cfg.regions.top_ymax    =  0.07;
        cfg.regions.bottom_ymin = -1.6;
        cfg.regions.bottom_ymax = -1.4;
        cfg.regions.bar_dx      =  0.1;
        cfg.regions.bar_ymin    = -1.45;
        cfg.regions.bar_ymax    = -0.05;

        // Content defaults for V
        cfg.content_V.scalar_name        = "";          // auto → potentialName
        cfg.content_V.cbar_title         = "V";
        cfg.content_V.color_map_name     = "viridis";
        cfg.content_V.show_contours      = true;
        cfg.content_V.show_streamlines   = false;
        cfg.content_V.n_contours         = 10;
        cfg.content_V.title_full         = "Potential V";
        cfg.content_V.title_stack_top    = "Potential V (top stack)";
        cfg.content_V.title_stack_bottom = "Potential V (bottom stack)";
        cfg.content_V.title_bar          = "Potential V (Field Cage)";

        // Content defaults for E
        cfg.content_E.scalar_name        = "";          // auto → "Enorm"
        cfg.content_E.cbar_title         = "E [V/m]";
        cfg.content_E.color_map_name     = "viridis";
        cfg.content_E.show_contours      = false;
        cfg.content_E.show_streamlines   = true;
        cfg.content_E.n_contours         = 0;
        cfg.content_E.title_full         = "|E|";
        cfg.content_E.title_stack_top    = "|E| (top stack)";
        cfg.content_E.title_stack_bottom = "|E| (bottom stack)";
        cfg.content_E.title_bar          = "|E| (Field Cage)";
    }
    else
    {
        // If image sizes or zoom are still 0, fall back to base_options/1.0
        if (cfg.frame_full.image_width  <= 0) cfg.frame_full.image_width  = base_options.image_width;
        if (cfg.frame_full.image_height <= 0) cfg.frame_full.image_height = base_options.image_height;
        if (cfg.frame_full.zoom_factor  <= 0.0) cfg.frame_full.zoom_factor = 1.0;

        if (cfg.frame_stack.image_width  <= 0) cfg.frame_stack.image_width  = cfg.frame_full.image_width;
        if (cfg.frame_stack.image_height <= 0) cfg.frame_stack.image_height = cfg.frame_full.image_height;
        if (cfg.frame_stack.zoom_factor  <= 0.0) cfg.frame_stack.zoom_factor = cfg.frame_full.zoom_factor;

        if (cfg.frame_bar.image_width  <= 0) cfg.frame_bar.image_width  = cfg.frame_full.image_width;
        if (cfg.frame_bar.image_height <= 0) cfg.frame_bar.image_height = cfg.frame_full.image_height;
        if (cfg.frame_bar.zoom_factor  <= 0.0) cfg.frame_bar.zoom_factor = cfg.frame_full.zoom_factor;
    }

    // ---------------------------------------------------------------------
    // Regions (combine config with grid bounds)
    // ---------------------------------------------------------------------
    double xmin = full_bounds[0];
    double xmax = full_bounds[1];
    double ymin = full_bounds[2];
    double ymax = full_bounds[3];
    double zmin = full_bounds[4];
    double zmax = full_bounds[5];

    double full[6];
    std::copy(std::begin(full_bounds), std::end(full_bounds), full);

    double top[6]    = { xmin, xmax,
                         cfg.regions.top_ymin, cfg.regions.top_ymax,
                         zmin, zmax };

    double bottom[6] = { xmin, xmax,
                         cfg.regions.bottom_ymin, cfg.regions.bottom_ymax,
                         zmin, zmax };

    double right_bar[6] = { xmax - cfg.regions.bar_dx, xmax,
                            cfg.regions.bar_ymin, cfg.regions.bar_ymax,
                            zmin, zmax };

    // ---------------------------------------------------------------------
    // Helper lambda to fill common text options from cfg.text
    // ---------------------------------------------------------------------
    auto apply_text_cfg = [&](ScalarViewRequest& r, const FrameConfig& fc)
    {
        r.title_font_size       = cfg.text.title_font_size;
        r.axis_title_font_size  = cfg.text.axis_title_font_size;
        r.axis_label_font_size  = cfg.text.axis_label_font_size;
        r.cbar_title_font_size  = cfg.text.cbar_title_font_size;
        r.cbar_label_font_size  = cfg.text.cbar_label_font_size;
        r.cbar_minmax_font_size = cfg.text.cbar_minmax_font_size;
        r.cbar_num_labels       = cfg.text.cbar_num_labels;

        if (fc.text.title_font_size       > 0) r.title_font_size       = fc.text.title_font_size;
        if (fc.text.axis_title_font_size  > 0) r.axis_title_font_size  = fc.text.axis_title_font_size;
        if (fc.text.axis_label_font_size  > 0) r.axis_label_font_size  = fc.text.axis_label_font_size;
        if (fc.text.cbar_title_font_size  > 0) r.cbar_title_font_size  = fc.text.cbar_title_font_size;
        if (fc.text.cbar_label_font_size  > 0) r.cbar_label_font_size  = fc.text.cbar_label_font_size;
        if (fc.text.cbar_minmax_font_size > 0) r.cbar_minmax_font_size = fc.text.cbar_minmax_font_size;
        if (fc.text.cbar_num_labels       > 0) r.cbar_num_labels       = fc.text.cbar_num_labels;
    };

    // ---------------------------------------------------------------------
    // Decide which scalar name to use for V plots:
    //  - If cfg.content_V.scalar_name exists on the dataset, use it.
    //  - Otherwise, always fall back to the detected potentialName.
    // ---------------------------------------------------------------------
    const ContentConfig& cv = cfg.content_V;

    std::string vScalarName;
    if (!cv.scalar_name.empty() && pd->GetArray(cv.scalar_name.c_str()))
    {
        vScalarName = cv.scalar_name;
    }
    else
    {
        if (!cv.scalar_name.empty() && !pd->GetArray(cv.scalar_name.c_str()))
        {
            std::cerr << "[PlotStandardViewsForPotentialAndEnorm] WARNING: "
                         "content_V.scalar_name='" << cv.scalar_name
                      << "' not found on dataset; using detected potential '"
                      << potentialName << "' instead.\n";
        }
        vScalarName = potentialName;
    }

    // ---------------------------------------------------------------------
    // V views
    // ---------------------------------------------------------------------

    // Full-domain V
    ScalarViewRequest r_full;
    r_full.scalar_name            = vScalarName;
    r_full.cbar_title             = cv.cbar_title.empty() ? "V" : cv.cbar_title;
    r_full.cbar_horizontal        = cfg.frame_full.cbar_horizontal;
    r_full.title                  = cv.title_full.empty() ? "Potential V" : cv.title_full;
    r_full.x_label                = cfg.axes.x_label;
    r_full.y_label                = cfg.axes.y_label;
    r_full.z_label                = cfg.axes.z_label;
    r_full.color_map_name         = cv.color_map_name.empty() ? "viridis" : cv.color_map_name;
    std::copy(full, full + 6, r_full.region_bounds);
    r_full.crop_to_region         = cfg.frame_full.crop_to_region;
    r_full.separate_cbar_viewport = cfg.frame_full.separate_cbar_viewport;
    r_full.show_contours          = cv.show_contours;
    r_full.show_streamlines       = cv.show_streamlines;
    r_full.n_contours             = cv.n_contours;
    r_full.output_path            = input.out_dir + "/V_full.png";
    r_full.image_width            = cfg.frame_full.image_width;
    r_full.image_height           = cfg.frame_full.image_height;
    r_full.zoom_factor            = cfg.frame_full.zoom_factor;
    apply_text_cfg(r_full, cfg.frame_full);
    r_full.x_num_labels = cfg.frame_full.x_num_labels;

    // Top V (stack)
    ScalarViewRequest r_top = r_full;
    r_top.cbar_horizontal        = cfg.frame_stack.cbar_horizontal;
    r_top.title                  = cv.title_stack_top.empty()
                                   ? "Potential V (top stack)"
                                   : cv.title_stack_top;
    std::copy(top, top + 6, r_top.region_bounds);
    r_top.crop_to_region         = cfg.frame_stack.crop_to_region;
    r_top.separate_cbar_viewport = cfg.frame_stack.separate_cbar_viewport;
    r_top.output_path            = input.out_dir + "/V_top.png";
    r_top.image_width            = cfg.frame_stack.image_width;
    r_top.image_height           = cfg.frame_stack.image_height;
    r_top.zoom_factor            = cfg.frame_stack.zoom_factor;
    apply_text_cfg(r_top, cfg.frame_stack);
    r_top.x_num_labels = cfg.frame_stack.x_num_labels;

    // Bottom V (stack)
    ScalarViewRequest r_bottom = r_full;
    r_bottom.cbar_horizontal        = cfg.frame_stack.cbar_horizontal;
    r_bottom.title                  = cv.title_stack_bottom.empty()
                                      ? "Potential V (bottom stack)"
                                      : cv.title_stack_bottom;
    std::copy(bottom, bottom + 6, r_bottom.region_bounds);
    r_bottom.crop_to_region         = cfg.frame_stack.crop_to_region;
    r_bottom.separate_cbar_viewport = cfg.frame_stack.separate_cbar_viewport;
    r_bottom.output_path            = input.out_dir + "/V_bottom.png";
    r_bottom.image_width            = cfg.frame_stack.image_width;
    r_bottom.image_height           = cfg.frame_stack.image_height;
    r_bottom.zoom_factor            = cfg.frame_stack.zoom_factor;
    apply_text_cfg(r_bottom, cfg.frame_stack);
    r_bottom.x_num_labels = cfg.frame_stack.x_num_labels;

    // Bar V (field cage)
    ScalarViewRequest r_bar = r_full;
    r_bar.cbar_horizontal        = cfg.frame_bar.cbar_horizontal;
    r_bar.title                  = cv.title_bar.empty()
                                   ? "Potential V (Field Cage)"
                                   : cv.title_bar;
    std::copy(right_bar, right_bar + 6, r_bar.region_bounds);
    r_bar.crop_to_region         = cfg.frame_bar.crop_to_region;
    r_bar.separate_cbar_viewport = cfg.frame_bar.separate_cbar_viewport;
    r_bar.output_path            = input.out_dir + "/V_bar.png";
    r_bar.image_width            = cfg.frame_bar.image_width;
    r_bar.image_height           = cfg.frame_bar.image_height;
    r_bar.zoom_factor            = cfg.frame_bar.zoom_factor;
    apply_text_cfg(r_bar, cfg.frame_bar);
    r_bar.x_num_labels = cfg.frame_bar.x_num_labels;

    PlotScalarFieldView(grid, r_full,   base_options);
    PlotScalarFieldView(grid, r_top,    base_options);
    PlotScalarFieldView(grid, r_bottom, base_options);
    PlotScalarFieldView(grid, r_bar,    base_options);

    // ---------------------------------------------------------------------
    // |E| views – only if Enorm is available
    // ---------------------------------------------------------------------
    if (haveEnorm)
    {
        const ContentConfig& ce = cfg.content_E;
        StreamlineConfig stream_cfg; // placeholder

        std::string eScalarName = !ce.scalar_name.empty()
                                  ? ce.scalar_name
                                  : std::string("Enorm");

        ScalarViewRequest r_fullE = r_full;
        r_fullE.scalar_name       = eScalarName;
        r_fullE.cbar_title        = ce.cbar_title.empty() ? "E [V/m]" : ce.cbar_title;
        r_fullE.title             = ce.title_full.empty() ? "|E|"     : ce.title_full;
        r_fullE.color_map_name    = ce.color_map_name.empty()
                                    ? r_full.color_map_name
                                    : ce.color_map_name;
        r_fullE.show_contours     = ce.show_contours;
        r_fullE.show_streamlines  = ce.show_streamlines;
        r_fullE.n_contours        = ce.n_contours;
        r_fullE.output_path       = input.out_dir + "/Enorm_full.png";
        apply_text_cfg(r_fullE, cfg.frame_full);

        ScalarViewRequest r_topE = r_top;
        r_topE.scalar_name       = eScalarName;
        r_topE.cbar_title        = r_fullE.cbar_title;
        r_topE.title             = ce.title_stack_top.empty()
                                   ? "|E| (top stack)"
                                   : ce.title_stack_top;
        r_topE.color_map_name    = r_fullE.color_map_name;
        r_topE.show_contours     = ce.show_contours;
        r_topE.show_streamlines  = ce.show_streamlines;
        r_topE.n_contours        = ce.n_contours;
        r_topE.output_path       = input.out_dir + "/Enorm_top.png";
        apply_text_cfg(r_topE, cfg.frame_stack);

        ScalarViewRequest r_bottomE = r_bottom;
        r_bottomE.scalar_name       = eScalarName;
        r_bottomE.cbar_title        = r_fullE.cbar_title;
        r_bottomE.title             = ce.title_stack_bottom.empty()
                                      ? "|E| (bottom stack)"
                                      : ce.title_stack_bottom;
        r_bottomE.color_map_name    = r_fullE.color_map_name;
        r_bottomE.show_contours     = ce.show_contours;
        r_bottomE.show_streamlines  = ce.show_streamlines;
        r_bottomE.n_contours        = ce.n_contours;
        r_bottomE.output_path       = input.out_dir + "/Enorm_bottom.png";
        apply_text_cfg(r_bottomE, cfg.frame_stack);

        ScalarViewRequest r_barE = r_bar;
        r_barE.scalar_name       = eScalarName;
        r_barE.cbar_title        = r_fullE.cbar_title;
        r_barE.title             = ce.title_bar.empty()
                                   ? "|E| (Field Cage)"
                                   : ce.title_bar;
        r_barE.color_map_name    = r_fullE.color_map_name;
        r_barE.show_contours     = ce.show_contours;
        r_barE.show_streamlines  = ce.show_streamlines;
        r_barE.n_contours        = ce.n_contours;
        r_barE.output_path       = input.out_dir + "/Enorm_bar.png";
        apply_text_cfg(r_barE, cfg.frame_bar);

        PlotScalarFieldView(grid, r_fullE,   base_options, &stream_cfg);
        PlotScalarFieldView(grid, r_topE,    base_options, &stream_cfg);
        PlotScalarFieldView(grid, r_bottomE, base_options, &stream_cfg);
        PlotScalarFieldView(grid, r_barE,    base_options, &stream_cfg);

        std::cout << "Standard V / |E| views written to: " << input.out_dir << "\n";
    }
    else
    {
        std::cerr << "[PlotStandardViewsForPotentialAndEnorm] NOTE: Enorm not available; "
                     "only V views were written.\n";
        std::cout << "Standard V views written to: " << input.out_dir << "\n";
    }
}



// -----------------------------------------------------------------------------
// Input preparation
// -----------------------------------------------------------------------------

PlotInput PreparePlotInput(const char* raw_path)
{
    PlotInput R;

    std::string path = raw_path ? std::string(raw_path) : std::string();

    if (!path.empty() && path[0] == '~')
    {
        const char* home = std::getenv("HOME");
        if (!home)
        {
            std::cerr << "[PreparePlotInput] ERROR: HOME not set.\n";
            std::exit(1);
        }
        path = std::string(home) + path.substr(1);
    }

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

    if (ends_with(path, ".pvd"))
    {
        R.base_dir = dirname_of(path);
        R.vtu_file = R.base_dir + "/Cycle000000/proc000000.vtu";
    }
    else
    {
        R.vtu_file = path;
        R.base_dir = dirname_of(path);
    }

    std::cout << "[PreparePlotInput] Reading VTU: " << R.vtu_file << "\n";

    R.grid = LoadGridFromVTU(R.vtu_file);
    if (!R.grid)
    {
        std::cerr << "[PreparePlotInput] ERROR: Failed to load grid from "
                  << R.vtu_file << "\n";
        std::exit(1);
    }

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

} // namespace FEMPlot

// -----------------------------------------------------------------------------
// Example standalone main entry
// -----------------------------------------------------------------------------

int make_plots(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << (argc > 0 ? argv[0] : "make_plots")
                  << " <path-to.vtu|.pvd>\n";
        return 1;
    }

    auto input = FEMPlot::PreparePlotInput(argv[1]);

    FEMPlot::PlottingOptions opts;
    FEMPlot::PlotStandardViewsForPotentialAndEnorm(input.grid, input, opts);

    return 0;
}
