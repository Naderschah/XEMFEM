// FEMPlotVTK.h
//
// High-level API for VTK-based FEM plotting.

#pragma once

#include <array>
#include <string>
#include <vector>

#include <vtkSmartPointer.h>

// Forward declarations for VTK types
class vtkDataSet;
class vtkRenderer;
class vtkRenderWindow;
class vtkUnstructuredGrid;

namespace FEMPlot
{
// Per-view configuration loaded from a text file
struct ViewConfig
{
    int    image_width  = 0;    // 0 -> keep compiled default
    int    image_height = 0;    // 0 -> keep compiled default
    double zoom_factor  = 0.0;  // <=0 -> keep compiled default
};
// -----------------------------------------------------------------------------
// Basic data structures
// -----------------------------------------------------------------------------
struct ScalarStats
{
    double min_roi  = 0.0;  // true min in the dataset / ROI
    double max_roi  = 0.0;  // true max
    double min_used = 0.0;  // range actually used for coloring
    double max_used = 0.0;
};

// Result container for vector fields at mesh points (e.g. E-field).
struct Result
{
    // One 3D vector per mesh point (indexed like VTK point ids).
    std::vector<std::array<double, 3>> E_at_points;
};

// Global plotting configuration shared by views.
struct PlottingOptions
{
    // Image resolution for offscreen rendering
    int image_width  = 1600;
    int image_height = 1000;

    // Mesh rendering
    bool  show_edges  = false;
    double edge_width = 1.0;

    // Legacy / optional scalar-bar geometry hints (kept for compatibility)
    double cbar_bar_ratio       = 0.5;
    double cbar_width           = 0.08;
    double cbar_height          = 0.6;
    double cbar_pos_x           = 0.90;
    double cbar_pos_y           = 0.15;
    double cbar_TitleRatio      = 0.3;
    double cbar_TitleSeperation = 0.5;
};

// Streamline configuration for E-field visualization.
struct StreamlineConfig
{
    int    n_seeds           = 500;   // number of seed points
    double max_propagation   = 1.0;   // integration length
    double initial_step      = 0.01;
    double min_step          = 1e-4;
    double max_step          = 0.1;
    int    max_steps         = 2000;
    double tube_radius_rel   = 0.005; // tube radius relative to domain size
};

// Input description for a plotting run.
struct PlotInput
{
    // Directory that contains the VTU or PVD.
    std::string base_dir;

    // The VTU file actually read by VTK.
    std::string vtu_file;

    // Loaded grid (unstructured mesh).
    vtkSmartPointer<vtkUnstructuredGrid> grid;

    // Output directory where PNGs are written.
    std::string out_dir;
};

// Description of a single scalar-field view to be rendered to PNG.
struct ScalarViewRequest
{
    // Which scalar to visualize (must exist as point-data array on the dataset)
    std::string scalar_name;       // e.g. "V", "Enorm"

    std::string cbar_title;
    bool   cbar_horizontal = false;

    // Human-readable title for the plot (shown at top)
    std::string title;

    // Axis labels (including units)
    std::string x_label;           // e.g. "x [m]"
    std::string y_label;           // e.g. "y [m]"
    std::string z_label;           // e.g. "z [m]"

    int x_num_labels = 0; 
    // Physical region to frame the camera (and optionally clip to)
    double region_bounds[6] = {0,0,0,0,0,0};

    // Name of the color map / LUT preset (e.g. "viridis", "gray", etc.)
    std::string color_map_name;

    // If true: clip visually to region_bounds using mapper clipping planes.
    // If false: show full geometry, but still use region_bounds for camera/framing.
    bool crop_to_region = true;

    // If true: place colorbar in its own viewport (no overlap with field view).
    bool separate_cbar_viewport = true;

    // Overlays:
    bool show_contours    = false;   // for V
    bool show_streamlines = false;   // for |E|

    // Number of contour levels if show_contours = true.
    int n_contours = 0;

    // Output file
    std::string output_path;

    // -------------------------------------------------------------------------
    // Per-view layout/appearance controls
    // -------------------------------------------------------------------------

    // Image size override for this view (pixels).
    // If <= 0, fall back to PlottingOptions.image_width/height.
    int image_width  = 0;
    int image_height = 0;

    // Camera zoom factor (1.0 = default; <1 => "zoom in", >1 => "zoom out").
    double zoom_factor = 1.0;

    // Font sizes (pixels). If <= 0, automatic scaling is used based on
    // the viewport size.
    int title_font_size       = 22;
    int axis_title_font_size  = 18;
    int axis_label_font_size  = 14;
    int cbar_title_font_size  = 18;
    int cbar_label_font_size  = 14;
    int cbar_minmax_font_size = 12;

    // Colorbar tick labels
    int cbar_num_labels       = 5;
};
// -----------------------------------------------------------------------------
// INI-based configuration
// -----------------------------------------------------------------------------
struct LayoutConfig
{
    // Horizontal colorbar layout (cbarHorizontal == true)
    double h_label_y0  = 0.00;  // min/max label band bottom
    double h_label_y1  = 0.06;  // min/max label band top
    double h_cbar_y0   = 0.06;  // cbar band bottom
    double h_cbar_y1   = 0.16;  // cbar band top
    double h_ctitle_y0 = 0.16;  // cbar title band bottom
    double h_ctitle_y1 = 0.22;  // cbar title band top
    double h_field_y0  = 0.22;  // field band bottom
    double h_field_y1  = 0.94;  // field band top
    double h_title_y0  = 0.94;  // plot title band bottom
    double h_title_y1  = 1.00;  // plot title band top

    double h_cbar_x0   = 0.10;  // horizontal cbar left
    double h_cbar_x1   = 0.90;  // horizontal cbar right

    // Vertical colorbar layout (cbarHorizontal == false)
    double v_field_x0      = 0.00; // field left
    double v_field_x1      = 0.75; // field right

    double v_label_bot_y0  = 0.00; // bottom label band bottom
    double v_label_bot_y1  = 0.06; // bottom label band top
    double v_cbar_y0       = 0.06; // cbar band bottom
    double v_cbar_y1       = 0.82; // cbar band top
    double v_label_top_y0  = 0.82; // top label band bottom
    double v_label_top_y1  = 0.88; // top label band top
    double v_title_y0      = 0.88; // title band bottom
    double v_title_y1      = 1.00; // title band top

    double v_cbar_x0       = 0.79; // vertical cbar left
    double v_cbar_x1       = 0.97; // vertical cbar right
};
struct TextConfig
{
    int title_font_size       = 0;
    int axis_title_font_size  = 0;
    int axis_label_font_size  = 0;
    int cbar_title_font_size  = 0;
    int cbar_label_font_size  = 0;
    int cbar_minmax_font_size = 0;
    int cbar_num_labels       = 0;
};
struct FrameConfig
{
    // View-geometry / layout: shared by V and E for the same frame
    int    image_width  = 0;    // 0 → use compiled default
    int    image_height = 0;    // 0 → use compiled default
    double zoom_factor  = 0.0;  // <=0 → use compiled default

    bool   cbar_horizontal        = false;
    bool   crop_to_region         = true;
    bool   separate_cbar_viewport = true;

    int x_num_labels = 0; 

    // NEW: optional per-frame font overrides (0 → use global / auto)
    TextConfig text;
};

struct ContentConfig
{
    // What is being plotted (V vs E)
    std::string scalar_name;      // e.g. "V" or "Enorm"
    std::string cbar_title;       // e.g. "V" or "E [V/m]"
    std::string color_map_name;   // e.g. "viridis"

    // Per-content display options
    bool show_contours    = false;
    bool show_streamlines = false;
    int  n_contours       = 0;

    // Titles per frame
    std::string title_full;
    std::string title_stack_top;
    std::string title_stack_bottom;
    std::string title_bar;
};

struct AxesConfig
{
    std::string x_label;  // e.g. "x [m]"
    std::string y_label;  // e.g. "y [m]"
    std::string z_label;  // e.g. "z [m]"
};


struct RegionConfig
{
    // Y-bounds for “stack” regions (absolute values in same units as mesh)
    double top_ymin    = 0.0;
    double top_ymax    = 0.0;
    double bottom_ymin = 0.0;
    double bottom_ymax = 0.0;

    // For bar region: width in x and y-extent
    double bar_dx      = 0.0; // distance from xmax: [xmax - bar_dx, xmax]
    double bar_ymin    = 0.0;
    double bar_ymax    = 0.0;
};



struct PlotConfig
{
    // Frame configs (shared by V/E)
    FrameConfig frame_full;
    FrameConfig frame_stack;
    FrameConfig frame_bar;

    // Content-specific configs
    ContentConfig content_V;
    ContentConfig content_E;

    // Shared across everything
    AxesConfig   axes;
    TextConfig   text;
    RegionConfig regions;

    // NEW: viewport layout configuration
    LayoutConfig layout;
};
// Load config from an INI file; returns true on success.
// Unspecified fields remain at their default values.
bool LoadPlotConfig(const std::string& path, PlotConfig& cfg);

// -----------------------------------------------------------------------------
// Viewport layout for a single scalar-field plot
// -----------------------------------------------------------------------------
struct ViewportRect
{
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 1.0;
    double y1 = 1.0;
};

struct PlotViewLayout
{
    // Region occupied by this plot in the render window (normalized).
    ViewportRect plot;

    // Main components inside that plot:
    ViewportRect field;              // main scalar field view
    ViewportRect cbar;               // colorbar view
    ViewportRect title;              // plot title area
    ViewportRect cbarTitle;          // colorbar title area
    ViewportRect cbarLabelTopOrLeft;
    ViewportRect cbarLabelBottomOrRight;

    bool cbarHorizontal = true;
};

// Compute a default layout (normalized viewports) for a single plot,
// given the requested colorbar orientation and whether it should have
// its own dedicated viewport bands.
PlotViewLayout MakeDefaultPlotViewLayout(bool cbarHorizontal,
                                         bool separateCbarViewport);

// -----------------------------------------------------------------------------
// Core IO / utilities
// -----------------------------------------------------------------------------

// Load an unstructured grid from a VTU file.
// Returns nullptr on failure.
vtkSmartPointer<vtkUnstructuredGrid>
LoadGridFromVTU(const std::string& filename);

// Attach a vector field E = -∇(potential_array_name) and its magnitude |E| to
// an unstructured grid as point-data arrays.
//
// On success:
//   - a 3-component vector array named e_array_name is added,
//   - a 1-component scalar array named emag_array_name is added.
//
// Returns false on error (missing potential, gradient failure, etc.).
bool AttachGradientEAndMagnitude(vtkUnstructuredGrid* grid,
                                 const std::string& potential_array_name,
                                 const std::string& e_array_name,
                                 const std::string& emag_array_name);

// -----------------------------------------------------------------------------
// Plotting primitives
// -----------------------------------------------------------------------------

// Add streamlines of the E-field to a renderer. The E-field at VTK points is
// provided in `result` and will be attached/updated as a vector array "E" on
// the given dataset.
void AddStreamlines(vtkDataSet* dataset,
                    vtkRenderer* renderer,
                    const Result& result,
                    const StreamlineConfig& cfg);

// Offscreen rendering helper for a single renderer.
// Creates an offscreen render window of (width x height), attaches the
// renderer, captures it as PNG, and writes to `filename`.
void RenderOffscreenPNG(vtkRenderer* renderer,
                        const std::string& filename,
                        int width,
                        int height);

// -----------------------------------------------------------------------------
// High-level view plotting
// -----------------------------------------------------------------------------

// Render a single scalar-field view (V or |E|) according to `request` and
// `plotting`, optionally adding streamlines configured by `stream_cfg`
// (used only if request.show_streamlines == true).
//
// This function:
//   - applies optional ROI cropping (crop_to_region),
//   - computes local scalar statistics for the region,
//   - chooses a display range (with internal clipping of outliers),
//   - builds the VTK pipeline (mapper, actors, axes, colorbar, contours,
//     streamlines as requested),
//   - sets up camera and view layout (including colorbar viewport),
//   - performs offscreen rendering,
//   - writes the PNG file to request.output_path.
//
// It does not modify the input grid except for reading its arrays.
void PlotScalarFieldView(vtkUnstructuredGrid* grid,
                         const ScalarViewRequest& request,
                         const PlottingOptions& plotting,
                         const StreamlineConfig* stream_cfg = nullptr);

// Convenience driver that produces the standard set of plots:
//   - V (full, top, bottom, bar region)
//   - |E| for the same views
// using the given `grid`, `input` (paths) and `base_options` (resolution, etc).
void PlotStandardViewsForPotentialAndEnorm(vtkUnstructuredGrid* grid,
                                           const PlotInput& input,
                                           const PlottingOptions& base_options);

// Example driver used by `make_plots`
PlotInput PreparePlotInput(const char* raw_path);

} // namespace FEMPlot

// Main entry point for a standalone tool (optional)
int make_plots(int argc, char** argv);

