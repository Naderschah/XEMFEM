/*
Contains the main plotting functions 
*/

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

// -----------------------------------------------------------------------------
// Basic data structures
// -----------------------------------------------------------------------------
struct ScalarStats
{
    double min_roi     = 0.0;  // true min in the dataset / ROI
    double max_roi     = 0.0;  // true max
    double min_used    = 0.0;  // range actually used for coloring
    double max_used    = 0.0;
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

    // Legacy / optional scalar-bar geometry hints (used by older paths).
    double cbar_bar_ratio      = 0.5;
    double cbar_width          = 0.08;
    double cbar_height         = 0.6;
    double cbar_pos_x          = 0.90;
    double cbar_pos_y          = 0.15;
    double cbar_TitleRatio     = 0.3;
    double cbar_TitleSeperation = 0.5;

    // Optional name of the scalar to visualize (used in some older helpers).
    std::string scalar_name;
};

// Streamline configuration for E-field visualization.
struct StreamlineConfig
{
    int    n_seeds           = 500;   // number of seed points
    double max_propagation   = 1.0;   // integration length in "space" units
    double initial_step      = 0.01;
    double min_step          = 1e-4;
    double max_step          = 0.1;
    int    max_steps         = 2000;
    double tube_radius_rel   = 0.005; // tube radius relative to domain size
};

// Input description for a plotting run (path handling, grid, output dir).
struct PlotInput
{
    // Base directory of the simulation (where .pvd/.vtu come from).
    std::string base_dir;

    // Concrete VTU file used.
    std::string vtu_file;

    // Mesh grid loaded from VTU.
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
    std::string title;             // e.g. "Electric potential V [V]"

    // Axis labels (including units)
    std::string x_label;           // e.g. "x [m]"
    std::string y_label;           // e.g. "y [m]"
    std::string z_label;           // e.g. "z [m]" or empty if unused

    // Name of the color map / palette to use (VTK preset or custom key)
    std::string color_map_name;    // e.g. "viridis", "coolwarm", ...

    // Region of interest in model coordinates:
    // [xmin, xmax, ymin, ymax, zmin, zmax]
    double region_bounds[6] = {0, 0, 0, 0, 0, 0};

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

    // ---------------------------------------------------------------------
    // New per-view layout/appearance controls
    // ---------------------------------------------------------------------

    // Image size override for this view (pixels).
    // If <= 0, fall back to PlottingOptions.image_width/height.
    int image_width  = 0;
    int image_height = 0;

    // Camera zoom factor (1.0 = default; <1 => "zoom in", >1 => "zoom out").
    double zoom_factor = 1.0;

    // Font sizes (pixels)
    int title_font_size         = 22;
    int axis_title_font_size    = 18;
    int axis_label_font_size    = 14;
    int cbar_title_font_size    = 18;
    int cbar_label_font_size    = 14;
    int cbar_minmax_font_size   = 12;

    // Colorbar tick labels
    int cbar_num_labels         = 5;
};

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
//
//   - For V:
//       * full domain
//       * top zoom
//       * bottom zoom
//       * right-side vertical bar
//
//   - For |E| (Enorm):
//       * the same four views, with streamlines (as support is added).
//
// It assumes that:
//   - the potential array "V" already exists on `grid`,
//   - AttachGradientEAndMagnitude(...) is called inside this
//     function to add "E" and "Enorm".
//
// Output filenames and view regions are constructed from `input` and
// `base_options` (using input.out_dir and grid bounds).
void PlotStandardViewsForPotentialAndEnorm(vtkUnstructuredGrid* grid,
                                           const PlotInput& input,
                                           const PlottingOptions& base_options);

// -----------------------------------------------------------------------------
// Example driver (used by your existing main)
// -----------------------------------------------------------------------------


} // namespace


// Main entry point
int make_plots(int argc, char** argv);