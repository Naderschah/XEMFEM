/*
Contains the main plotting functions 
*/

#ifndef FEM_PLOT_VTK_H
#define FEM_PLOT_VTK_H

#include <array>
#include <string>
#include <utility>
#include <vector>
#include <limits>
#include <vtkSmartPointer.h>

// Forward declarations of VTK classes used in the public API.
// Full definitions are included in the .cpp.
class vtkUnstructuredGrid;
class vtkRenderer;
class vtkActor;
class vtkDataSetMapper;

namespace FEMPlot
{
// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

// E field sampled at VTK points, in VTK point order.
struct Result
{
    std::vector<std::array<double, 3>> E_at_points;
};

// Parameters controlling streamline generation and appearance.
struct StreamlineConfig
{
    int    n_seeds           = 200;
    double max_propagation   = 1000.0;
    double initial_step      = 0.01;
    double min_step          = 0.001;
    double max_step          = 1.0;
    int    max_steps         = 10000;
    double tube_radius_rel   = 0.002; // relative to domain radius
};

// ---------------------------------------------------------------------------
// Core helpers
// ---------------------------------------------------------------------------

// Read an unstructured grid from a VTU file.
vtkSmartPointer<vtkUnstructuredGrid>
LoadGridFromVTU(const std::string& filename);

// Create a mapper for a named point-data scalar on the grid.
// If range_min and range_max are NaN, the global scalar range is used.
vtkSmartPointer<vtkDataSetMapper>
CreateScalarMapper(vtkUnstructuredGrid* grid,
                   const std::string& scalar_name,
                   double range_min = std::numeric_limits<double>::quiet_NaN(),
                   double range_max = std::numeric_limits<double>::quiet_NaN());

// Create a mesh actor from a mapper, optionally showing edges.
vtkSmartPointer<vtkActor>
CreateMeshActor(vtkSmartPointer<vtkDataSetMapper> mapper,
                bool show_edges   = true,
                double line_width = 1.0);

// Create a renderer with a single actor and background color.
// If background is nullptr, white (1,1,1) is used.
vtkSmartPointer<vtkRenderer>
CreateRendererWithActor(vtkSmartPointer<vtkActor> actor,
                        const double* background = nullptr);

// Compute scalar range for a named point-data array within region_bounds.
// region_bounds = {xmin,xmax, ymin,ymax, zmin,zmax}.
std::pair<double, double>
ComputeLocalScalarRange(vtkUnstructuredGrid* grid,
                        const std::string& scalar_name,
                        const double region_bounds[6]);

// ---------------------------------------------------------------------------
// High-level plot builders
// ---------------------------------------------------------------------------

// Full-domain potential renderer with scalar bar (global range).
vtkSmartPointer<vtkRenderer>
CreatePotentialRenderer(vtkUnstructuredGrid* grid,
                        const std::string& scalar_name);

// Zoomed potential renderer with local scalar range.
// region_bounds = {xmin,xmax, ymin,ymax, zmin,zmax}.
vtkSmartPointer<vtkRenderer>
CreateZoomedRenderer(vtkUnstructuredGrid* grid,
                     const std::string& scalar_name,
                     const double region_bounds[6]);

// ---------------------------------------------------------------------------
// Decorations and overlays
// ---------------------------------------------------------------------------

// Add x/y(/z) axes for given bounds (typically full or zoom bounds).
void AddXYAxes(vtkRenderer* renderer,
               const double bounds[6],
               const std::string& x_label = "x",
               const std::string& y_label = "y",
               const std::string& z_label = "z");

// Ensure a scalar bar is present; attach to the first mapper with a LUT.
// If title is non-empty, use it as the scalar bar title.
void EnsureScalarBar(vtkRenderer* renderer,
                     const std::string& title = "");

// Add streamlines of E(x) on top of an existing renderer.
// E field is taken from Result::E_at_points in VTK point order.
void AddStreamlines(vtkUnstructuredGrid* grid,
                    vtkRenderer* renderer,
                    const Result& result,
                    const StreamlineConfig& cfg = StreamlineConfig{});

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

// Render a single renderer offscreen into a PNG file.
void RenderOffscreenPNG(vtkRenderer* renderer,
                        const std::string& filename,
                        int width  = 1920,
                        int height = 1080);

// Render three renderers side by side into a single PNG
// (3 columns, 1 row).
void RenderThreeSideBySide(vtkRenderer* left,
                           vtkRenderer* middle,
                           vtkRenderer* right,
                           const std::string& filename,
                           int width  = 2400,
                           int height = 800);

} // namespace FEMPlot

int make_plots(int, char**);



#endif // FEM_PLOT_VTK_H