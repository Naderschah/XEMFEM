#include "FEMPlotVTK.h"

int main(int argc, char** argv) {
    make_plots(argc, argv);  
    return 0;
}
/*

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
// Mesh IO
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>

// Scalars
#include <vtkPointData.h>
#include <vtkDataArray.h>


int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <mesh.vtu>\n";
    return EXIT_FAILURE;
  }

  const char* meshFile = argv[1];

  vtkNew<vtkNamedColors> colors;

  // Setup offscreen rendering (same as example)
  vtkNew<vtkGraphicsFactory> graphics_factory;
  graphics_factory->SetOffScreenOnlyMode(1);
  graphics_factory->SetUseMesaClasses(1);

  // ------------------------------------------------------------------
  // Read unstructured mesh from VTU
  // ------------------------------------------------------------------
  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(meshFile);
  reader->Update();

  vtkUnstructuredGrid* grid = reader->GetOutput();
  if (!grid)
  {
    std::cerr << "Failed to read mesh from " << meshFile << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "Grid: " << grid->GetNumberOfPoints() << " points, "
            << grid->GetNumberOfCells() << " cells\n";

  // ------------------------------------------------------------------
  // Create a mapper and actor for the mesh
  // ------------------------------------------------------------------
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputData(grid);

  // ------------------------------------------------------------------
  // Configure mapper to show scalar field "V" (voltage) on points
  // ------------------------------------------------------------------
  const char* scalarName = "V";

  vtkPointData* pd = grid->GetPointData();
  vtkDataArray* V = pd ? pd->GetArray(scalarName) : nullptr;
  if (!V)
  {
    std::cerr << "Warning: point-data array \"" << scalarName
              << "\" not found. Rendering geometry only.\n";
    mapper->ScalarVisibilityOff();
  }
  else
  {
    double range[2];
    V->GetRange(range);

    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(scalarName);
    mapper->ScalarVisibilityOn();
    mapper->SetScalarRange(range);

    std::cout << "Using scalar \"" << scalarName << "\" with range ["
              << range[0] << ", " << range[1] << "]\n";
  }

  // ------------------------------------------------------------------
  // Actor for the mesh
  // ------------------------------------------------------------------
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  // If no scalars, use a solid color
  if (!V)
  {
    actor->GetProperty()->SetColor(colors->GetColor3d("White").GetData());
  }

  // ------------------------------------------------------------------
  // A renderer and render window (offscreen)
  // ------------------------------------------------------------------
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->SetOffScreenRendering(1);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(800, 600);

  // Add the mesh actor to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

  // Frame the mesh
  double bounds[6];
  grid->GetBounds(bounds);
  renderer->ResetCamera(bounds);

  renderWindow->Render();

  // ------------------------------------------------------------------
  // Capture to PNG via WindowToImageFilter
  // ------------------------------------------------------------------
  vtkNew<vtkWindowToImageFilter> windowToImageFilter;
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetInputBufferTypeToRGBA();
  windowToImageFilter->ReadFrontBufferOff();
  windowToImageFilter->Update();

  vtkNew<vtkPNGWriter> writer;
  writer->SetFileName("mesh_voltage.png");
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();

  return EXIT_SUCCESS;
}
*/
