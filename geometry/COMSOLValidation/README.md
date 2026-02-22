This folder holds configuration files and matched COMSOL files to validate the function of XEMFEM.

This is required as matching the more complex TPC geometry has prooved prohibitatively complex.

The examples include in 2D only (quicker to implement in both frameworks)

1. Parallel plate capacitor in volume – verifies uniform-field accuracy, fringing fields, and sensitivity to outer boundary conditions.

2. Eccentric coaxial cylinder – tests asymmetric gap handling, surface charge redistribution, and peak-field localization due to misalignment.

3. Two opposite-polarity wires – probes small-gap field amplification and resolution of near-contact electrostatic interactions.

4. Conducting wedge facing a plane – evaluates solver behavior near true corner singularities and field-scaling with distance from the edge.

To mesh these examples:

```bash 
./make_mesh.sh tui --config_path COMSOLValidation/TwoWires/
```

Then proceed to run XEMFEM pointing the config to the relevant ones


In paraview then to compare, this script can be used (View->python shell):

```python
study_to_load_from = 'TwoWires'

# trace generated using paraview version 5.13.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
cOMSOLvtu = XMLUnstructuredGridReader(registrationName='COMSOL.vtu', FileName=['/home/felix/MFEMElectrostatics/sim_results/{}/COMSOL.vtu'.format(study_to_load_from)])

# create a new 'PVD Reader'
simulationpvd = PVDReader(registrationName='Simulation.pvd', FileName='/home/felix/MFEMElectrostatics/sim_results/{}/Simulation/Simulation.pvd'.format(study_to_load_from))

# Properties modified on cOMSOLvtu
cOMSOLvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
cOMSOLvtuDisplay = Show(cOMSOLvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cOMSOLvtuDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# show data in view
simulationpvdDisplay = Show(simulationpvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
simulationpvdDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(cOMSOLvtu)

# set active source
SetActiveSource(simulationpvd)

# set active source
SetActiveSource(cOMSOLvtu)

# set active source
SetActiveSource(simulationpvd)

# create a new 'Resample With Dataset'
resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=simulationpvd,
    DestinationMesh=cOMSOLvtu)

RenameProxy(resampleWithDataset1, 'sources', 'ResampleXEMFEM')

# rename source object
RenameSource('ResampleXEMFEM', resampleWithDataset1)

# show data in view
resampleWithDataset1Display = Show(resampleWithDataset1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
resampleWithDataset1Display.Representation = 'Surface'

# hide data in view
Hide(cOMSOLvtu, renderView1)

# hide data in view
Hide(simulationpvd, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Resample With Dataset'
resampleWithDataset1_1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=resampleWithDataset1,
    DestinationMesh=simulationpvd)

RenameProxy(resampleWithDataset1_1, 'sources', 'SampleToSource')

# rename source object
RenameSource('SampleToSource', resampleWithDataset1_1)

# show data in view
resampleWithDataset1_1Display = Show(resampleWithDataset1_1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
resampleWithDataset1_1Display.Representation = 'Surface'

# hide data in view
Hide(simulationpvd, renderView1)

# hide data in view
Hide(resampleWithDataset1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(simulationpvd)

# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1', Input=[resampleWithDataset1_1, simulationpvd])

RenameProxy(appendAttributes1, 'sources', 'XEMFEM_And_Resampled')

# rename source object
RenameSource('XEMFEM_And_Resampled', appendAttributes1)

# show data in view
appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1Display.Representation = 'Surface'

# hide data in view
Hide(simulationpvd, renderView1)

# hide data in view
Hide(resampleWithDataset1_1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(resampleWithDataset1)

# set active source
SetActiveSource(cOMSOLvtu)

# set active source
SetActiveSource(resampleWithDataset1)

# create a new 'Append Attributes'
appendAttributes1_1 = AppendAttributes(registrationName='AppendAttributes1', Input=resampleWithDataset1)

# set active source
SetActiveSource(resampleWithDataset1)

# destroy appendAttributes1_1
Delete(appendAttributes1_1)
del appendAttributes1_1

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(resampleWithDataset1_1)

# set active source
SetActiveSource(resampleWithDataset1)

# set active source
SetActiveSource(cOMSOLvtu)

# create a new 'Append Attributes'
appendAttributes1_1 = AppendAttributes(registrationName='AppendAttributes1', Input=[resampleWithDataset1, cOMSOLvtu])

RenameProxy(appendAttributes1_1, 'sources', 'XEMFEM_And_COMSOL')

# rename source object
RenameSource('XEMFEM_And_COMSOL', appendAttributes1_1)

# show data in view
appendAttributes1_1Display = Show(appendAttributes1_1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1_1Display.Representation = 'Surface'

# hide data in view
Hide(resampleWithDataset1, renderView1)

# hide data in view
Hide(cOMSOLvtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(appendAttributes1)

# set active source
SetActiveSource(simulationpvd)

# hide data in view
Hide(appendAttributes1, renderView1)

# show data in view
simulationpvdDisplay = Show(simulationpvd, renderView1, 'UnstructuredGridRepresentation')

# destroy appendAttributes1
Delete(appendAttributes1)
del appendAttributes1

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(resampleWithDataset1_1)

# set active source
SetActiveSource(simulationpvd)

# set active source
SetActiveSource(resampleWithDataset1_1)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=resampleWithDataset1_1)

# rename source object
RenameSource('Rename', calculator1)

# Properties modified on calculator1
calculator1.ResultArrayName = 'V_resampled'
calculator1.Function = 'V'

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'

# hide data in view
Hide(resampleWithDataset1_1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'V_resampled'
v_resampledLUT = GetColorTransferFunction('V_resampled')

# get opacity transfer function/opacity map for 'V_resampled'
v_resampledPWF = GetOpacityTransferFunction('V_resampled')

# get 2D transfer function for 'V_resampled'
v_resampledTF2D = GetTransferFunction2D('V_resampled')

# set active source
SetActiveSource(simulationpvd)

# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1', Input=[calculator1, simulationpvd])

RenameProxy(appendAttributes1, 'sources', 'XEMFEM_And_Resampled')

# rename source object
RenameSource('XEMFEM_And_Resampled', appendAttributes1)

# show data in view
appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1Display.Representation = 'Surface'

# hide data in view
Hide(calculator1, renderView1)

# hide data in view
Hide(simulationpvd, renderView1)

# show color bar/color legend
appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Calculator'
calculator1_1 = Calculator(registrationName='Calculator1', Input=appendAttributes1)

# Properties modified on calculator1_1
calculator1_1.Function = 'V-V_resampled'

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1_1Display.Representation = 'Surface'

# hide data in view
Hide(appendAttributes1, renderView1)

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'Result'
resultLUT = GetColorTransferFunction('Result')

# get opacity transfer function/opacity map for 'Result'
resultPWF = GetOpacityTransferFunction('Result')

# get 2D transfer function for 'Result'
resultTF2D = GetTransferFunction2D('Result')

# Properties modified on calculator1_1
calculator1_1.Function = 'abs(2*(V-V_resampled)/(V+V_resampled))'

# update the view to ensure updated data information
renderView1.Update()

# rename source object
RenameSource('SymErrResample', calculator1_1)

# set active source
SetActiveSource(appendAttributes1_1)

# create a new 'Calculator'
calculator1_2 = Calculator(registrationName='Calculator1', Input=appendAttributes1_1)

# rename source object
RenameSource('SymErrFrameworks', calculator1_2)

# set active source
SetActiveSource(appendAttributes1_1)

# set active source
SetActiveSource(calculator1_2)

# Properties modified on calculator1_2
calculator1_2.Function = 'abs(2*(V-Electric_potential)/(V+Electric_potential))'

# show data in view
calculator1_2Display = Show(calculator1_2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1_2Display.Representation = 'Surface'

# hide data in view
Hide(appendAttributes1_1, renderView1)

# show color bar/color legend
calculator1_2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_1)

# Properties modified on calculator1_1
calculator1_1.ResultArrayName = 'SampleErr'

# Properties modified on calculator1_2
calculator1_2.ResultArrayName = 'SymErr'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_2)

# set scalar coloring
ColorBy(calculator1_2Display, ('POINTS', 'SymErr'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(resultLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
calculator1_2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
calculator1_2Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'SymErr'
symErrLUT = GetColorTransferFunction('SymErr')

# get opacity transfer function/opacity map for 'SymErr'
symErrPWF = GetOpacityTransferFunction('SymErr')

# get 2D transfer function for 'SymErr'
symErrTF2D = GetTransferFunction2D('SymErr')

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(appendAttributes1_1)

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(calculator1_1)

# create a new 'Append Attributes'
appendAttributes1_2 = AppendAttributes(registrationName='AppendAttributes1', Input=[calculator1_2, calculator1_1])

# rename source object
RenameSource('CombineErrs', appendAttributes1_2)

# show data in view
appendAttributes1_2Display = Show(appendAttributes1_2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1_2Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_1, renderView1)

# hide data in view
Hide(calculator1_2, renderView1)

# show color bar/color legend
appendAttributes1_2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(appendAttributes1)

# set active source
SetActiveSource(appendAttributes1_2)

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(appendAttributes1_2)

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(appendAttributes1_2)

# set active source
SetActiveSource(calculator1_1)

# hide data in view
Hide(appendAttributes1_2, renderView1)

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# destroy appendAttributes1_2
Delete(appendAttributes1_2)
del appendAttributes1_2

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(calculator1_2)

# create a new 'Append Attributes'
appendAttributes1_2 = AppendAttributes(registrationName='AppendAttributes1', Input=[calculator1_1, calculator1_2])

# rename source object
RenameSource('ErrorsCombined', appendAttributes1_2)

# show data in view
appendAttributes1_2Display = Show(appendAttributes1_2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1_2Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_1, renderView1)

# hide data in view
Hide(calculator1_2, renderView1)

# show color bar/color legend
appendAttributes1_2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'SampleErr'
sampleErrLUT = GetColorTransferFunction('SampleErr')

# get opacity transfer function/opacity map for 'SampleErr'
sampleErrPWF = GetOpacityTransferFunction('SampleErr')

# get 2D transfer function for 'SampleErr'
sampleErrTF2D = GetTransferFunction2D('SampleErr')

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(appendAttributes1_2)

# set active source
SetActiveSource(calculator1_1)

# hide data in view
Hide(appendAttributes1_2, renderView1)

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# destroy appendAttributes1_2
Delete(appendAttributes1_2)
del appendAttributes1_2

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Calculator'
calculator1_3 = Calculator(registrationName='Calculator1', Input=calculator1_1)

# rename source object
RenameSource('ExtractErrOnly', calculator1_3)

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(calculator1_3)

# Properties modified on calculator1_3
calculator1_3.ResultArrayName = 'SampleErr'
calculator1_3.Function = 'SampleErr'

# show data in view
calculator1_3Display = Show(calculator1_3, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1_3Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_1, renderView1)

# show color bar/color legend
calculator1_3Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_2)

# create a new 'Calculator'
calculator1_4 = Calculator(registrationName='Calculator1', Input=calculator1_2)

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(calculator1_4)

# rename source object
RenameSource('ExtracErr', calculator1_4)

# Properties modified on calculator1_4
calculator1_4.ResultArrayName = 'SymErr'
calculator1_4.Function = 'SymErr'

# show data in view
calculator1_4Display = Show(calculator1_4, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1_4Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_2, renderView1)

# show color bar/color legend
calculator1_4Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_3)

# create a new 'Append Attributes'
appendAttributes1_2 = AppendAttributes(registrationName='AppendAttributes1', Input=[calculator1_4, calculator1_3])

# rename source object
RenameSource('CombineErrs', appendAttributes1_2)

# show data in view
appendAttributes1_2Display = Show(appendAttributes1_2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1_2Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_3, renderView1)

# hide data in view
Hide(calculator1_4, renderView1)

# show color bar/color legend
appendAttributes1_2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_3)

# set active source
SetActiveSource(appendAttributes1_2)

# set active source
SetActiveSource(calculator1_3)

# hide data in view
Hide(appendAttributes1_2, renderView1)

# show data in view
calculator1_3Display = Show(calculator1_3, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_3Display.SetScalarBarVisibility(renderView1, True)

# destroy appendAttributes1_2
Delete(appendAttributes1_2)
del appendAttributes1_2

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(calculator1_1)

# hide data in view
Hide(calculator1_3, renderView1)

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# destroy calculator1_3
Delete(calculator1_3)
del calculator1_3

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(calculator1_4)

# set active source
SetActiveSource(calculator1_2)

# destroy calculator1_4
Delete(calculator1_4)
del calculator1_4

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(calculator1_1)

# create a new 'Append Datasets'
appendDatasets1 = AppendDatasets(registrationName='AppendDatasets1', Input=[calculator1_2, calculator1_1])

# rename source object
RenameSource('CombineErrs', appendDatasets1)

# show data in view
appendDatasets1Display = Show(appendDatasets1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendDatasets1Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_1, renderView1)

# hide data in view
Hide(calculator1_2, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_1)

# hide data in view
Hide(appendDatasets1, renderView1)

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# destroy appendDatasets1
Delete(appendDatasets1)
del appendDatasets1

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(calculator1_1)

# set active source
SetActiveSource(calculator1_2)

# create a new 'Append Attributes'
appendAttributes1_2 = AppendAttributes(registrationName='AppendAttributes1', Input=[calculator1_1, calculator1_2])

# show data in view
appendAttributes1_2Display = Show(appendAttributes1_2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
appendAttributes1_2Display.Representation = 'Surface'

# hide data in view
Hide(calculator1_1, renderView1)

# hide data in view
Hide(calculator1_2, renderView1)

# show color bar/color legend
appendAttributes1_2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(appendAttributes1_2)

# set active source
SetActiveSource(calculator1_1)

# hide data in view
Hide(appendAttributes1_2, renderView1)

# show data in view
calculator1_1Display = Show(calculator1_1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# destroy appendAttributes1_2
Delete(appendAttributes1_2)
del appendAttributes1_2

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Histogram'
histogram1 = Histogram(registrationName='Histogram1', Input=calculator1_1)

# Properties modified on histogram1
histogram1.BinCount = 100

# Create a new 'Bar Chart View'
barChartView1 = CreateView('XYBarChartView')

# show data in view
histogram1Display = Show(histogram1, barChartView1, 'XYBarChartRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=barChartView1, layout=layout1, hint=0)

# Properties modified on histogram1Display
histogram1Display.SeriesOpacity = ['bin_extents', '1', 'bin_values', '1']

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator1_1)

# hide data in view
Hide(histogram1, barChartView1)

# destroy histogram1
Delete(histogram1)
del histogram1

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# destroy barChartView1
Delete(barChartView1)
del barChartView1

# close an empty frame
layout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# set scalar coloring
ColorBy(calculator1_1Display, ('POINTS', 'SampleErr'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(resultLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
calculator1_1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
calculator1_1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(appendAttributes1)

# hide data in view
Hide(calculator1_1, renderView1)

# show data in view
appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

# destroy calculator1_1
Delete(calculator1_1)
del calculator1_1

# set active source
SetActiveSource(calculator1)

# hide data in view
Hide(appendAttributes1, renderView1)

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# destroy appendAttributes1
Delete(appendAttributes1)
del appendAttributes1

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(resampleWithDataset1_1)

# hide data in view
Hide(calculator1, renderView1)

# show data in view
resampleWithDataset1_1Display = Show(resampleWithDataset1_1, renderView1, 'UnstructuredGridRepresentation')

# destroy calculator1
Delete(calculator1)
del calculator1

# set active source
SetActiveSource(simulationpvd)

# hide data in view
Hide(resampleWithDataset1_1, renderView1)

# show data in view
simulationpvdDisplay = Show(simulationpvd, renderView1, 'UnstructuredGridRepresentation')

# destroy resampleWithDataset1_1
Delete(resampleWithDataset1_1)
del resampleWithDataset1_1

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(appendAttributes1_1)

# set active source
SetActiveSource(calculator1_2)

# set active source
SetActiveSource(calculator1_2)

# show data in view
calculator1_2Display = Show(calculator1_2, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator1_2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(simulationpvd, renderView1)

# create a new 'Histogram'
histogram1 = Histogram(registrationName='Histogram1', Input=calculator1_2)

# Properties modified on histogram1
histogram1.BinCount = 1000

# Create a new 'Bar Chart View'
barChartView1 = CreateView('XYBarChartView')

# show data in view
histogram1Display = Show(histogram1, barChartView1, 'XYBarChartRepresentation')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=barChartView1, layout=layout1, hint=0)

# Properties modified on histogram1Display
histogram1Display.SeriesOpacity = ['bin_extents', '1', 'bin_values', '1']

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on histogram1
histogram1.BinCount = 100

# update the view to ensure updated data information
barChartView1.Update()

# set active source
SetActiveSource(simulationpvd)

# set active source
SetActiveSource(cOMSOLvtu)

# set active source
SetActiveSource(simulationpvd)

# set active source
SetActiveSource(calculator1_2)

# create a new 'Descriptive Statistics'
descriptiveStatistics1 = DescriptiveStatistics(registrationName='DescriptiveStatistics1', Input=calculator1_2,
    ModelInput=None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
descriptiveStatistics1Display = Show(descriptiveStatistics1, spreadSheetView1, 'SpreadSheetRepresentation')

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=2)

# update the view to ensure updated data information
spreadSheetView1.Update()

# Properties modified on descriptiveStatistics1
descriptiveStatistics1.VariablesofInterest = ['SymErr']

# update the view to ensure updated data information
spreadSheetView1.Update()

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(descriptiveStatistics1)

# set active source
SetActiveSource(calculator1_2)

# change representation type
calculator1_2Display.SetRepresentationType('Surface With Edges')

# change representation type
calculator1_2Display.SetRepresentationType('Surface')

# get color legend/bar for symErrLUT in view renderView1
symErrLUTColorBar = GetScalarBar(symErrLUT, renderView1)

symErrLUTColorBar.WindowLocation = 'Any Location'
symErrLUTColorBar.Position = [0.7672134935304991, 0.5534709193245779]
symErrLUTColorBar.ScalarBarLength = 0.32999999999999996

layout1.SetSplitFraction(0, 0.8340929808568824)

layout1.SetSplitFraction(0, 0.9279854147675478)

layout1.SetSize(2419, 1066)

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4570489784651149, -0.49411320893318217, 4.505813621664205]
renderView1.CameraFocalPoint = [0.4570489784651149, -0.49411320893318217, 0.0]
renderView1.CameraParallelScale = 0.7965237203531584

```