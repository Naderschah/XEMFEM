cmake_minimum_required(VERSION 3.15)
project(PLOT LANGUAGES CXX)

# Point CMake at your VTK install
set(vtk_DIR "$ENV{HOME}/vtk/install/lib/cmake/vtk-9.5" CACHE PATH "Path to vtk-config.cmake")

# Use the new-style config package "vtk"
find_package(vtk CONFIG REQUIRED)

# Library with your plotting implementation
add_library(plotting STATIC
    FEMPlotVTK.cpp
)

target_include_directories(plotting
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link against the VTK imported targets you actually use.
# Make sure these match the modules you enabled when building VTK.
target_link_libraries(plotting
    PUBLIC
        vtk::CommonCore
        vtk::CommonDataModel
        vtk::IOXML
        vtk::IOImage
        vtk::FiltersCore
        vtk::FiltersGeneral
        vtk::FiltersGeometry
        vtk::FiltersSources
        vtk::FiltersFlowPaths
        vtk::RenderingCore
        vtk::RenderingOpenGL2
        vtk::RenderingFreeType
        vtk::RenderingContext2D
        vtk::RenderingContextOpenGL2
        vtk::RenderingAnnotation
        vtk::ChartsCore
        vtk::InteractionStyle
)

# Executable that uses the plotting library
add_executable(PLOT
    main.cpp
)

target_link_libraries(PLOT
    PRIVATE
        plotting
)
