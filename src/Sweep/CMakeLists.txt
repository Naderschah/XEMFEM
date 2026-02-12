# src/Sweep/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

set(IS_SWEEP_TOPLEVEL OFF)
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IS_SWEEP_TOPLEVEL ON)
    project(SWEEP LANGUAGES CXX)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

find_package(yaml-cpp REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

if (IS_SWEEP_TOPLEVEL)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../config
                     ${CMAKE_CURRENT_BINARY_DIR}/config)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cmdLineInteraction
                     ${CMAKE_CURRENT_BINARY_DIR}/cmdLineInteraction)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../SolverCore
                     ${CMAKE_CURRENT_BINARY_DIR}/SolverCore)
endif()

# ------ Nelder Mead (sometimes network fails for this for some reason) ---------------------------
set(NELDERMEAD_LOCAL "/home/felix/third_party/nelder-mead")
if(EXISTS "${NELDERMEAD_LOCAL}/nelder-mead.h")
  add_library(neldermead_headers INTERFACE)
  target_include_directories(neldermead_headers INTERFACE "${NELDERMEAD_LOCAL}")
else()
  include(FetchContent)
  FetchContent_Declare(
      neldermead
      GIT_REPOSITORY https://github.com/YibaiMeng/nelder-mead.git
      GIT_TAG master
  )
  FetchContent_GetProperties(neldermead)
  if(NOT neldermead_POPULATED)
      FetchContent_Populate(neldermead)
  endif()
  add_library(neldermead_headers INTERFACE)
  target_include_directories(neldermead_headers INTERFACE "${neldermead_SOURCE_DIR}")
endif()

# ------ Indicators (header-only progress bar) --------------------------------
set(INDICATORS_LOCAL "/home/felix/third_party/indicators")

if(EXISTS "${INDICATORS_LOCAL}/indicators/progress_bar.hpp")
  add_library(indicators_headers INTERFACE)
  target_include_directories(indicators_headers INTERFACE "${INDICATORS_LOCAL}")
else()
  include(FetchContent)
  FetchContent_Declare(
      indicators
      GIT_REPOSITORY https://github.com/p-ranav/indicators.git
      GIT_TAG master
  )
  FetchContent_GetProperties(indicators)
  if(NOT indicators_POPULATED)
      FetchContent_Populate(indicators)
  endif()

  add_library(indicators_headers INTERFACE)
  target_include_directories(indicators_headers INTERFACE
      "${indicators_SOURCE_DIR}/include"
  )
endif()

# ------ Stepping (header-only) ---------------------------
# We need Boost Odeint
include(FetchContent)
FetchContent_Declare(
    boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.tar.gz
)
FetchContent_GetProperties(boost)
if(NOT boost_POPULATED)
    FetchContent_Populate(boost)
endif()

# ---------------------------------------------------------
# --- Sweep implementation as an OBJECT library ------------
# ---------------------------------------------------------
add_library(sweep_objects OBJECT
  config_modification.cpp
  optimization.cpp
  field_spread.cpp
  sweeps.cpp
  trace_fieldlines.cpp
  trace_fieldlines_MPI.cpp
  CIV.cpp
  common_tracing/tracing_objects.cpp
)
set(VTK_DIR "$ENV{HOME}/vtk/install/lib/cmake/vtk-9.5" CACHE PATH "Path to VTK config")

find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonColor
    CommonDataModel
    IOXML
    IOImage
    FiltersCore
    FiltersGeneral
    FiltersGeometry
    FiltersExtraction
    FiltersSources
    FiltersFlowPaths
    RenderingCore
    RenderingOpenGL2
    RenderingFreeType
    RenderingContext2D
    RenderingContextOpenGL2
    RenderingAnnotation
    ChartsCore
    InteractionStyle
    FiltersPoints
    RenderingGridAxes
    FiltersExtraction
)

if (NOT VTK_FOUND)
    message(FATAL_ERROR "PLOT: Unable to find the VTK build folder.")
endif()


target_include_directories(sweep_objects PRIVATE
  ${MFEM_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..          # src/
  ${CMAKE_CURRENT_SOURCE_DIR}/../SolverCore
  ${CMAKE_CURRENT_SOURCE_DIR}/../SolverCore/common_tracing
  ${neldermead_SOURCE_DIR}
  ${boost_SOURCE_DIR}
)

target_link_libraries(sweep_objects PRIVATE
  ${VTK_LIBRARIES}
  SolverCore
  neldermead_headers
  interpolator
  indicators_headers
  mfem
  path_handler
)

