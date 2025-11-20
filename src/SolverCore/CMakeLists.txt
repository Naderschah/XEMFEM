cmake_minimum_required(VERSION 3.15)
project(SOLVER LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---------------------------------------------------------
# --- Third-party packages (found via hardcoded prefix) ---
# ---------------------------------------------------------
find_package(yaml-cpp REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(MFEM_DIR "/home/felix/mfem-install/lib/cmake/mfem")
find_package(MFEM CONFIG REQUIRED)

# ---------------------------------------------------------
# --- Explicitly include Hypre (Ubuntu dev package) -------
#     Ensures HYPRE_utilities.h is always visible
# ---------------------------------------------------------
find_path(HYPRE_INCLUDE_DIR HYPRE_utilities.h
          PATHS /usr/include PATH_SUFFIXES hypre REQUIRED)
find_library(HYPRE_LIBRARY NAMES HYPRE REQUIRED)

# ---------------------------------------------------------

# ---------------------------------------------------------
# --- Config (YAML loader) --------------------------------
# ---------------------------------------------------------
add_subdirectory(../config)
add_subdirectory(../cmdLineParser)

# ---------------------------------------------------------
# --- Solver executable (no geometry coupling) ------------
# ---------------------------------------------------------
add_executable(SOLVER
  main.cpp
  load_mesh.cpp
  boundary_conditions.cpp
  solver.cpp
  ComputeElectricField.cpp
  cmdLineParser.cpp
)

# Our own includes
target_include_directories(SOLVER PRIVATE
  ${CMAKE_SOURCE_DIR}
  ${HYPRE_INCLUDE_DIR}      # ensure HYPRE_utilities.h visible during compile
  ${MFEM_INCLUDE_DIRS}
)

target_link_libraries(SOLVER PRIVATE
  ${HYPRE_LIBRARY}
  mfem
  config
  cmdLineParser
  yaml-cpp
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
)

# ---------------------------------------------------------
# --- Geometry (optional, separate build) -----------------
# ---------------------------------------------------------
# TODO Not sure how to handle this properly. Its not part of the solver
# And will likely only be done to generate geometries
# add_subdirectory(geometry)
