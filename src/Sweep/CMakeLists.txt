cmake_minimum_required(VERSION 3.15)
project(SWEEP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---------------------------------------------------------
# --- Third-party packages (same as SolverCore) ----------
# ---------------------------------------------------------
find_package(yaml-cpp REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(MFEM_DIR "/home/felix/mfem-install/lib/cmake/mfem")
find_package(MFEM CONFIG REQUIRED)

# ---------------------------------------------------------
# --- Config (YAML loader) --------------------------------
# ---------------------------------------------------------
# Reuse the config and cmdLineParser subprojects from src/
# Adjust relative paths if your layout is different.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../config
                 ${CMAKE_CURRENT_BINARY_DIR}/config)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cmdLineParser
                 ${CMAKE_CURRENT_BINARY_DIR}/cmdLineParser)

# ---------------------------------------------------------
# --- Sweep executable ------------------------------------
# ---------------------------------------------------------
# main.cpp is your sweep main in src/Sweep
# The other sources are the shared solver implementation from SolverCore.
add_executable(SWEEP
  main.cpp

  # From SolverCore
  ../SolverCore/load_mesh.cpp
  ../SolverCore/parallelization.cpp
  ../SolverCore/boundary_conditions.cpp
  ../SolverCore/solver.cpp
  ../SolverCore/solver_api.cpp
  ../SolverCore/ComputeElectricField.cpp

  # From this module
  config_modification.cpp
  optimization.cpp
  trace_fieldlines.cpp
)

# ------ Nelder Mead (header-only) ---------------------------
include(FetchContent)

FetchContent_Declare(
    neldermead
    GIT_REPOSITORY https://github.com/YibaiMeng/nelder-mead.git
    GIT_TAG        master   # or a specific commit/tag
)

FetchContent_GetProperties(neldermead)
if(NOT neldermead_POPULATED)
    FetchContent_Populate(neldermead)
endif()


# Include paths:
#  - MFEM headers
#  - the src root (so that #include "Config.h" etc. work)
#  - SolverCore (for solver_api.h etc.)
target_include_directories(SWEEP PRIVATE
  ${MFEM_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..          # src/
  ${CMAKE_CURRENT_SOURCE_DIR}/../SolverCore
  ${neldermead_SOURCE_DIR}
)

# Link against the same libs as SOLVER
target_link_libraries(SWEEP PRIVATE
  mfem
  config
  cmdLineParser
  yaml-cpp
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
)
