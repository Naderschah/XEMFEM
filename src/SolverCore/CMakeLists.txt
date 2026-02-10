# src/SolverCore/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# ---------------------------------------------------------
# --- Standalone or subproject detection ------------------
# ---------------------------------------------------------
set(IS_SOLVERCORE_TOPLEVEL OFF)
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IS_SOLVERCORE_TOPLEVEL ON)
    project(SolverCoreProj LANGUAGES CXX)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# ---------------------------------------------------------
# --- Third-party packages (needed in any mode) -----------
# ---------------------------------------------------------
find_package(yaml-cpp REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(MFEM_DIR "/home/felix/mfem-install/lib/cmake/mfem")
find_package(MFEM CONFIG REQUIRED)

# ---------------------------------------------------------
# --- Dependencies only when building SolverCore alone ----
# ---------------------------------------------------------
if (IS_SOLVERCORE_TOPLEVEL)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../config
                     ${CMAKE_CURRENT_BINARY_DIR}/config)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cmdLineParser
                     ${CMAKE_CURRENT_BINARY_DIR}/cmdLineParser)
endif()

# ---------------------------------------------------------
# --- SolverCore library ----------------------------------
# ---------------------------------------------------------
add_library(SolverCore STATIC
  mesh_tools.cpp
  parallelization.cpp
  boundary_conditions.cpp
  solver.cpp
  solver_api.cpp
  ComputeElectricField.cpp
)

target_include_directories(SolverCore
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MFEM_INCLUDE_DIRS}
)

target_link_libraries(SolverCore
  PUBLIC
    mfem
    config
    cmdLineParser
    yaml-cpp
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)
