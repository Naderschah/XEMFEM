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
# --- Config (YAML loader) --------------------------------
# ---------------------------------------------------------
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../config
                 ${CMAKE_CURRENT_BINARY_DIR}/config)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cmdLineParser
                 ${CMAKE_CURRENT_BINARY_DIR}/cmdLineParser)
                 

# ---------------------------------------------------------
# --- Solver executable (no geometry coupling) ------------
# ---------------------------------------------------------
add_executable(SOLVER
  main.cpp
  load_mesh.cpp
  parallelization.cpp
  boundary_conditions.cpp
  solver.cpp
  solver_api.cpp
  ComputeElectricField.cpp
)

target_include_directories(SOLVER PRIVATE
  ${CMAKE_SOURCE_DIR}
  ${MFEM_INCLUDE_DIRS}
)

target_link_libraries(SOLVER PRIVATE
  mfem
  config
  cmdLineParser
  yaml-cpp
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
)
