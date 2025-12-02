# src/Sweep/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

# ---------------------------------------------------------
# --- Standalone or subproject detection ------------------
# ---------------------------------------------------------
set(IS_SWEEP_TOPLEVEL OFF)
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(IS_SWEEP_TOPLEVEL ON)
    project(SWEEP LANGUAGES CXX)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# ---------------------------------------------------------
# --- Third-party packages (any mode) ---------------------
# ---------------------------------------------------------
find_package(yaml-cpp REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(MFEM_DIR "/home/felix/mfem-install/lib/cmake/mfem")
find_package(MFEM CONFIG REQUIRED)

# ---------------------------------------------------------
# --- Dependencies when Sweep is standalone ---------------
# ---------------------------------------------------------
if (IS_SWEEP_TOPLEVEL)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../config
                     ${CMAKE_CURRENT_BINARY_DIR}/config)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cmdLineParser
                     ${CMAKE_CURRENT_BINARY_DIR}/cmdLineParser)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../SolverCore
                     ${CMAKE_CURRENT_BINARY_DIR}/SolverCore)
endif()

# ---------------------------------------------------------
# --- Sweep executable ------------------------------------
# ---------------------------------------------------------
add_executable(SWEEP
  main.cpp
  config_modification.cpp
  optimization.cpp
  trace_fieldlines.cpp
  field_spread.cpp
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

target_include_directories(SWEEP PRIVATE
  ${MFEM_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..          # src/
  ${CMAKE_CURRENT_SOURCE_DIR}/../SolverCore
  ${neldermead_SOURCE_DIR}
)

target_link_libraries(SWEEP PRIVATE
  SolverCore
  # SolverCore already links PUBLIC to:
  # mfem, config, cmdLineParser, yaml-cpp, MPI::MPI_CXX, OpenMP::OpenMP_CXX
)
