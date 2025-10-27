cmake_minimum_required(VERSION 3.15)
project(TPC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ---------------------------------------------------------
# --- Build geometry generator (geometry/main.cpp) --------
# ---------------------------------------------------------
add_executable(build_geometry
    ${CMAKE_SOURCE_DIR}/geometry/main.cpp
    ${CMAKE_SOURCE_DIR}/geometry/partition_tools.cpp
    ${CMAKE_SOURCE_DIR}/geometry/debugging_things.cpp
)
# Link against gmsh (must be installed system-wide or via package)
find_library(GMSH_LIB gmsh REQUIRED)
target_link_libraries(build_geometry PRIVATE ${GMSH_LIB})
target_include_directories(build_geometry PRIVATE 
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/geometry)

# ---------------------------------------------------------
# --- Run geometry generator to produce .msh file ----------
# ---------------------------------------------------------
set(GEOM_SRC_DIR ${CMAKE_SOURCE_DIR}/geometry)
set(GEOM_MSH_SRC ${GEOM_SRC_DIR}/tpc_occ.msh)
set(GEOM_MSH_BIN ${CMAKE_BINARY_DIR}/geometry.msh)

add_custom_command(
  OUTPUT ${GEOM_MSH_BIN}
  COMMAND build_geometry                              # run the geometry builder
  COMMAND ${CMAKE_COMMAND} -E copy_if_different        # copy .msh into build directory
          ${GEOM_MSH_SRC} ${GEOM_MSH_BIN}
  WORKING_DIRECTORY ${GEOM_SRC_DIR}                    # run inside geometry folder
  DEPENDS build_geometry 
          ${CMAKE_SOURCE_DIR}/geometry/main.cpp
          ${CMAKE_SOURCE_DIR}/geometry/partition_tools.cpp
          ${CMAKE_SOURCE_DIR}/geometry/partition_tools.h
          ${CMAKE_SOURCE_DIR}/geometry/debugging_things.cpp
          ${CMAKE_SOURCE_DIR}/geometry/debugging_things.h
          ${CMAKE_SOURCE_DIR}/constants.h
  BYPRODUCTS ${GEOM_MSH_SRC}
  COMMENT "Running geometry generator and copying tpc_occ.msh into build directory"
  VERBATIM
)

# convenience target so you can run: make mesh
add_custom_target(mesh ALL DEPENDS ${GEOM_MSH_BIN})

# ---------------------------------------------------------
# --- Main TPC simulation executable ----------------------
# ---------------------------------------------------------
add_executable(TPC
  main.cpp
  geometry.cpp
  boundary_conditions.cpp
  solver.cpp
  ComputeElectricField.cpp
)

target_link_libraries(TPC PRIVATE mfem)

# Make sure the mesh is ready before building/running TPC
add_dependencies(TPC mesh)
