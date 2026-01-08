# src/xemfem/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)

add_executable(XEMFEM main.cpp)

set_target_properties(XEMFEM PROPERTIES
  OUTPUT_NAME "XEMFEM"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

target_include_directories(XEMFEM PRIVATE
  ${MFEM_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..          # src/
  ${CMAKE_CURRENT_SOURCE_DIR}/../Sweep/
  ${CMAKE_CURRENT_SOURCE_DIR}/../config/
)

target_link_libraries(XEMFEM PRIVATE
  SolverCore
  plotting
  sweep_objects
)
