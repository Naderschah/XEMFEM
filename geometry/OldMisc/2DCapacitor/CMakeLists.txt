cmake_minimum_required(VERSION 3.14)
project(2DCapacitor LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Your sources
add_executable(2DCapacitor
  main.cpp
  ../../src/config/Config.cpp
)

target_include_directories(2DCapacitor PRIVATE
  ../../src/config
)

# yaml-cpp if you use it
find_package(yaml-cpp REQUIRED)
find_library(GMSH_LIB gmsh REQUIRED)

target_link_libraries(2DCapacitor PRIVATE yaml-cpp ${GMSH_LIB})

# --- Embed config.yaml into a header
set(CONFIG_YAML ${CMAKE_CURRENT_SOURCE_DIR}/config.yaml)
file(READ "${CONFIG_YAML}" CONFIG_YAML_CONTENT)

# Generate header with the YAML embedded as a raw string
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/embedded_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/embedded_config.h
  @ONLY
)

# Make the generated header visible
target_include_directories(2DCapacitor PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
