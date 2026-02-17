# src/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

project(MFEMElectrostatics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Put all build artifacts into the top-level build dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR})

# Common libraries
add_subdirectory(cmdLineInteraction)
add_subdirectory(config)

# Modules / executables
add_subdirectory(SolverCore)
#add_subdirectory(plotting)
add_subdirectory(interpolator)
add_subdirectory(Sweep)
add_subdirectory(xemfem) 
add_subdirectory(path_handler) 