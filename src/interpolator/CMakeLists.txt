# src/config
add_library(interpolator STATIC
    interpolator.cpp
    interpolator.h
)

find_package(HDF5 REQUIRED COMPONENTS CXX)

target_include_directories(interpolator
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        config
        SolverCore
        HDF5::HDF5
)

target_link_libraries(interpolator
    PUBLIC
        config
        SolverCore
        HDF5::HDF5
)