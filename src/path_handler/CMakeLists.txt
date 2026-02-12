# src/config
add_library(path_handler STATIC
    path_handler.cpp
    path_handler.h
)

find_package(yaml-cpp REQUIRED)

target_include_directories(path_handler
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        config
)

target_link_libraries(path_handler
    PUBLIC
        config
)