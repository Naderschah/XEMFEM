# src/config
add_library(config STATIC
    Config.cpp
    Config.h
)

target_include_directories(config
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)