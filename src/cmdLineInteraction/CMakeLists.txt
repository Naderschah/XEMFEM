# src/cmdLineInteraction
add_library(cmdLineInteraction STATIC
    cmdLineInteraction.cpp
    cmdLineInteraction.h
)

target_include_directories(cmdLineInteraction
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)