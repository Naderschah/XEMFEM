# src/cmdLineParser
add_library(cmdLineParser STATIC
    cmdLineParser.cpp
    cmdLineParser.h
)

target_include_directories(cmdLineParser
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)