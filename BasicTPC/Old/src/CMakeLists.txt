cmake_minimum_required(VERSION 3.15)
project(TPC)

set(CMAKE_CXX_STANDARD 17)

# --- Populate the geo.in file with constants from constants.h and produce the geo file---

# Load constants and expose the constexpr 
file(READ ${CMAKE_SOURCE_DIR}/constants.h CONST_CONTENTS)
# Match lines like: constexpr double plate_gap = 0.5;
string(REGEX MATCHALL
    "constexpr[ \t]+[A-Za-z0-9_<>:]+[ \t]+([A-Za-z0-9_]+)[ \t]*=[ \t]*([-+0-9.eE]+)"
    CONST_MATCHES "${CONST_CONTENTS}"
)

# Format the values 
foreach(match ${CONST_MATCHES})
    string(REGEX REPLACE
        "constexpr[ \t]+[A-Za-z0-9_<>:]+[ \t]+([A-Za-z0-9_]+)[ \t]*=.*"
        "\\1" name "${match}")
    string(REGEX REPLACE
        ".*=[ \t]*([-+0-9.eE]+).*"
        "\\1" value "${match}")
    set(${name} "${value}" CACHE INTERNAL "")
endforeach()

# Read the .geo.in file 
set(GEO_TEMPLATE_FILE ${CMAKE_SOURCE_DIR}/geometry.geo.in)
file(READ ${GEO_TEMPLATE_FILE} GEO_TEMPLATE_CONTENTS)

# Validate that each token required exists in constants
foreach(token ${GEO_TOKENS})
	string(REGEX REPLACE "@" "" var "${token}")
	if (NOT DEFINED ${var})
		message(FATAL_ERROR
			"Variable ${token} used in the geometry file is not defined in constants.h")
	endif()
endforeach()

# Configure the .geo file
configure_file(${GEO_TEMPLATE_FILE} ${CMAKE_BINARY_DIR}/geometry.geo @ONLY)
# --- End Geo Matching ---

# --- Generate the mesh file from the geo file ---

set(GMSH_OUTPUT ${CMAKE_BINARY_DIR}/geometry.msh)
set(GMSH_INPUT  ${CMAKE_BINARY_DIR}/geometry.geo)

find_program(GMSH_EXECUTABLE gmsh REQUIRED)

add_custom_command(
  OUTPUT ${GMSH_OUTPUT}
  COMMAND ${GMSH_EXECUTABLE} -3 ${GMSH_INPUT} -o ${GMSH_OUTPUT}
  DEPENDS ${GMSH_INPUT}
  COMMENT "Generating mesh: ${GMSH_OUTPUT}"
  VERBATIM
)

# In case we only want to make the mesh
add_custom_target(mesh ALL DEPENDS ${GMSH_OUTPUT})

# --- End generating mesh file ---


add_executable(TPC
	main.cpp
	geometry.cpp
	boundary_conditions.cpp
	solver.cpp
	ComputeElectricField.cpp
)

target_link_libraries(TPC mfem )

add_dependencies(TPC mesh)