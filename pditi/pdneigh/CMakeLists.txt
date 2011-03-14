# Library for Peridynamics neighborhood creation

add_subdirectory(unit_test)

# include this path
include_directories(${UTILITIES_DIR} ${MESH_INPUT_DIR} ${TRILINOS_DIR}/include )
add_library(PdNeigh NeighborhoodList.cxx PdZoltan.cxx BondFilter.cxx)
