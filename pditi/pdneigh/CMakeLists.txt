# Library for Peridynamics neighborhood creation


# include this path
include_directories(${UTILITIES_DIR} ${TRILINOS_DIR}/include )
add_library(PdNeigh NeighborhoodList.cxx)
