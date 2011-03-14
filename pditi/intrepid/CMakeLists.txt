#
# Intrepid Interface
#

# add unit_test directory
add_subdirectory(unit_test)

include_directories(${TRILINOS_DIR}/include)
add_library(IIntrepid FIFE_2D_IntrepidQuadrature.C)
