#
# Initial Conditions Library for Peridigm
#

# Library of Initial conditions
include_directories(${CMAKE_SOURCE_DIR}/io/utilities)
add_library(PdIC InitialCondition.cpp Q2Cylinder.cpp)


