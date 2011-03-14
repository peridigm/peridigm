
include_directories(${CMAKE_SOURCE_DIR}/io ${TRILINOS_DIR}/include)

#
add_library(magOp2D 
quadrulator.C
Mag2D_Operator.C
)

add_executable(utQuadrulator utQuadrulator.cxx)
target_link_libraries(utQuadrulator magOp2D ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

