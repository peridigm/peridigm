
include_directories(${UTILITIES_DIR} ${MESH_OUTPUT_DIR} ${VTK_INCLUDE_DIRS} ${TRILINOS_DIR}/include)

#
add_library(magOp2D 
quadrulator.C
Mag2D_Operator.C
)

add_executable(utQuadrulator utQuadrulator.cxx)
target_link_libraries(utQuadrulator magOp2D ${MESH_OUTPUT_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utQuadrulator python run_unit_test.py ./utQuadrulator)
