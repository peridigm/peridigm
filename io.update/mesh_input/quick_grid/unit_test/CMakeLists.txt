#
# QuickGrid Unit Tests
#

include_directories(${UTILITIES_DIR} ${PDNEIGH_DIR} ${MESH_INPUT_DIR} ${MESH_OUTPUT_DIR} ${TRILINOS_DIR}/include)

add_executable(ut_Q2CylinderRingHorizon ut_Q2CylinderRingHorizon.cxx)
target_link_libraries(ut_Q2CylinderRingHorizon  QuickGrid ${UT_REQUIRED_LIBS})
add_test (ut_Q2CylinderRingHorizon python run_unit_test.py ./ut_Q2CylinderRingHorizon)

add_executable(ut_QuickGridHorizon ut_QuickGridHorizon.cxx)
target_link_libraries(ut_QuickGridHorizon  QuickGrid ${UT_REQUIRED_LIBS})
add_test (ut_QuickGridHorizon python run_unit_test.py ./ut_QuickGridHorizon)

add_executable(ut_QuickGrid_solidCylinder_np1 ut_QuickGrid_solidCylinder_np1.cxx)
target_link_libraries(ut_QuickGrid_solidCylinder_np1  QuickGrid PdVTK ${VTK_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGrid_solidCylinder_np1 python run_unit_test.py ./ut_QuickGrid_solidCylinder_np1)

add_executable(ut_QuickGrid_solidCylinder_np2 ut_QuickGrid_solidCylinder_np2.cxx)
target_link_libraries(ut_QuickGrid_solidCylinder_np2  QuickGrid PdVTK PdNeigh ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGrid_solidCylinder_np2 python run_unit_test.py mpiexec -np 2 ./ut_QuickGrid_solidCylinder_np2)

add_executable(ut_QuickGridMPI_np2_3x3x2 ut_QuickGridMPI_np2_3x3x2.cxx)
target_link_libraries(ut_QuickGridMPI_np2_3x3x2  QuickGrid PdNeigh ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGridMPI_np2_3x3x2 python run_unit_test.py mpiexec -np 2 ./ut_QuickGridMPI_np2_3x3x2)

add_executable(ut_QuickGridMPI_np3 ut_QuickGridMPI_np3.cxx)
target_link_libraries(ut_QuickGridMPI_np3  QuickGrid PdNeigh ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGridMPI_np3 python run_unit_test.py mpiexec -np 3 ./ut_QuickGridMPI_np3)

add_executable(ut_SmallMeshCylinder_np4 ut_SmallMeshCylinder_np4.cxx)
target_link_libraries(ut_SmallMeshCylinder_np4  QuickGrid PdNeigh ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_SmallMeshCylinder_np4 python run_unit_test.py mpiexec -np 4 ./ut_SmallMeshCylinder_np4)

add_executable(ut_SmallMeshCylinder ut_SmallMeshCylinder.cxx)
target_link_libraries(ut_SmallMeshCylinder  QuickGrid PdNeigh ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_SmallMeshCylinder python run_unit_test.py ./ut_SmallMeshCylinder)


