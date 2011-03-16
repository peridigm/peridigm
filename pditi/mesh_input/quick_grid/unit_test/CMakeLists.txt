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
