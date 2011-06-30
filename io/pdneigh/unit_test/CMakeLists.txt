
#
# unit tests for neighborhood library
# 

# include this path
include_directories(${MESH_INPUT_DIR} ${MESH_OUTPUT_DIR} ${UTILITIES_DIR} ${TRILINOS_DIR}/include)

#
add_executable(ut_frameset_2x2x1_np4 ut_frameset_2x2x1_np4.cxx)
target_link_libraries(ut_frameset_2x2x1_np4  PdNeigh QuickGrid PdVTK Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_frameset_2x2x1_np4 python run_unit_test.py mpiexec -np 4 ./ut_frameset_2x2x1_np4)

#
add_executable(ut_neighborhood_list ut_neighborhood_list.cxx)
target_link_libraries(ut_neighborhood_list  PdNeigh QuickGrid PdVTK Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_neighborhood_list python run_unit_test.py ./ut_neighborhood_list)

#
add_executable(ut_Y-Z_crack ut_Y-Z_crack.cxx)
target_link_libraries(ut_Y-Z_crack  PdNeigh QuickGrid PdVTK Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_Y-Z_crack python run_unit_test.py ./ut_Y-Z_crack)

#
add_executable(ut_Y-Z_crack_np4 ut_Y-Z_crack_np4.cxx)
target_link_libraries(ut_Y-Z_crack_np4  PdNeigh QuickGrid PdVTK Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_Y-Z_crack_np4 python run_unit_test.py mpiexec -np 4 ./ut_Y-Z_crack_np4)

#
add_executable(ut_Y-Z_crack_jacobian_np4 ut_Y-Z_crack_jacobian_np4.cxx)
target_link_libraries(ut_Y-Z_crack_jacobian_np4  PdNeigh QuickGrid PdVTK Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_Y-Z_crack_jacobian_np4 python run_unit_test.py mpiexec -np 4 ./ut_Y-Z_crack_jacobian_np4)

#
add_executable(ut_axisAlignedBoundaryPoints ut_axisAlignedBoundaryPoints.cxx)
target_link_libraries(ut_axisAlignedBoundaryPoints  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_axisAlignedBoundaryPoints python run_unit_test.py ./ut_axisAlignedBoundaryPoints)

#
add_executable(ut_NeighborhoodProcessorFrame ut_NeighborhoodProcessorFrame.cxx)
target_link_libraries(ut_NeighborhoodProcessorFrame  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_NeighborhoodProcessorFrame python run_unit_test.py ./ut_NeighborhoodProcessorFrame)

#
add_executable(ut_QuickGrid_loadBal_np2_3x1x1 ut_QuickGrid_loadBal_np2_3x1x1.cxx)
target_link_libraries(ut_QuickGrid_loadBal_np2_3x1x1  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGrid_loadBal_np2_3x1x1 python run_unit_test.py mpiexec -np 2 ./ut_QuickGrid_loadBal_np2_3x1x1)

#
add_executable(ut_QuickGrid_loadBal_np2_4x4x4 ut_QuickGrid_loadBal_np2_4x4x4.cxx)
target_link_libraries(ut_QuickGrid_loadBal_np2_4x4x4  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGrid_loadBal_np2_4x4x4 python run_unit_test.py mpiexec -np 2 ./ut_QuickGrid_loadBal_np2_4x4x4)

#
add_executable(ut_QuickGrid_loadBal_np8_4x4x4 ut_QuickGrid_loadBal_np8_4x4x4.cxx)
target_link_libraries(ut_QuickGrid_loadBal_np8_4x4x4  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_QuickGrid_loadBal_np8_4x4x4 python run_unit_test.py mpiexec -np 8 ./ut_QuickGrid_loadBal_np8_4x4x4)

#
add_executable(utPdVtkKdTree utPdVtkKdTree.cxx)
target_link_libraries(utPdVtkKdTree  QuickGrid PdVTK Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdVtkKdTree python run_unit_test.py ./utPdVtkKdTree)

#
add_executable(utPdVtkPlane utPdVtkPlane.cxx)
target_link_libraries(utPdVtkPlane  PdNeigh Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdVtkPlane python run_unit_test.py ./utPdVtkPlane)

add_executable(utReloadBalance_np4 utReloadBalance_np4.cxx)
target_link_libraries(utReloadBalance_np4  PdNeigh QuickGrid Utilities PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utReloadBalance_np4 python run_unit_test.py mpiexec -np 4 ./utReloadBalance_np4)

add_executable(ut_twoPointReloadBalance ut_twoPointReloadBalance.cxx)
target_link_libraries(ut_twoPointReloadBalance  PdNeigh QuickGrid Utilities PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_twoPointReloadBalance python run_unit_test.py ./ut_twoPointReloadBalance)


add_executable(utZoltanQuery_pointSizeInBytes utZoltanQuery_pointSizeInBytes.cxx)
target_link_libraries(utZoltanQuery_pointSizeInBytes  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utZoltanQuery_pointSizeInBytes python run_unit_test.py ./utZoltanQuery_pointSizeInBytes)

add_executable(utZoltanQuery_packUnpackPointsMultiFunction utZoltanQuery_packUnpackPointsMultiFunction.cxx)
target_link_libraries(utZoltanQuery_packUnpackPointsMultiFunction PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utZoltanQuery_packUnpackPointsMultiFunction python run_unit_test.py ./utZoltanQuery_packUnpackPointsMultiFunction)




