#
# CMakeLists.txt file for Peridigm Unit Tests
#
# Create symbolic links to files in Test directory
#
#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for units tests for out-of-source build\n")
#   file (GLOB testFiles "${CMAKE_SOURCE_DIR}/test/unit/*.*")
#   foreach (testFile ${testFiles})
#      get_filename_component(fileName ${testFile} NAME)
#	  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${testFile} ${CMAKE_BINARY_DIR}/test/unit/${fileName})  
#   endforeach (testFile)
#endif (NOT ${buildInSource})
#
# Add the tests
# 
# Peridigm unit tests:
add_test (utPeridigm_PdQuickGridDiscretization          python run_unit_test.py ../../src/utPeridigm_PdQuickGridDiscretization)
add_test (utPeridigm_PdQuickGridDiscretization_MPI_np2  python run_unit_test.py mpiexec -np 2 ../../src/utPeridigm_PdQuickGridDiscretization_MPI_np2)
add_test (utPeridigm_LinearElasticIsotropicMaterial     python run_unit_test.py ../../src/utPeridigm_LinearElasticIsotropicMaterial)
add_test (utPHAL_UpdateForceState                       python run_unit_test.py ../../src/utPHAL_UpdateForceState)
#add_test (utPHAL_UpdateForceState_MPI_np2              python run_unit_test.py mpiexec -np 2 ../../src/utPHAL_UpdateForceState_MPI_np2)
add_test (utPHAL_EvaluateForce                          python run_unit_test.py ../../src/utPHAL_EvaluateForce)
#add_test (utPHAL_EvaluateForce_MPI_np2                 python run_unit_test.py mpiexec -np 2 ../../src/utPHAL_EvaluateForce_MPI_np2)
add_test (utPeridigm_decompositionStates                python run_unit_test.py ../../src/utPeridigm_decompositionStates)
add_test (utPeridigm                                    python run_unit_test.py ../../src/utPeridigm)
# io unit tests:
add_test (utPdVtkKdTree                                python run_unit_test.py ../../io/unitTest/utPdVtkKdTree)
add_test (utPdQuickGridHorizon                         python run_unit_test.py ../../io/unitTest/utPdQuickGridHorizon)
add_test (utPdSmallMeshCylinder                        python run_unit_test.py ../../io/unitTest/utPdSmallMeshCylinder)
add_test (utPdSmallMeshCylinder_np4                    python run_unit_test.py mpiexec -np 4 ../../io/unitTest/utPdSmallMeshCylinder_np4)
add_test (utPdQ2CylinderRingHorizon                    python run_unit_test.py ../../io/unitTest/utPdQ2CylinderRingHorizon)
add_test (utPdQuickGrid_loadBal_np2_3x1x1              python run_unit_test.py mpiexec -np 2 ../../io/unitTest/utPdQuickGrid_loadBal_np2_3x1x1)
add_test (utPdQuickGrid_loadBal_np2_4x4x4              python run_unit_test.py mpiexec -np 2 ../../io/unitTest/utPdQuickGrid_loadBal_np2_4x4x4)
add_test (utPdQuickGrid_loadBal_np8_4x4x4              python run_unit_test.py mpiexec -np 8 ../../io/unitTest/utPdQuickGrid_loadBal_np8_4x4x4)
add_test (utPdQuickGridMPI_np3                         python run_unit_test.py mpiexec -np 3 ../../io/unitTest/utPdQuickGridMPI_np3)
add_test (utPdQuickGridMPI_np2_3x3x2                   python run_unit_test.py mpiexec -np 2 ../../io/unitTest/utPdQuickGridMPI_np2_3x3x2)
add_test (utZoltanQuery_pointSizeInBytes               python run_unit_test.py ../../io/unitTest/utZoltanQuery_pointSizeInBytes)
add_test (utZoltanQuery_packUnpackPointsMultiFunction  python run_unit_test.py ../../io/unitTest/utZoltanQuery_packUnpackPointsMultiFunction)
add_test (utPdNeighborhoodSort                         python run_unit_test.py ../../io/unitTest/utPdNeighborhoodSort)
add_test (utPdNeighborhoodProcessorFrame               python run_unit_test.py ../../io/unitTest/utPdNeighborhoodProcessorFrame)
add_test (utPdNeighborhoodList                         python run_unit_test.py ../../io/unitTest/utPdNeighborhoodList)
add_test (utPd_axisAlignedBoundaryPoints               python run_unit_test.py ../../io/unitTest/utPd_axisAlignedBoundaryPoints)
add_test (utPdQuickGrid_solidCylinder_np1              python run_unit_test.py ../../io/unitTest/utPdQuickGrid_solidCylinder_np1)
add_test (utPdQuickGrid_solidCylinder_np2              python run_unit_test.py mpiexec -np 2 ../../io/unitTest/utPdQuickGrid_solidCylinder_np2)
add_test (utReloadBalance_np4                          python run_unit_test.py mpiexec -np 4 ../../io/unitTest/utReloadBalance_np4)
add_test (utPd_frameSet_2x2x1_np4                      python run_unit_test.py mpiexec -np 4 ../../io/unitTest/utPd_frameSet_2x2x1_np4)
add_test (utVTK_ioExample_np4                          python run_unit_test.py mpiexec -np 4 ../../io/unitTest/utVTK_ioExample_np4)

add_custom_target( utest
   COMMAND ctest
)
