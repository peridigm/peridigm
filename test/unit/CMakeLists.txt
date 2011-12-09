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
add_test (utPeridigm_PdQuickGridDiscretization          python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ../../src/utPeridigm_PdQuickGridDiscretization)
add_test (utPeridigm_PdQuickGridDiscretization_MPI_np2  python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ../../src/utPeridigm_PdQuickGridDiscretization_MPI_np2)
add_test (utPeridigm_Compute_Force                      python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ../../src/utPeridigm_Compute_Force)
add_test (utPeridigm_Compute_Force_MPI_np2              python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ../../src/utPeridigm_Compute_Force)

add_custom_target( utest
   COMMAND ctest
)
