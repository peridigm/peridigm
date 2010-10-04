# List tests to be run by "make test"
#
# Create symbolic links to files in Test directory
#

### MESSAGE TO DAVE

# I replaced the cmake script below with a python script -- its already hooked up and running.
# Its should work out of the box when you turn on these tests.  
#



#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for regression tests for out-of-source build\n")
#   file (GLOB testFiles "${CMAKE_SOURCE_DIR}/test/regression/*.*")
#   foreach (testFile ${testFiles})
#      get_filename_component(fileName ${testFile} NAME)
#	  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${testFile} ${CMAKE_BINARY_DIR}/test/regression/${fileName})  
#   endforeach (testFile)
#endif (NOT ${buildInSource})



#
# Add the tests
#
##add_test (no_args_serial python no_args_serial.py)
##add_test (no_args_parallel python no_args_parallel.py)
##add_test (PdQuickGrid_3x3_serial python PdQuickGrid_3x3_serial.py)
##add_test (PdQuickGrid_3x3_mpi2 python PdQuickGrid_3x3_mpi2.py)

add_custom_target( rtest
   COMMAND ctest
)
