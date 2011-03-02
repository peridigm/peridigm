#
# Unit tests for PdITI -- Peridynamics Implicit
#
#
# Create symbolic links to files in Test directory
#
#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for units tests for out-of-source build\n")
#   file (GLOB testFiles "${CMAKE_SOURCE_DIR}/operators/unitTest/eigen1D/*.*")
#   foreach (testFile ${testFiles})
#      get_filename_component(fileName ${testFile} NAME)
#	  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${testFile} ${CMAKE_BINARY_DIR}/operators/unitTest/eigen1D/${fileName})  
#   endforeach (testFile)
#endif (NOT ${buildInSource})
#


#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

# Add test
add_executable(utPdImp_eigen1D utPdImp_eigen1D.cxx)
target_link_libraries(utPdImp_eigen1D PdITIOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (utPdImp_eigen1D python run_unit_test.py ./utPdImp_eigen1D)

