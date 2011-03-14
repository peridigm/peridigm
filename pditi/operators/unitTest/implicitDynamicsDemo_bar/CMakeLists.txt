#
# Unit tests for PdITI -- Peridynamics Implicit
#
#
# Create symbolic links to files in Test directory
#
#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for units tests for out-of-source build\n")
#   file (GLOB testFiles "${CMAKE_SOURCE_DIR}/operators/unitTest/implicitDynamicsDemo_bar/*.*")
#   foreach (testFile ${testFiles})
#      get_filename_component(fileName ${testFile} NAME)
#	  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${testFile} ${CMAKE_BINARY_DIR}/operators/unitTest/implicitDynamicsDemo_bar/${fileName})  
#   endforeach (testFile)
#endif (NOT ${buildInSource})
#

# Add test
include_directories(${UTILITIES_DIR} ${PDNEIGH_DIR} ${TRILINOS_DIR}/include )

add_executable(implicitLinearDynamicsDemo_bar implicitLinearDynamicsDemo_bar.cxx)
target_link_libraries(implicitLinearDynamicsDemo_bar PdITIOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (implicitLinearDynamicsDemo_bar python run_unit_test.py ./implicitLinearDynamicsDemo_bar)


