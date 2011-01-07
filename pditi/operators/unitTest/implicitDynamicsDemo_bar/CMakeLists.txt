#
# Unit tests for PdImp -- Peridynamics Implicit
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
#
# Collection of libs that we need to link with
#include_directories(${TRILINOS_DIR}/include  ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${PERIDIGM_IO_SRC} ../.. ${PERIDIGM_MAT_UTILS_SRC})
#link_directories(${TRILINOS_DIR}/lib ${VTK_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${PERIDIGM_IO_LIB} ../.. ${PERIDIGM_MAT_UTILS_LIB})
#
#
# Add test
#add_executable(implicitLinearDynamicsDemo_bar implicitLinearDynamicsDemo_bar.cxx)
#target_link_libraries(implicitLinearDynamicsDemo_bar PdImpOperators ${PdMaterialUtilities_LIBRARY} PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#add_test (implicitLinearDynamicsDemo_bar python run_unit_test.py ./implicitLinearDynamicsDemo_bar)
#
#
#
#set (REQUIRED_LIBS
#liblapack.a
#libblas.a
#libmpich.a
#${PTHREAD_LIB}
#libgfortran.a
#libm.a
#)


#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

# Add test
add_executable(implicitLinearDynamicsDemo_bar implicitLinearDynamicsDemo_bar.cxx)
target_link_libraries(implicitLinearDynamicsDemo_bar PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (implicitLinearDynamicsDemo_bar python run_unit_test.py ./implicitLinearDynamicsDemo_bar)


