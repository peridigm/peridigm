#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

include_directories(${UTILITIES_DIR} ${PDNEIGH_DIR} ${TRILINOS_DIR}/include )

# Add test
add_executable(utCrackOpeningDemo_npX utCrackOpeningDemo_npX.cxx)
target_link_libraries(utCrackOpeningDemo_npX utPdITI PdITIOperators PdNeigh ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
