#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

# Add test
add_executable(utCrackOpeningDemo utCrackOpeningDemo.cxx)
target_link_libraries(utCrackOpeningDemo utPdITI PdITIOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
