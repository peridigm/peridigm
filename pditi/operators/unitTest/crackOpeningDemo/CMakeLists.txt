#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

include_directories(${PD_MATERIAL_UTILITIES_DIR} ${PDNEIGH_DIR} ${MESH_INPUT_DIR} ${MESH_OUTPUT_DIR} ${UTILITIES_DIR} ${TRILINOS_DIR}/include)

# Add test
add_executable(utCrackOpeningDemo_npX utCrackOpeningDemo_npX.cxx)
target_link_libraries(utCrackOpeningDemo_npX utPdITI PdITIOperators QuickGrid PdNeigh PdVTK Utilities PdMaterialUtilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
