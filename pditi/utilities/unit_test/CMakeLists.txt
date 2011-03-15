#
# unit test for 'utilities' library
#


# include this path
include_directories(${MESH_INPUT_DIR} ${PDNEIGH_DIR} ${UTILITIES_DIR} ${TRILINOS_DIR}/include)

add_executable(ut_Array ut_Array.cxx)
target_link_libraries(ut_Array  ${UT_REQUIRED_LIBS})
add_test (ut_Array python run_unit_test.py ./ut_Array)

add_executable(ut_Sortable ut_Sortable.cxx)
target_link_libraries(ut_Sortable  ${UT_REQUIRED_LIBS})
add_test (ut_Sortable python run_unit_test.py ./ut_Sortable)

add_executable(ut_axisAlignedBoundaryPoints ut_axisAlignedBoundaryPoints.cxx)
target_link_libraries(ut_axisAlignedBoundaryPoints  PdNeigh QuickGrid Utilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_axisAlignedBoundaryPoints python run_unit_test.py ./ut_axisAlignedBoundaryPoints)
