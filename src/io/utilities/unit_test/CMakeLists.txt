#
# unit test for 'utilities' library
#


# include this path
include_directories(${MESH_INPUT_DIR} ${PDNEIGH_DIR} ${UTILITIES_DIR} ${TRILINOS_DIR}/include)

add_executable(ut_Array ut_Array.cxx)
target_link_libraries(ut_Array  ${UT_REQUIRED_LIBS})
add_test (ut_Array python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_Array)

add_executable(ut_Sortable ut_Sortable.cxx)
target_link_libraries(ut_Sortable  ${UT_REQUIRED_LIBS})
add_test (ut_Sortable python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_Sortable)

add_executable(ut_Bits ut_Bits.cxx)
target_link_libraries(ut_Bits  ${UT_REQUIRED_LIBS})

