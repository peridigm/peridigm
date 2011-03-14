#
# unit test for 'utilities' library
#


# include this path
include_directories(${UTILITIES_DIR})


add_executable(ut_Array ut_Array.cxx)
target_link_libraries(ut_Array  ${UT_REQUIRED_LIBS})
add_test (ut_Array python run_unit_test.py ./ut_Array)

add_executable(ut_Sortable ut_Sortable.cxx)
target_link_libraries(ut_Sortable  ${UT_REQUIRED_LIBS})
add_test (ut_Sortable python run_unit_test.py ./ut_Sortable)

