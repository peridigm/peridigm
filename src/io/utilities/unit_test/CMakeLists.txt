#
# unit test for 'utilities' library
#

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d utilities/unit_test)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${IO_DIR} -b ${CMAKE_BINARY_DIR}/src/io ${UT_DIRS} ) 
endif (NOT ${buildInSource})

add_executable(ut_Array ut_Array.cxx)
target_link_libraries(ut_Array ${Peridigm_LIBRARY} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
add_test (ut_Array python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_Array)

add_executable(ut_Sortable ut_Sortable.cxx)
target_link_libraries(ut_Sortable ${Peridigm_LIBRARY} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
add_test (ut_Sortable python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_Sortable)

add_executable(ut_Bits ut_Bits.cxx)
target_link_libraries(ut_Bits  ${Peridigm_LIBRARY} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})

