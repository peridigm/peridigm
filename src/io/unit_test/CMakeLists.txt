
# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/io/unit_test/input_files)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS})
endif (NOT ${buildInSource})

add_executable(utPeridigm_SearchTree ./utPeridigm_SearchTree.cpp)
target_link_libraries(utPeridigm_SearchTree
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${Zoltan_LIBRARY}
  ${REQUIRED_LIBS}
  ${Boost_LIBRARIES}
)
add_test (utPeridigm_SearchTree python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_SearchTree)

# This test is not run as part of the normal test suite
add_executable(utPeridigm_SearchTree_Performance ./utPeridigm_SearchTree_Performance.cpp)
target_link_libraries(utPeridigm_SearchTree_Performance
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${Zoltan_LIBRARY}
  ${REQUIRED_LIBS}
  ${Boost_LIBRARIES}
)
