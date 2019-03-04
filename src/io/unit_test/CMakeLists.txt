
# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/io/unit_test/input_files)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS})
endif (NOT ${buildInSource})

add_executable(utPeridigm_ProximitySearch ./utPeridigm_ProximitySearch.cpp)
target_link_libraries(utPeridigm_ProximitySearch
  ${Peridigm_LIBRARY}
  ${PDNEIGH_LIBS}
  ${MESH_INPUT_LIBS}
  ${UTILITIES_LIBS}
  ${Trilinos_LIBRARIES}
  ${REQUIRED_LIBS}
)
add_test (utPeridigm_ProximitySearch_np1 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ProximitySearch)
add_test (utPeridigm_ProximitySearch_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_ProximitySearch)
add_test (utPeridigm_ProximitySearch_np3 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 3 ./utPeridigm_ProximitySearch)
add_test (utPeridigm_ProximitySearch_np4 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 4 ./utPeridigm_ProximitySearch)
add_test (utPeridigm_ProximitySearch_np5 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 5 ./utPeridigm_ProximitySearch)

add_executable(utPeridigm_SearchTree ./utPeridigm_SearchTree.cpp)
target_link_libraries(utPeridigm_SearchTree
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${REQUIRED_LIBS}
)
add_test (utPeridigm_SearchTree python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_SearchTree)

# This test is not run as part of the normal test suite
add_executable(utPeridigm_SearchTree_Performance ./utPeridigm_SearchTree_Performance.cpp)
target_link_libraries(utPeridigm_SearchTree_Performance
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${REQUIRED_LIBS}
)

add_executable(ut_kdtree_nn_search ut_kdtree_nn_search.cpp)
target_link_libraries(ut_kdtree_nn_search QuickGrid
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${REQUIRED_LIBS}
)
add_test (ut_kdtree_nn_search python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_kdtree_nn_search )

