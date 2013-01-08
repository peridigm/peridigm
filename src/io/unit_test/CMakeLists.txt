
add_executable(utPeridigm_SearchTree ./utPeridigm_SearchTree.cpp)
target_link_libraries(utPeridigm_SearchTree
  ${Peridigm_LIBRARY}
  ${Zoltan_LIBRARY}
  ${REQUIRED_LIBS}
  ${Boost_LIBRARIES}
)
add_test (utPeridigm_SearchTree python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_SearchTree)
