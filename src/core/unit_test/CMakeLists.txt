#
# Unit Tests
#

add_executable(utPeridigm utPeridigm.cpp)
target_link_libraries(utPeridigm ${Peridigm_LIBRARY} ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})
add_test (utPeridigm python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm)

add_executable(utPeridigm_MPI_np2 utPeridigm_MPI_np2.cpp)
target_link_libraries(utPeridigm_MPI_np2 ${Peridigm_LIBRARY} ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})
add_test (utPeridigm_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_MPI_np2)

add_executable(utPeridigm_State ./utPeridigm_State.cpp)
target_link_libraries(utPeridigm_State ${Peridigm_LIBRARY} ${Trilinos_LIBRARIES} ${REQUIRED_LIBS} ${Boost_LIBRARIES})
add_test (utPeridigm_State python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_State)
add_test (utPeridigm_State_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_State)

