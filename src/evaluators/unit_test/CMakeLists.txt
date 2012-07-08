#add_executable(utPHAL_UpdateForceState ./utPHAL_UpdateForceState.cpp)
#target_link_libraries(utPHAL_UpdateForceState
#                      ${Peridigm_LIBRARY}
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)
#add_test (utPHAL_UpdateForceState python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPHAL_UpdateForceState)

###add_executable(utPHAL_UpdateForceState_MPI_np2 ./utPHAL_UpdateForceState_MPI_np2.cpp)
###target_link_libraries(utPHAL_UpdateForceState_MPI_np2
###                      ${Peridigm_LIBRARY}
###                      ${Trilinos_LIBRARIES}
###                      ${PdMaterialUtilitiesLib}
###                      ${REQUIRED_LIBS}
###                      ${Boost_LIBRARIES}
###)
###add_test (utPHAL_UpdateForceState_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPHAL_UpdateForceState_MPI_np2)

#add_executable(utPHAL_EvaluateForce ./utPHAL_EvaluateForce.cpp)
#target_link_libraries(utPHAL_EvaluateForce
#                      ${Peridigm_LIBRARY}
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)
#add_test (utPHAL_EvaluateForce python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPHAL_EvaluateForce)

###add_executable(utPHAL_EvaluateForce_MPI_np2 ./utPHAL_EvaluateForce_MPI_np2.cpp)
###target_link_libraries(utPHAL_EvaluateForce_MPI_np2
###                      ${Peridigm_LIBRARY}
###                      ${Trilinos_LIBRARIES}
###                      ${PdMaterialUtilitiesLib}
###                      ${REQUIRED_LIBS}
###                      ${Boost_LIBRARIES}
###)
###add_test (utPHAL_EvaluateForce_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPHAL_EvaluateForce_MPI_np2)


