
#add_executable(utPHPD_EvaluateForce ./utPHPD_EvaluateForce.cpp)
#target_link_libraries(utPHPD_EvaluateForce
#                      ${Peridigm_LIBRARY}
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)
#add_test (utPHPD_EvaluateForce python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPHPD_EvaluateForce)

###add_executable(utPHPD_EvaluateForce_MPI_np2 ./utPHPD_EvaluateForce_MPI_np2.cpp)
###target_link_libraries(utPHPD_EvaluateForce_MPI_np2
###                      ${Peridigm_LIBRARY}
###                      ${Trilinos_LIBRARIES}
###                      ${PdMaterialUtilitiesLib}
###                      ${REQUIRED_LIBS}
###                      ${Boost_LIBRARIES}
###)
###add_test (utPHPD_EvaluateForce_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPHPD_EvaluateForce_MPI_np2)


