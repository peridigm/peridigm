add_executable(utPHAL_UpdateForceState
               ${CORE_DIR}/PHAL_PeridigmTraits.cpp
               ${CORE_DIR}/PHAL_Dimension.cpp
               ${CORE_DIR}/Peridigm_Timer.cpp
               ${CORE_DIR}/Peridigm_Block.cpp
               ${CORE_DIR}/Peridigm_State.cpp
               ${CORE_DIR}/Peridigm_DataManager.cpp
               ${CORE_DIR}/Peridigm_SerialMatrix.cpp
               ${MATERIALS_DIR}/Peridigm_Material.cpp
               ${MATERIALS_DIR}/Peridigm_LinearElasticIsotropicMaterial.cpp
               ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
               ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
               ${MESH_OUTPUT_DIR}/Field.cxx
               ./utPHAL_UpdateForceState.cpp)
target_link_libraries(utPHAL_UpdateForceState
                      ${Trilinos_LIBRARIES}
                      ${PdMaterialUtilitiesLib}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)
add_test (utPHAL_UpdateForceState python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPHAL_UpdateForceState)

#add_executable(utPHAL_UpdateForceState_MPI_np2
#               ${CORE_DIR}/PHAL_PeridigmTraits.cpp
#               ${CORE_DIR}/PHAL_Dimension.cpp
#               ${MATERIALS_DIR}/Peridigm_LinearElasticIsotropicMaterial.cpp
#               ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
#               ./utPHAL_UpdateForceState_MPI_np2.cpp)
#target_link_libraries(utPHAL_UpdateForceState_MPI_np2
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)
#add_test (utPHAL_UpdateForceState_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPHAL_UpdateForceState_MPI_np2)

add_executable(utPHAL_EvaluateForce
               ${CORE_DIR}/PHAL_PeridigmTraits.cpp
               ${CORE_DIR}/PHAL_Dimension.cpp
               ${CORE_DIR}/Peridigm_Timer.cpp
               ${CORE_DIR}/Peridigm_State.cpp
               ${CORE_DIR}/Peridigm_DataManager.cpp
               ${CORE_DIR}/Peridigm_Block.cpp
               ${CORE_DIR}/Peridigm_SerialMatrix.cpp
               ${MESH_OUTPUT_DIR}/Field.cxx
               ${MATERIALS_DIR}/Peridigm_Material.cpp
               ${MATERIALS_DIR}/Peridigm_LinearElasticIsotropicMaterial.cpp
               ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
               ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
               ./utPHAL_EvaluateForce.cpp
)
target_link_libraries(utPHAL_EvaluateForce
                      ${Trilinos_LIBRARIES}
                      ${PdMaterialUtilitiesLib}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)
add_test (utPHAL_EvaluateForce python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPHAL_EvaluateForce)

#add_executable(utPHAL_EvaluateForce_MPI_np2
#               PHAL_PeridigmTraits.cpp
#               PHAL_Dimension.cpp
#               ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
#               ./materials/Peridigm_CriticalStretchDamageModel.cpp
#               ./evaluators/utPHAL_EvaluateForce_MPI_np2.cpp)
#target_link_libraries(utPHAL_EvaluateForce_MPI_np2
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)
#add_test (utPHAL_EvaluateForce_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPHAL_EvaluateForce_MPI_np2)


