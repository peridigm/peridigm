add_executable(utPeridigm_IsotropicElasticPlasticMaterial
               ${MATERIALS_DIR}/Peridigm_Material.cpp
               ${MATERIALS_DIR}/Peridigm_IsotropicElasticPlasticMaterial.cpp
               ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
               ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
               ${CORE_DIR}/Peridigm_Timer.cpp
               ${CORE_DIR}/Peridigm_State.cpp
               ${CORE_DIR}/Peridigm_DataManager.cpp
               ${CORE_DIR}/Peridigm_SerialMatrix.cpp
               ./utPeridigm_IsotropicElasticPlasticMaterial.cpp)

target_link_libraries(utPeridigm_IsotropicElasticPlasticMaterial
  ${Trilinos_LIBRARIES}
  QuickGrid 
  PdNeigh 
  PdVTK 
  Utilities 
  PdMaterialUtilities 
  ${VTK_LIBS} 
  ${UT_REQUIRED_LIBS}
  ${TRILINOS_LIBS} 
)

add_test (utPeridigm_IsotropicElasticPlasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_IsotropicElasticPlasticMaterial)
