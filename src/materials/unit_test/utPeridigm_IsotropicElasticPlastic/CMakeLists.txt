add_executable(utPeridigm_IsotropicHardeningPlasticMaterial
               ${MATERIALS_DIR}/Peridigm_Material.cpp
               ${MATERIALS_DIR}/Peridigm_IsotropicHardeningPlasticMaterial.cpp
               ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
               ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
               ${CORE_DIR}/Peridigm_Timer.cpp
               ${CORE_DIR}/Peridigm_State.cpp
               ${CORE_DIR}/Peridigm_DataManager.cpp
               ${CORE_DIR}/Peridigm_SerialMatrix.cpp
               ./utPeridigm_IsotropicHardeningPlasticMaterial.cpp)

target_link_libraries(utPeridigm_IsotropicHardeningPlasticMaterial
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

add_test (utPeridigm_IsotropicHardeningPlasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_IsotropicHardeningPlasticMaterial)
