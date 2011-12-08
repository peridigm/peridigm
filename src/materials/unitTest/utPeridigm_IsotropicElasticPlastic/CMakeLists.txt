add_executable(utPeridigm_IsotropicElasticPlasticMaterial
               ../../Peridigm_Material.cpp
               ../../Peridigm_IsotropicElasticPlasticMaterial.cpp
               ../../Peridigm_DamageModelFactory.cpp
               ../../Peridigm_CriticalStretchDamageModel.cpp
               ../../../core/Peridigm_Timer.cpp
               ../../../core/Peridigm_State.cpp
               ../../../core/Peridigm_DataManager.cpp
               ../../../core/Peridigm_SerialMatrix.cpp
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

