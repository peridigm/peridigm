add_executable(utPeridigm_IsotropicElasticPlasticMaterial
               ./utPeridigm_IsotropicElasticPlasticMaterial.cpp)

target_link_libraries(utPeridigm_IsotropicElasticPlasticMaterial
  ${Peridigm_LIBRARY}
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
