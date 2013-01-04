add_executable(utPeridigm_ElasticPlasticMaterial
               ./utPeridigm_ElasticPlasticMaterial.cpp)

target_link_libraries(utPeridigm_ElasticPlasticMaterial
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  QuickGrid 
  PdNeigh 
  PdVTK 
  Utilities 
  PdMaterialUtilities 
  ${VTK_LIBS}
  ${PARSER_LIBS}
  ${UT_REQUIRED_LIBS}
  ${TRILINOS_LIBS} 
)

add_test (utPeridigm_ElasticPlasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ElasticPlasticMaterial)
