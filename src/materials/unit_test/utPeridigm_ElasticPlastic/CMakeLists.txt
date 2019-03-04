add_executable(utPeridigm_ElasticPlasticMaterial
               ./utPeridigm_ElasticPlasticMaterial.cpp)

target_link_libraries(utPeridigm_ElasticPlasticMaterial
  ${Peridigm_LIBRARY}
  QuickGrid
  PdNeigh
  Utilities
  PdMaterialUtilities
  PdField
  ${UT_REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)

add_test (utPeridigm_ElasticPlasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ElasticPlasticMaterial)
