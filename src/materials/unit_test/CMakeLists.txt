add_subdirectory (twoPoint_SLS_Relaxation/)
add_subdirectory (utPeridigm_ElasticPlastic/)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unit_test/twoPoint_SLS_Relaxation)
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unit_test/utPeridigm_ElasticPlastic)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})

add_executable(utPeridigm_ElasticMaterial ./utPeridigm_ElasticMaterial.cpp
#is the ) correct here?
)
target_link_libraries(utPeridigm_ElasticMaterial
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${PdMaterialUtilitiesLib}
  PdField
  ${PARSER_LIBS}
  ${REQUIRED_LIBS}
  ${Boost_LIBRARIES}
)
add_test (utPeridigm_ElasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ElasticMaterial)


add_executable(utPeridigm_MultiphysicsElasticMaterial ./utPeridigm_MultiphysicsElasticMaterial.cpp)
target_link_libraries(utPeridigm_MultiphysicsElasticMaterial
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  ${PdMaterialUtilitiesLib}
  PdField
  ${PARSER_LIBS}
  ${REQUIRED_LIBS}
  ${Boost_LIBRARIES}
)
add_test (utPeridigm_MultiphysicsElasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_MultiphysicsElasticMaterial)
