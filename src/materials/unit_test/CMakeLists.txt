add_subdirectory (twoPoint_SLS_Relaxation/)
add_subdirectory (utPeridigm_IsotropicElasticPlastic/)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unitTest/twoPoint_SLS_Relaxation)
   set(UT_DIRS ${UT_DIRS}  -d
src/materials/unitTest/utPeridigm_IsotropicElasticPlastic)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s
${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})


add_executable(utPeridigm_LinearElasticIsotropicMaterial
               ${CORE_DIR}/Peridigm_Timer.cpp
               ${CORE_DIR}/Peridigm_State.cpp
               ${CORE_DIR}/Peridigm_DataManager.cpp
               ${CORE_DIR}/Peridigm_SerialMatrix.cpp
               ${MESH_OUTPUT_DIR}/Field.cxx
               ${MATERIALS_DIR}/Peridigm_Material.cpp
               ${MATERIALS_DIR}/Peridigm_LinearElasticIsotropicMaterial.cpp
               ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
               ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
               ./utPeridigm_LinearElasticIsotropicMaterial.cpp
)
target_link_libraries(utPeridigm_LinearElasticIsotropicMaterial
  ${Trilinos_LIBRARIES}
  ${MESH_OUTPUT_LIBS}
  ${PdMaterialUtilitiesLib}
  ${REQUIRED_LIBS}
  ${Boost_LIBRARIES}
)
add_test (utPeridigm_LinearElasticIsotropicMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_LinearElasticIsotropicMaterial)
