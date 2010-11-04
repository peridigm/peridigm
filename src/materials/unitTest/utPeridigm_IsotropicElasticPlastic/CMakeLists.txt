add_executable(utPeridigm_IsotropicElasticPlasticMaterial
		../../Peridigm_IsotropicElasticPlasticMaterial.cpp
		../../Peridigm_DecompositionStates.cpp
        ../../Peridigm_CriticalStretchDamageModel.cpp
        ../../../Peridigm_DataManager.cpp
		./utPeridigm_IsotropicElasticPlasticMaterial.cpp)
target_link_libraries(utPeridigm_IsotropicElasticPlasticMaterial
${TRILINOS_LIBS}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_test (utPeridigm_IsotropicElasticPlasticMaterial python run_unit_test.py ./utPeridigm_IsotropicElasticPlasticMaterial)

