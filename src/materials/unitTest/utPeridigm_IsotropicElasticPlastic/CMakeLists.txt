add_executable(utPeridigm_IsotropicElasticPlasticMaterial
		../../Peridigm_Material.cpp
		../../Peridigm_IsotropicElasticPlasticMaterial.cpp
        ../../Peridigm_CriticalStretchDamageModel.cpp
        ../../../Peridigm_State.cpp
        ../../../Peridigm_DataManager.cpp
        ../../../Peridigm_SerialMatrix.cpp
		./utPeridigm_IsotropicElasticPlasticMaterial.cpp)
target_link_libraries(utPeridigm_IsotropicElasticPlasticMaterial
${Trilinos_LIBRARIES}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_test (utPeridigm_IsotropicElasticPlasticMaterial python run_unit_test.py ./utPeridigm_IsotropicElasticPlasticMaterial)

