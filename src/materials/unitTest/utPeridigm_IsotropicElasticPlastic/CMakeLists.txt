add_executable(utPeridigm_IsotropicElasticPlasticMaterial
		../../Peridigm_Material.cpp
		../../Peridigm_IsotropicElasticPlasticMaterial.cpp
		../../Peridigm_DamageModelFactory.cpp
        ../../Peridigm_CriticalStretchDamageModel.cpp
        ../../../Peridigm_Timer.cpp
        ../../../Peridigm_State.cpp
        ../../../Peridigm_DataManager.cpp
        ../../../Peridigm_SerialMatrix.cpp
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

add_test (utPeridigm_IsotropicElasticPlasticMaterial python run_unit_test.py ./utPeridigm_IsotropicElasticPlasticMaterial)

