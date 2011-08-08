add_executable(utPeridigm_twoPointMaxwellRelaxation
		../../Peridigm_Material.cpp
		../../Peridigm_ViscoelasticStandardLinearSolid.cpp
		../../Peridigm_DamageModelFactory.cpp
        ../../Peridigm_CriticalStretchDamageModel.cpp
        ../../../Peridigm_State.cpp
        ../../../Peridigm_DataManager.cpp
        ../../../Peridigm_SerialMatrix.cpp
		./utPeridigm_twoPointMaxwellRelaxation.cpp)
target_link_libraries(utPeridigm_twoPointMaxwellRelaxation
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

add_test (utPeridigm_twoPointMaxwellRelaxation python run_unit_test.py ./utPeridigm_twoPointMaxwellRelaxation)

