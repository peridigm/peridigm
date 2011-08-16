add_executable(twoPoint_SLS_Relaxation
		../../Peridigm_Material.cpp
		../../Peridigm_ViscoelasticStandardLinearSolid.cpp
		../../Peridigm_DamageModelFactory.cpp
        ../../Peridigm_CriticalStretchDamageModel.cpp
        ../../../Peridigm_State.cpp
        ../../../Peridigm_DataManager.cpp
        ../../../Peridigm_SerialMatrix.cpp
		./twoPoint_SLS_Relaxation.cpp)
target_link_libraries(twoPoint_SLS_Relaxation
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

add_test (twoPoint_SLS_Relaxation python run_unit_test.py ./twoPoint_SLS_Relaxation)

