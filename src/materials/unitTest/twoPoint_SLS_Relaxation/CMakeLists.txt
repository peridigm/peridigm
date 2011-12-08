add_executable(twoPoint_SLS_Relaxation
               ../../Peridigm_Material.cpp
               ../../Peridigm_ViscoelasticStandardLinearSolid.cpp
               ../../Peridigm_DamageModelFactory.cpp
               ../../Peridigm_CriticalStretchDamageModel.cpp
               ../../../core/Peridigm_State.cpp
               ../../../core/Peridigm_DataManager.cpp
               ../../../core/Peridigm_SerialMatrix.cpp
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

add_test (twoPoint_SLS_Relaxation python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./twoPoint_SLS_Relaxation)

