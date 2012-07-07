add_executable(twoPoint_SLS_Relaxation
               ${MATERIALS_DIR}/Peridigm_Material.cpp
               ${MATERIALS_DIR}/Peridigm_ViscoelasticStandardLinearSolid.cpp
               ${CORE_DIR}/Peridigm_State.cpp
               ${CORE_DIR}/Peridigm_DataManager.cpp
               ${CORE_DIR}/Peridigm_SerialMatrix.cpp
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

