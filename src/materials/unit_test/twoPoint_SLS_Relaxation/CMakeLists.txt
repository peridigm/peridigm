add_executable(twoPoint_SLS_Relaxation
               ./twoPoint_SLS_Relaxation.cpp)

target_link_libraries(twoPoint_SLS_Relaxation
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  QuickGrid 
  PdNeigh 
  PdVTK 
  Utilities 
  PdMaterialUtilities 
  ${VTK_LIBS} 
  ${PARSER_LIBS}
  ${UT_REQUIRED_LIBS}
  ${TRILINOS_LIBS} 
)

add_test (twoPoint_SLS_Relaxation python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./twoPoint_SLS_Relaxation)

