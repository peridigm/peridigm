add_executable(twoPoint_SLS_Relaxation
               ./twoPoint_SLS_Relaxation.cpp)

set_property(TARGET twoPoint_SLS_Relaxation
             PROPERTY IMPORTED_LINK_INTERFACE_MULTIPLICITY 4
)

target_link_libraries(twoPoint_SLS_Relaxation
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
  QuickGrid 
  PdNeigh 
  Utilities 
  PdMaterialUtilities
  PdField
  ${PARSER_LIBS}
  ${UT_REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)

add_test (twoPoint_SLS_Relaxation python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./twoPoint_SLS_Relaxation)

