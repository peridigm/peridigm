
#
add_executable(ut_frameset_2x2x1_np4 ut_frameset_2x2x1_np4.cxx)
target_link_libraries(ut_frameset_2x2x1_np4  PdNeigh ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_frameset_2x2x1_np4 python run_unit_test.py mpiexec -np 4 ./ut_frameset_2x2x1_np4)

#
add_executable(ut_neighborhood_list ut_neighborhood_list.cxx)
target_link_libraries(ut_neighborhood_list  PdNeigh ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_neighborhood_list python run_unit_test.py ./ut_neighborhood_list)
