#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#
# Create symbolic links to files in Test directory
#
#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for units tests for out-of-source build\n")
#   file (GLOB testFiles "${CMAKE_SOURCE_DIR}/operators/unitTest/*.*")
#   foreach (testFile ${testFiles})
#      get_filename_component(fileName ${testFile} NAME)
#	  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${testFile} ${CMAKE_BINARY_DIR}/operators/unitTest/${fileName})  
#   endforeach (testFile)
#endif (NOT ${buildInSource})
#

# Add test
add_executable(utPdImpMatrix_np2_11x1x1 utPdImpMatrix_np2_11x1x1.cxx)
target_link_libraries(utPdImpMatrix_np2_11x1x1 PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImpMatrix_np2_11x1x1 python run_unit_test.py mpiexec -np 2 utPdImpMatrix_np2_11x1x1)

# Add test
add_executable(utPdImpMatrix_11x1x1 utPdImpMatrix_11x1x1.cxx)
target_link_libraries(utPdImpMatrix_11x1x1 PdImpOperators ${PERIDIGM_LIBS}  ${TRILINOS_LIBS} ${VTK_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImpMatrix_11x1x1 python run_unit_test.py ./utPdImpMatrix_11x1x1)

# Add test
add_executable(utPdImp_twoPointJacobian utPdImp_twoPointJacobian.cxx)
target_link_libraries(utPdImp_twoPointJacobian PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_twoPointJacobian python run_unit_test.py ./utPdImp_twoPointJacobian)

# Add test
add_executable(utPdImp_twoPointJacobian_np2 utPdImp_twoPointJacobian_np2.cxx)
target_link_libraries(utPdImp_twoPointJacobian_np2 PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_twoPointJacobian_np2 python run_unit_test.py mpiexec -np 2 ./utPdImp_twoPointJacobian_np2)

# Add test
#PdZoltan PdQuickGrid PdNeighborhood  PdVTK
add_executable(utPdImp_probeJacobianNu.0.0 utPdImp_probeJacobianNu.0.0.cxx)
target_link_libraries(utPdImp_probeJacobianNu.0.0 PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_probeJacobianNu.0.0 python run_unit_test.py ./utPdImp_probeJacobianNu.0.0)

# Add test
add_executable(utPdImp_linearSolve_barGravity_npX utPdImp_linearSolve_barGravity_npX.cxx)
target_link_libraries(utPdImp_linearSolve_barGravity_npX PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_linearSolve_barGravity_npX python run_unit_test.py mpiexec -np 5 ./utPdImp_linearSolve_barGravity_npX)

# Add test
add_executable(utPdImp_graph2x2x1_np4 utPdImp_graph2x2x1_np4.cxx)
target_link_libraries(utPdImp_graph2x2x1_np4 PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_graph2x2x1_np4 python run_unit_test.py mpiexec -np 4 ./utPdImp_graph2x2x1_np4)

# Add test
add_executable(utPdImp_componentDirichletBcSpec utPdImp_componentDirichletBcSpec.cxx)
target_link_libraries(utPdImp_componentDirichletBcSpec PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_componentDirichletBcSpec python run_unit_test.py ./utPdImp_componentDirichletBcSpec)

# Add test
add_executable(utPdImp_constantBodyLoad utPdImp_constantBodyLoad.cxx)
target_link_libraries(utPdImp_constantBodyLoad PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_constantBodyLoad python run_unit_test.py ./utPdImp_constantBodyLoad)

# Add test
add_executable(utPdImp_pullCylinder_npX utPdImp_pullCylinder_npX.cxx)
target_link_libraries(utPdImp_pullCylinder_npX PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (utPdImp_pullCylinder_npX python run_unit_test.py mpiexec -np 4 ./utPdImp_pullCylinder_npX)

# Add test
add_executable(utPdImp_smallCylinderPull_npX utPdImp_smallCylinderPull_npX.cxx)
target_link_libraries(utPdImp_smallCylinderPull_npX PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (utPdImp_smallCylinderPull_npX python run_unit_test.py mpiexec -np 4 ./utPdImp_smallCylinderPull_npX)

# Add test
add_executable(utPdImp_isotropicElasticPlasticSpec utPdImp_isotropicElasticPlasticSpec.cxx)
target_link_libraries(utPdImp_isotropicElasticPlasticSpec PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_isotropicElasticPlasticSpec python run_unit_test.py ./utPdImp_isotropicElasticPlasticSpec)

# Add test
add_executable(utPdImp_isotropicElasticPlastic_twoPointProbe utPdImp_isotropicElasticPlastic_twoPointProbe.cxx)
target_link_libraries(utPdImp_isotropicElasticPlastic_twoPointProbe PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_isotropicElasticPlastic_twoPointProbe python run_unit_test.py ./utPdImp_isotropicElasticPlastic_twoPointProbe)

# Add test
add_executable(utPdImp_implicitLinearDynamicsDemo_twoPoint utPdImp_implicitLinearDynamicsDemo_twoPoint.cxx)
target_link_libraries(utPdImp_implicitLinearDynamicsDemo_twoPoint PdImpOperators ${PERIDIGM_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (utPdImp_implicitLinearDynamicsDemo_twoPoint python run_unit_test.py ./utPdImp_implicitLinearDynamicsDemo_twoPoint)

#add_custom_target( utest
#   COMMAND ctest
#)

