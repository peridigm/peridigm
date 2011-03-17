#
# Operators for PdITI -- Peridynamics Implicit
#
#add_subdirectory(unitTest)

include_directories(${MESH_INPUT_DIR} ${MESH_OUTPUT_DIR} ${UTILITIES_DIR} ${PDNEIGH_DIR} ${TRILINOS_DIR}/include )

#PdITIOperators library
add_library(PdITIOperators 
NewmarkBetaIntegrator.cxx
#ImplicitDynamicsIntegrator.cxx
ImplicitLinearDynamicsIntegrator.cxx
#PdITI_Utilities.cxx
#BodyLoadSpec.cxx
#ComponentDirichletBcSpec.cxx
#StageComponentDirichletBc.cxx
#IsotropicElasticConstitutiveModel.cxx
#IsotropicElastic_No_DSF.cxx
#IsotropicElasticPlasticModel.cxx
)

#
#add_executable(epetraDemo epetraDemo.cxx)
#target_link_libraries(epetraDemo 
#PdITIOperators
#PdNeigh
#PdVTK 
#${VTK_LIBS}
#${TRILINOS_LIBS}
#${REQUIRED_LIBS}
#)


