#
# Operators for PdITI -- Peridynamics Implicit
#


include_directories(${UTILITIES_DIR} ${PDNEIGH_DIR} ${TRILINOS_DIR}/include )


#PimpOperators library
add_library(PdITIOperators 
PdImpOperator.cxx
PdITI_Operator.cxx
NewmarkBetaIntegrator.cxx
ImplicitDynamicsIntegrator.cxx
ImplicitLinearDynamicsIntegrator.cxx
PdITI_Utilities.cxx
BodyLoadSpec.cxx
ComponentDirichletBcSpec.cxx
StageComponentDirichletBc.cxx
IsotropicElasticConstitutiveModel.cxx
IsotropicElastic_No_DSF.cxx
IsotropicElasticPlasticModel.cxx
)

#
add_executable(epetraDemo epetraDemo.cxx)
target_link_libraries(epetraDemo 
PdITIOperators
PdNeigh
${PERIDIGM_LIBS}
${VTK_LIBS}
${TRILINOS_LIBS}
${REQUIRED_LIBS}
)

#set (REQUIRED_LIBS
#liblapack.a
#libblas.a
#libmpich.a
#${PTHREAD_LIB}
#libgfortran.a
#libm.a
#)
#
