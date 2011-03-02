#
# Operators for PdITI -- Peridynamics Implicit
#


#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

#PimpOperators library
add_library(PdITIOperators 
PdImpOperator.cxx
PdITI_Operator.cxx
NewmarkBetaIntegrator.cxx
ImplicitDynamicsIntegrator.cxx
ImplicitLinearDynamicsIntegrator.cxx
PdImpOperatorUtilities.cxx
BodyLoadSpec.cxx
ComponentDirichletBcSpec.cxx
StageComponentDirichletBc.cxx
IsotropicElasticConstitutiveModel.cxx
IsotropicElasticPlasticModel.cxx
)

#
add_executable(epetraDemo epetraDemo.cxx)
target_link_libraries(epetraDemo 
PdITIOperators
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
