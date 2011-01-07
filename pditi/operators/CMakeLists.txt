#
# Operators for Pimp -- Peridynamics Implicit
#
#
# Collection of libs that we need to link with
#include_directories(${TRILINOS_DIR}/include  ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${PERIDIGM_IO_SRC} ${PERIDIGM_MAT_UTILS_SRC})
#link_directories(${TRILINOS_DIR}/lib ${VTK_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ${PERIDIGM_IO_LIB} ${PERIDIGM_MAT_UTILS_LIB})
#
#PimpOperators library
#add_library(PdImpOperators 
#PdImpOperator.cxx
#NewmarkBetaIntegrator.cxx
#ImplicitDynamicsIntegrator.cxx
#ImplicitLinearDynamicsIntegrator.cxx
#PdImpOperatorUtilities.cxx
#BodyLoadSpec.cxx
#ComponentDirichletBcSpec.cxx
#StageComponentDirichletBc.cxx
#IsotropicElasticConstitutiveModel.cxx
#IsotropicElasticPlasticModel.cxx
#)
#
#
#add_executable(epetraDemo epetraDemo.cxx)
#target_link_libraries(epetraDemo 
#PdImpOperators
#PdNeighborhood
#${PdMaterialUtilities_LIBRARY}
#${PdVTK_LIBRARY}
#${PdZoltan_LIBRARY}
#${PdQuickGrid_LIBRARY}
#${VTK_LIBS}
#${TRILINOS_LIBS}
#${REQUIRED_LIBS}
#${Boost_LIBRARIES}
#)
#

#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

#PimpOperators library
add_library(PdImpOperators 
PdImpOperator.cxx
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
PdImpOperators
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
