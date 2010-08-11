#
# Main executable for Peridigm
#

include_directories (../io ${TRILINOS_DIR}/include ${Boost_INCLUDE_DIRS} ./evaluators ./materials ./contact ./ ../PdMaterialUtilities)
link_directories (${TRILINOS_DIR}/lib ../io ../PdMaterialUtilities ${VTK_LIBRARY_DIRS})

##add_executable (${EXECUTABLE_NAME}
##ENAT_RythmosSolver.cpp
##PHAL_Dimension.cpp
##PHAL_PeridigmTraits.cpp
##Peridigm.cpp
##Peridigm_DiscretizationFactory.cpp
##Peridigm_ModelEvaluator.cpp
##Peridigm_SolverFactory.cpp
##Peridigm_PdQuickGridDiscretization.cpp
##Peridigm_RythmosObserver.cpp
##Peridigm_OutputManager_VTK_XML.cpp
##Peridigm_VerletSolver.cpp
##Peridigm_VerletIntegrator.cpp
##Peridigm_VerletObserver.cpp
##./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
##./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
##./materials/Peridigm_DecompositionStates.cpp
##./materials/Peridigm_CriticalStretchDamageModel.cpp
##./contact/Peridigm_ShortRangeForceContactModel.cpp
##)

##target_link_libraries (${EXECUTABLE_NAME}
##${TRILINOS_LIBS}
##${PdZoltanLibs}
##${Zoltan_LIBRARY}
##${PdQuickGridLibs}
##${PdVTKLibs}
##${VTK_LIBS}
##${PdMaterialUtilitiesLib}
##${REQUIRED_LIBS}
##)

#
# Unit Tests
#

##add_executable(utPeridigm_PdQuickGridDiscretization
##		utPeridigm_PdQuickGridDiscretization.cpp)
##target_link_libraries(utPeridigm_PdQuickGridDiscretization
##${TRILINOS_LIBS}
##${PdQuickGridLibs}
##${REQUIRED_LIBS}
##${Boost_LIBRARIES}
##)

add_executable(utPeridigm_LinearElasticIsotropicMaterial
		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
		./materials/Peridigm_DecompositionStates.cpp
         ./materials/Peridigm_CriticalStretchDamageModel.cpp
		./materials/utPeridigm_LinearElasticIsotropicMaterial.cpp)
target_link_libraries(utPeridigm_LinearElasticIsotropicMaterial
${TRILINOS_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPeridigm_IsotropicElasticPlasticMaterial
		./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
		./materials/Peridigm_DecompositionStates.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
		./materials/utPeridigm_IsotropicElasticPlasticMaterial.cpp)
target_link_libraries(utPeridigm_IsotropicElasticPlasticMaterial
${TRILINOS_LIBS}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPeridigm_decompositionStates
		./materials/Peridigm_DecompositionStates.cpp
		./materials/utPeridigm_decompositionStates.cpp)
target_link_libraries(utPeridigm_decompositionStates
${TRILINOS_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPHAL_UpdateForceState
		PHAL_PeridigmTraits.cpp
		PHAL_Dimension.cpp
		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
		./materials/Peridigm_DecompositionStates.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
		./evaluators/utPHAL_UpdateForceState.cpp)
target_link_libraries(utPHAL_UpdateForceState
${TRILINOS_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

#add_executable(utPHAL_UpdateForceState_MPI_np2
#		PHAL_PeridigmTraits.cpp
#		PHAL_Dimension.cpp
#		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
#       ./materials/Peridigm_DecompositionStates.cpp
#       ./materials/Peridigm_CriticalStretchDamageModel.cpp
#		./evaluators/utPHAL_UpdateForceState_MPI_np2.cpp)
#target_link_libraries(utPHAL_UpdateForceState_MPI_np2
#${TRILINOS_LIBS}
#${PdMaterialUtilitiesLib}
#${REQUIRED_LIBS}
#${Boost_LIBRARIES}
#)

add_executable(utPHAL_EvaluateForce
		PHAL_PeridigmTraits.cpp
		PHAL_Dimension.cpp
		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
		./materials/Peridigm_DecompositionStates.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
		./evaluators/utPHAL_EvaluateForce.cpp)
target_link_libraries(utPHAL_EvaluateForce
${TRILINOS_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

#add_executable(utPHAL_EvaluateForce_MPI_np2
#		PHAL_PeridigmTraits.cpp
#        PHAL_Dimension.cpp
#		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
#       ./materials/Peridigm_DecompositionStates.cpp
#       ./materials/Peridigm_CriticalStretchDamageModel.cpp
#		./evaluators/utPHAL_EvaluateForce_MPI_np2.cpp)
#target_link_libraries(utPHAL_EvaluateForce_MPI_np2
#${TRILINOS_LIBS}
#${PdMaterialUtilitiesLib}
#${REQUIRED_LIBS}
#${Boost_LIBRARIES}
#)
