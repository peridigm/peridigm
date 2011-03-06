#
# Main executable for Peridigm
#

add_executable (${EXECUTABLE_NAME}
        Peridigm_Main.cpp
        Peridigm_Factory.cpp
        Peridigm.cpp
        Peridigm_DiscretizationFactory.cpp
        Peridigm_PdQuickGridDiscretization.cpp
        Peridigm_DataManager.cpp
        Peridigm_DenseMatrix.cpp
        Peridigm_ModelEvaluator.cpp
        PHAL_Dimension.cpp
        PHAL_PeridigmTraits.cpp
        Peridigm_OutputManager_VTK_XML.cpp
#        Peridigm_VerletSolver.cpp
#        Peridigm_VerletIntegrator.cpp
#        Peridigm_VerletObserver.cpp
        Peridigm_Timer.cpp
        ./materials/Peridigm_Material.cpp
        ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
        ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
        ./contact/Peridigm_ShortRangeForceContactModel.cpp
)

target_link_libraries (${EXECUTABLE_NAME}
${TRILINOS_LIBS}
${PdZoltanLibs}
${Zoltan_LIBRARY}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${PdMaterialUtilitiesLib}
PdIC
${REQUIRED_LIBS}
)

#
# Unit Tests
#

add_executable(utPeridigm
		utPeridigm.cpp
        Peridigm.cpp
        Peridigm_DiscretizationFactory.cpp
        Peridigm_PdQuickGridDiscretization.cpp
        Peridigm_DataManager.cpp
        Peridigm_DenseMatrix.cpp
        Peridigm_ModelEvaluator.cpp
        PHAL_Dimension.cpp
        PHAL_PeridigmTraits.cpp
        Peridigm_OutputManager_VTK_XML.cpp
        Peridigm_Timer.cpp
        ../io/Field.cxx
        ./materials/Peridigm_Material.cpp
        ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
        ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
        ./contact/Peridigm_ShortRangeForceContactModel.cpp
)
target_link_libraries(utPeridigm
${TRILINOS_LIBS}
${PdZoltanLibs}
${Zoltan_LIBRARY}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${REQUIRED_LIBS}
${PdMaterialUtilitiesLib}
PdIC
${Boost_LIBRARIES}
)

add_executable(utPeridigm_MPI_np2
		utPeridigm_MPI_np2.cpp
        Peridigm.cpp
        Peridigm_DiscretizationFactory.cpp
        Peridigm_PdQuickGridDiscretization.cpp
        Peridigm_DataManager.cpp
        Peridigm_DenseMatrix.cpp
        Peridigm_ModelEvaluator.cpp
        PHAL_Dimension.cpp
        PHAL_PeridigmTraits.cpp
        Peridigm_OutputManager_VTK_XML.cpp
        Peridigm_Timer.cpp
        ../io/Field.cxx
        ./materials/Peridigm_Material.cpp
        ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
        ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
        ./contact/Peridigm_ShortRangeForceContactModel.cpp
)
target_link_libraries(utPeridigm_MPI_np2
${TRILINOS_LIBS}
${PdZoltanLibs}
${Zoltan_LIBRARY}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${REQUIRED_LIBS}
${PdMaterialUtilitiesLib}
PdIC
${Boost_LIBRARIES}
)

add_executable(utPeridigm_PdQuickGridDiscretization
		utPeridigm_PdQuickGridDiscretization.cpp
        Peridigm_PdQuickGridDiscretization.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization
${TRILINOS_LIBS}
${PdZoltanLibs}
${Zoltan_LIBRARY}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)


add_executable(utPeridigm_PdQuickGridDiscretization_MPI_np2
		utPeridigm_PdQuickGridDiscretization_MPI_np2.cpp
        Peridigm_PdQuickGridDiscretization.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization_MPI_np2
${TRILINOS_LIBS}
${PdZoltanLibs}
${Zoltan_LIBRARY}
${PdQuickGridLibs}
${PdVTKLibs}
${VTK_LIBS}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPeridigm_LinearElasticIsotropicMaterial
        Peridigm_DataManager.cpp
        Peridigm_DenseMatrix.cpp
        ../io/Field.cxx
        ./materials/Peridigm_Material.cpp
		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
        ./materials/Peridigm_CriticalStretchDamageModel.cpp
		./materials/utPeridigm_LinearElasticIsotropicMaterial.cpp)
target_link_libraries(utPeridigm_LinearElasticIsotropicMaterial
${TRILINOS_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPHAL_UpdateForceState
		PHAL_PeridigmTraits.cpp
		PHAL_Dimension.cpp
        Peridigm_DataManager.cpp
        Peridigm_DenseMatrix.cpp
        ../io/Field.cxx
        ./materials/Peridigm_Material.cpp
		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
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
        Peridigm_DataManager.cpp
        Peridigm_DenseMatrix.cpp
        ../io/Field.cxx
        ./materials/Peridigm_Material.cpp
		./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
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
#       ./materials/Peridigm_CriticalStretchDamageModel.cpp
#		./evaluators/utPHAL_EvaluateForce_MPI_np2.cpp)
#target_link_libraries(utPHAL_EvaluateForce_MPI_np2
#${TRILINOS_LIBS}
#${PdMaterialUtilitiesLib}
#${REQUIRED_LIBS}
#${Boost_LIBRARIES}
#)

add_subdirectory (materials/)
add_subdirectory (initial_conditions/)
