#
# Main executable for Peridigm
#
include_directories(${TRILINOS_DIR}/include)

# Glob the compute source files
file (GLOB computeSOURCES "${CMAKE_SOURCE_DIR}/src/compute/*.cpp")
# Tell CMake a header file won't be available until build time.
set_source_files_properties(./compute/compute_includes.hpp PROPERTIES GENERATED 1)

set(Peridigm_SOURCES
    Peridigm_Main.cpp
    Peridigm_Factory.cpp
    Peridigm.cpp
    Peridigm_DiscretizationFactory.cpp
    Peridigm_PdQuickGridDiscretization.cpp
    Peridigm_State.cpp
    Peridigm_DataManager.cpp
    Peridigm_SerialMatrix.cpp
    Peridigm_ModelEvaluator.cpp
    PHAL_Dimension.cpp
    PHAL_PeridigmTraits.cpp
    Peridigm_OutputManager_VTK_XML.cpp
    Peridigm_ComputeManager.cpp
    Peridigm_Timer.cpp
    ${computeSOURCES}
    ./materials/Peridigm_Material.cpp
    ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
    ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
    ./materials/Peridigm_ViscoelasticStandardLinearSolid.cpp
    ./materials/Peridigm_CriticalStretchDamageModel.cpp
    ./contact/Peridigm_ShortRangeForceContactModel.cpp
)
if(PERIDIGM_STK)
  set(Peridigm_SOURCES Peridigm_STKDiscretization.cpp ${Peridigm_SOURCES})
endif()

set(Peridigm_LINK_LIBRARIES
    ${Trilinos_LIBRARIES}
    ${PdZoltanLibs}
    ${Zoltan_LIBRARY}
    ${PdQuickGridLibs}
    ${PdVTKLibs}
    ${VTK_LIBS}
    ${REQUIRED_LIBS}
    ${PdMaterialUtilitiesLib}
    PdIC
)

add_executable(Peridigm ${Peridigm_SOURCES})
target_link_libraries(Peridigm ${Peridigm_LINK_LIBRARIES})
# Make sure the header is generated before the executable builds.
add_dependencies(Peridigm generate_compute_includes)

#
# Unit Tests
#
set(utPeridigm_SOURCES
    Peridigm.cpp
    Peridigm_DiscretizationFactory.cpp
    Peridigm_PdQuickGridDiscretization.cpp
    Peridigm_State.cpp
    Peridigm_DataManager.cpp
    Peridigm_SerialMatrix.cpp
    Peridigm_ModelEvaluator.cpp
    PHAL_Dimension.cpp
    PHAL_PeridigmTraits.cpp
    Peridigm_OutputManager_VTK_XML.cpp
    Peridigm_ComputeManager.cpp
    Peridigm_Timer.cpp
    ${computeSOURCES}
    ../io/Field.cxx
    ./materials/Peridigm_Material.cpp
    ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
    ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
    ./materials/Peridigm_ViscoelasticStandardLinearSolid.cpp
    ./materials/Peridigm_CriticalStretchDamageModel.cpp
    ./contact/Peridigm_ShortRangeForceContactModel.cpp
)
if(PERIDIGM_STK)
  set(utPeridigm_SOURCES Peridigm_STKDiscretization.cpp ${utPeridigm_SOURCES})
endif()

add_executable(utPeridigm utPeridigm.cpp ${utPeridigm_SOURCES})
target_link_libraries(utPeridigm ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})
add_dependencies(utPeridigm generate_compute_includes)

add_executable(utPeridigm_MPI_np2 utPeridigm_MPI_np2.cpp ${utPeridigm_SOURCES})
target_link_libraries(utPeridigm_MPI_np2 ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})
add_dependencies(utPeridigm_MPI_np2 generate_compute_includes)

add_executable(utPeridigm_PdQuickGridDiscretization
	       utPeridigm_PdQuickGridDiscretization.cpp
               Peridigm_PdQuickGridDiscretization.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization
                      ${Trilinos_LIBRARIES}
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
                      ${Trilinos_LIBRARIES}
                      ${PdZoltanLibs}
                      ${Zoltan_LIBRARY}
                      ${PdQuickGridLibs}
                      ${PdVTKLibs}
                      ${VTK_LIBS}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

add_executable(utPeridigm_State
               Peridigm_State.cpp
               utPeridigm_State.cpp
               ../io/Field.cxx)
target_link_libraries(utPeridigm_State
                      ${Trilinos_LIBRARIES}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

add_executable(utPeridigm_LinearElasticIsotropicMaterial
               Peridigm_State.cpp
               Peridigm_DataManager.cpp
               Peridigm_SerialMatrix.cpp
               ../io/Field.cxx
               ./materials/Peridigm_Material.cpp
		       ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
               ./materials/Peridigm_CriticalStretchDamageModel.cpp
		       ./materials/utPeridigm_LinearElasticIsotropicMaterial.cpp
)
target_link_libraries(utPeridigm_LinearElasticIsotropicMaterial
${Trilinos_LIBRARIES}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPHAL_UpdateForceState
	           PHAL_PeridigmTraits.cpp
	           PHAL_Dimension.cpp
               Peridigm_State.cpp
               Peridigm_DataManager.cpp
               Peridigm_SerialMatrix.cpp
               ../io/Field.cxx
               ./materials/Peridigm_Material.cpp
	       ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
               ./materials/Peridigm_CriticalStretchDamageModel.cpp
	       ./evaluators/utPHAL_UpdateForceState.cpp)
target_link_libraries(utPHAL_UpdateForceState
                      ${Trilinos_LIBRARIES}
                      ${PdMaterialUtilitiesLib}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

#add_executable(utPHAL_UpdateForceState_MPI_np2
#	        PHAL_PeridigmTraits.cpp
#	        PHAL_Dimension.cpp
#	        ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
#               ./materials/Peridigm_CriticalStretchDamageModel.cpp
#	        ./evaluators/utPHAL_UpdateForceState_MPI_np2.cpp)
#target_link_libraries(utPHAL_UpdateForceState_MPI_np2
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)

add_executable(utPHAL_EvaluateForce
	           PHAL_PeridigmTraits.cpp
	           PHAL_Dimension.cpp
               Peridigm_State.cpp
               Peridigm_DataManager.cpp
               Peridigm_SerialMatrix.cpp
               ../io/Field.cxx
               ./materials/Peridigm_Material.cpp
	           ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
               ./materials/Peridigm_CriticalStretchDamageModel.cpp
	           ./evaluators/utPHAL_EvaluateForce.cpp
)
target_link_libraries(utPHAL_EvaluateForce
                      ${Trilinos_LIBRARIES}
                      ${PdMaterialUtilitiesLib}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

#add_executable(utPHAL_EvaluateForce_MPI_np2
#	        PHAL_PeridigmTraits.cpp
#               PHAL_Dimension.cpp
#	        ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
#               ./materials/Peridigm_CriticalStretchDamageModel.cpp
#	        ./evaluators/utPHAL_EvaluateForce_MPI_np2.cpp)
#target_link_libraries(utPHAL_EvaluateForce_MPI_np2
#                      ${Trilinos_LIBRARIES}
#                      ${PdMaterialUtilitiesLib}
#                      ${REQUIRED_LIBS}
#                      ${Boost_LIBRARIES}
#)

add_subdirectory (materials/)
add_subdirectory (initial_conditions/)
add_subdirectory (compute/)
