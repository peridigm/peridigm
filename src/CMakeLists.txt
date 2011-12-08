#
# Main executable for Peridigm
#
include_directories(${TRILINOS_DIR}/include)

# Glob the compute source files
file (GLOB computeSOURCES 
  "${CMAKE_SOURCE_DIR}/src/compute/*.cpp"
)
# Tell CMake a header file won't be available until build time.
# set_source_files_properties(./compute/compute_includes.hpp PROPERTIES GENERATED 1)

set(Peridigm_SOURCES
    ./core/Peridigm_Main.cpp
    ./core/Peridigm_Factory.cpp
    ./core/Peridigm.cpp
    ./io/discretization/Peridigm_AbstractDiscretization.cpp
    ./io/discretization/Peridigm_DiscretizationFactory.cpp
    ./io/discretization/Peridigm_PdQuickGridDiscretization.cpp
    ./core/Peridigm_State.cpp
    ./core/Peridigm_DataManager.cpp
    ./core/Peridigm_SerialMatrix.cpp
    ./core/Peridigm_Block.cpp
    ./core/Peridigm_ModelEvaluator.cpp
    ./core/PHAL_Dimension.cpp
    ./core/PHAL_PeridigmTraits.cpp
    ./core/Peridigm_ComputeManager.cpp
    ./core/Peridigm_Timer.cpp
    ${computeSOURCES}
    ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp
    ./io/Peridigm_OutputManager_VTK_XML.cpp
    ./materials/Peridigm_Material.cpp
    ./materials/Peridigm_MaterialFactory.cpp
    ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
    ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
    ./materials/Peridigm_ViscoelasticStandardLinearSolid.cpp
    ./materials/Peridigm_DamageModelFactory.cpp
    ./materials/Peridigm_CriticalStretchDamageModel.cpp
    ./contact/Peridigm_ContactModelFactory.cpp
    ./contact/Peridigm_ShortRangeForceContactModel.cpp
)
if(PERIDIGM_STK)
  set(Peridigm_SOURCES ./io/discretization/Peridigm_STKDiscretization.cpp ${Peridigm_SOURCES})
endif()

set(Peridigm_LINK_LIBRARIES
    ${Trilinos_LIBRARIES}
    ${PDNEIGH_LIBS}
    ${Zoltan_LIBRARY}
    ${PARSER_LIBS}
    ${MESH_INPUT_LIBS}
    ${UTILITIES_LIBS}
    ${VTK_LIBS}
    ${MESH_OUTPUT_LIBS}
    ${REQUIRED_LIBS}
    ${PdMaterialUtilitiesLib}
    PdIC
)

add_executable(Peridigm ${Peridigm_SOURCES})
target_link_libraries(Peridigm ${Peridigm_LINK_LIBRARIES})

#
# Install target for Peridigm main executable
#
install(TARGETS Peridigm
    RUNTIME DESTINATION bin
)
set_property(TARGET Peridigm
    PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE
)

#
# Unit Tests
#
set(utPeridigm_SOURCES
    ./core/Peridigm.cpp
    ./io/discretization/Peridigm_AbstractDiscretization.cpp
    ./io/discretization/Peridigm_DiscretizationFactory.cpp
    ./io/discretization/Peridigm_PdQuickGridDiscretization.cpp
    ./core/Peridigm_State.cpp
    ./core/Peridigm_DataManager.cpp
    ./core/Peridigm_SerialMatrix.cpp
    ./core/Peridigm_Block.cpp
    ./core/Peridigm_ModelEvaluator.cpp
    ./core/PHAL_Dimension.cpp
    ./core/PHAL_PeridigmTraits.cpp
    ./core/Peridigm_ComputeManager.cpp
    ./core/Peridigm_Timer.cpp
    ${computeSOURCES}
    ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp
    ./io/Peridigm_OutputManager_VTK_XML.cpp
    ./io/mesh_output/Field.cxx
    ./materials/Peridigm_Material.cpp
    ./materials/Peridigm_MaterialFactory.cpp
    ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
    ./materials/Peridigm_IsotropicElasticPlasticMaterial.cpp
    ./materials/Peridigm_ViscoelasticStandardLinearSolid.cpp
    ./materials/Peridigm_DamageModelFactory.cpp
    ./materials/Peridigm_CriticalStretchDamageModel.cpp
    ./contact/Peridigm_ContactModelFactory.cpp
    ./contact/Peridigm_ShortRangeForceContactModel.cpp
)
if(PERIDIGM_STK)
  set(utPeridigm_SOURCES ./io/discretization/Peridigm_STKDiscretization.cpp ${utPeridigm_SOURCES})
endif()

add_executable(utPeridigm ./core/unit_test/utPeridigm.cpp ${utPeridigm_SOURCES})
target_link_libraries(utPeridigm ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(utPeridigm_MPI_np2 ./core/unit_test/utPeridigm_MPI_np2.cpp ${utPeridigm_SOURCES})
target_link_libraries(utPeridigm_MPI_np2 ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(utPeridigm_PdQuickGridDiscretization
               ./io/discretization/Peridigm_AbstractDiscretization.cpp
               ./io/discretization/Peridigm_PdQuickGridDiscretization.cpp
               ./core/unit_test/utPeridigm_PdQuickGridDiscretization.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization
                      ${Trilinos_LIBRARIES}
                      ${PDNEIGH_LIBS}
                      ${Zoltan_LIBRARY}
                      ${MESH_INPUT_LIBS}
                      ${MESH_OUTPUT_LIBS}
                      ${VTK_LIBS}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

add_executable(utPeridigm_PdQuickGridDiscretization_MPI_np2
               ./io/discretization/Peridigm_AbstractDiscretization.cpp
               ./io/discretization/Peridigm_PdQuickGridDiscretization.cpp
               ./core/unit_test/utPeridigm_PdQuickGridDiscretization_MPI_np2.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization_MPI_np2
                      ${Trilinos_LIBRARIES}
                      ${PDNEIGH_LIBS}
                      ${Zoltan_LIBRARY}
                      ${MESH_INPUT_LIBS}
                      ${MESH_OUTPUT_LIBS}
                      ${VTK_LIBS}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

add_executable(utPeridigm_State
               ./core/Peridigm_State.cpp
               ./core/unit_test/utPeridigm_State.cpp
               ./io/mesh_output/Field.cxx)
target_link_libraries(utPeridigm_State
                      ${Trilinos_LIBRARIES}
                      ${MESH_OUTPUT_LIBS}
                      ${REQUIRED_LIBS}
                      ${Boost_LIBRARIES}
)

add_executable(utPeridigm_LinearElasticIsotropicMaterial
               ./core/Peridigm_Timer.cpp
               ./core/Peridigm_State.cpp
               ./core/Peridigm_DataManager.cpp
               ./core/Peridigm_SerialMatrix.cpp
               ./io/mesh_output/Field.cxx
               ./materials/Peridigm_Material.cpp
	       ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
               ./materials/Peridigm_DamageModelFactory.cpp
               ./materials/Peridigm_CriticalStretchDamageModel.cpp
	       ./materials/utPeridigm_LinearElasticIsotropicMaterial.cpp
)
target_link_libraries(utPeridigm_LinearElasticIsotropicMaterial
${Trilinos_LIBRARIES}
${MESH_OUTPUT_LIBS}
${PdMaterialUtilitiesLib}
${REQUIRED_LIBS}
${Boost_LIBRARIES}
)

add_executable(utPHAL_UpdateForceState
               ./core/PHAL_PeridigmTraits.cpp
	       ./core/PHAL_Dimension.cpp
               ./core/Peridigm_Timer.cpp
               ./core/Peridigm_Block.cpp
               ./core/Peridigm_State.cpp
               ./core/Peridigm_DataManager.cpp
               ./core/Peridigm_SerialMatrix.cpp
               ./io/mesh_output/Field.cxx
               ./materials/Peridigm_Material.cpp
               ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
               ./materials/Peridigm_DamageModelFactory.cpp
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
	       ./core/PHAL_PeridigmTraits.cpp
	       ./core/PHAL_Dimension.cpp
               ./core/Peridigm_Timer.cpp
               ./core/Peridigm_State.cpp
               ./core/Peridigm_DataManager.cpp
               ./core/Peridigm_Block.cpp
               ./core/Peridigm_SerialMatrix.cpp
               ./io/mesh_output/Field.cxx
               ./materials/Peridigm_Material.cpp
               ./materials/Peridigm_LinearElasticIsotropicMaterial.cpp
               ./materials/Peridigm_DamageModelFactory.cpp
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

set(utPeridigm_Compute_Force_SOURCES
   ${CMAKE_SOURCE_DIR}/src/compute/unitTest/utPeridigm_Compute_Force.cpp
   ./io/discretization/Peridigm_AbstractDiscretization.cpp
   ./io/discretization/Peridigm_DiscretizationFactory.cpp
   ./io/discretization/Peridigm_PdQuickGridDiscretization.cpp
   ./core/Peridigm_DataManager.cpp
   ./core/Peridigm_State.cpp
   ./core/Peridigm_ComputeManager.cpp
   ${computeSOURCES}
   ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp
)
if(PERIDIGM_STK)
   set(utPeridigm_Compute_Force_SOURCES ./io/discretization/Peridigm_STKDiscretization.cpp ${utPeridigm_Compute_Force_SOURCES})
endif()
add_executable(utPeridigm_Compute_Force 
   ${utPeridigm_Compute_Force_SOURCES}
)
target_link_libraries(utPeridigm_Compute_Force 
   ${Peridigm_LINK_LIBRARIES} 
   ${Boost_LIBRARIES}
)

# 
# Add subdirectories
#
add_subdirectory (core/)
add_subdirectory (io/)
add_subdirectory (materials/)
add_subdirectory (initial_conditions/)
add_subdirectory (compute/)
add_subdirectory (muParser/)

#
# Custom Commands
# 
#file (GLOB Peridigm_COMPUTE_HEADERS "${CMAKE_SOURCE_DIR}/src/compute/Peridigm*.hpp")
#list(REMOVE_ITEM Peridigm_COMPUTE_HEADERS ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp)
#add_custom_command(
#   OUTPUT  ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.hpp
#   DEPENDS ${Peridigm_COMPUTE_HEADERS}
#   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/compute
#   COMMAND sh ${CMAKE_SOURCE_DIR}/src/compute/compute_includes.sh
#)

