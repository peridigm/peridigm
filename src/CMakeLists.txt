#
# Main executable for Peridigm
#

# Glob the compute source files
file (GLOB computeSOURCES 
  "${COMPUTE_DIR}/*.cpp"
)

set(Peridigm_SOURCES
    ${COMPUTE_DIR}/compute_includes.hpp
    ${CONTACT_DIR}/Peridigm_ContactModelFactory.cpp
    ${CONTACT_DIR}/Peridigm_ShortRangeForceContactModel.cpp
#    ${CORE_DIR}/Peridigm_Main.cpp
    ${CORE_DIR}/Peridigm_Factory.cpp
    ${CORE_DIR}/Peridigm.cpp
    ${CORE_DIR}/Peridigm_State.cpp
    ${CORE_DIR}/Peridigm_DataManager.cpp
    ${CORE_DIR}/Peridigm_SerialMatrix.cpp
    ${CORE_DIR}/Peridigm_Block.cpp
    ${CORE_DIR}/Peridigm_BoundaryAndInitialConditionManager.cpp
    ${CORE_DIR}/Peridigm_CriticalTimeStep.cpp
    ${CORE_DIR}/Peridigm_ModelEvaluator.cpp
    ${CORE_DIR}/PHAL_Dimension.cpp
    ${CORE_DIR}/PHAL_PeridigmTraits.cpp
    ${CORE_DIR}/Peridigm_ComputeManager.cpp
    ${CORE_DIR}/Peridigm_RandomNumber.cpp
    ${CORE_DIR}/Peridigm_Timer.cpp
    ${DISCRETIZATION_DIR}/Peridigm_AbstractDiscretization.cpp
    ${DISCRETIZATION_DIR}/Peridigm_DiscretizationFactory.cpp
    ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
    ${DISCRETIZATION_DIR}/Peridigm_PartialVolumeCalculator.cpp
    ${DISCRETIZATION_DIR}/Peridigm_STKDiscretization.cpp
    ${IO_DIR}/Peridigm_OutputManager_VTK_XML.cpp
    ${IO_DIR}/Peridigm_OutputManager_ExodusII.cpp
    ${MATERIALS_DIR}/Peridigm_Material.cpp
    ${MATERIALS_DIR}/Peridigm_MaterialFactory.cpp
    ${MATERIALS_DIR}/Peridigm_LinearElasticIsotropicMaterial.cpp
    ${MATERIALS_DIR}/Peridigm_IsotropicElasticPlasticMaterial.cpp
    ${MATERIALS_DIR}/Peridigm_IsotropicHardeningPlasticMaterial.cpp
    ${MATERIALS_DIR}/Peridigm_ViscoelasticStandardLinearSolid.cpp
    ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
    ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
    ${computeSOURCES}
)

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
)

add_executable(Peridigm ${CORE_DIR}/Peridigm_Main.cpp ${Peridigm_SOURCES})
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
# Add subdirectories
#
add_subdirectory (core/)
add_subdirectory (io/)
add_subdirectory (materials/)
add_subdirectory (compute/)
add_subdirectory (muParser/)
