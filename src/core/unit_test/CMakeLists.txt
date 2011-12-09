#
# Unit Tests
#

set(utPeridigm_SOURCES
    ${CORE_DIR}/Peridigm.cpp
    ${CORE_DIR}/Peridigm_State.cpp
    ${CORE_DIR}/Peridigm_DataManager.cpp
    ${CORE_DIR}/Peridigm_SerialMatrix.cpp
    ${CORE_DIR}/Peridigm_Block.cpp
    ${CORE_DIR}/Peridigm_ModelEvaluator.cpp
    ${CORE_DIR}/PHAL_Dimension.cpp
    ${CORE_DIR}/PHAL_PeridigmTraits.cpp
    ${CORE_DIR}/Peridigm_ComputeManager.cpp
    ${CORE_DIR}/Peridigm_Timer.cpp
    ${CORE_DIR}/Peridigm_RandomNumber.cpp
    ${COMPUTE_DIR}/compute_includes.hpp
    ${CONTACT_DIR}/Peridigm_ContactModelFactory.cpp
    ${CONTACT_DIR}/Peridigm_ShortRangeForceContactModel.cpp
    ${DISCRETIZATION_DIR}/Peridigm_AbstractDiscretization.cpp
    ${DISCRETIZATION_DIR}/Peridigm_DiscretizationFactory.cpp
    ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
    ${IO_DIR}/Peridigm_OutputManager_VTK_XML.cpp
    ${MATERIALS_DIR}/Peridigm_Material.cpp
    ${MATERIALS_DIR}/Peridigm_CriticalStretchDamageModel.cpp
    ${MATERIALS_DIR}/Peridigm_DamageModelFactory.cpp
    ${MATERIALS_DIR}/Peridigm_MaterialFactory.cpp
    ${MATERIALS_DIR}/Peridigm_LinearElasticIsotropicMaterial.cpp
    ${MATERIALS_DIR}/Peridigm_IsotropicElasticPlasticMaterial.cpp
    ${MATERIALS_DIR}/Peridigm_ViscoelasticStandardLinearSolid.cpp
    ${MESH_OUTPUT_DIR}/Field.cxx
    ${computeSOURCES}
)
if(PERIDIGM_STK)
  set(utPeridigm_SOURCES ${DISCRETIZATION_DIR}/Peridigm_STKDiscretization.cpp ${utPeridigm_SOURCES})
endif()

add_executable(utPeridigm utPeridigm.cpp ${utPeridigm_SOURCES})
target_link_libraries(utPeridigm ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})
add_test (utPeridigm python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm)

add_executable(utPeridigm_MPI_np2 utPeridigm_MPI_np2.cpp ${utPeridigm_SOURCES})
target_link_libraries(utPeridigm_MPI_np2 ${Peridigm_LINK_LIBRARIES} ${Boost_LIBRARIES})
add_test (utPeridigm_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_MPI_np2)

add_executable(utPeridigm_State ${CORE_DIR}/Peridigm_State.cpp ${MESH_OUTPUT_DIR}/Field.cxx ./utPeridigm_State.cpp)
target_link_libraries(utPeridigm_State ${Trilinos_LIBRARIES} ${MESH_OUTPUT_LIBS} ${REQUIRED_LIBS} ${Boost_LIBRARIES})
add_test (utPeridigm_State python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_State)
add_test (utPeridigm_State_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_State)

