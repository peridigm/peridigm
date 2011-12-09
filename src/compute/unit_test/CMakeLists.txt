set(utPeridigm_Compute_Force_SOURCES
   ${DISCRETIZATION_DIR}/Peridigm_AbstractDiscretization.cpp
   ${DISCRETIZATION_DIR}/Peridigm_DiscretizationFactory.cpp
   ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
   ${CORE_DIR}/Peridigm_DataManager.cpp
   ${CORE_DIR}/Peridigm_State.cpp
   ${CORE_DIR}/Peridigm_ComputeManager.cpp
   ${COMPUTE_DIR}/compute_includes.hpp
   ${computeSOURCES}
   ./utPeridigm_Compute_Force.cpp
)
if(PERIDIGM_STK)
   set(utPeridigm_Compute_Force_SOURCES ${DISCRETIZATION_DIR}/Peridigm_STKDiscretization.cpp ${utPeridigm_Compute_Force_SOURCES})
endif()

add_executable(utPeridigm_Compute_Force ${utPeridigm_Compute_Force_SOURCES})
target_link_libraries(utPeridigm_Compute_Force
   ${Peridigm_LINK_LIBRARIES}
   ${Boost_LIBRARIES}
)


add_test (utPeridigm_Compute_Force python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Compute_Force)
add_test (utPeridigm_Compute_Force_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_Compute_Force)

