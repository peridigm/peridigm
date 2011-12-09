add_executable(utPeridigm_PdQuickGridDiscretization
               ${DISCRETIZATION_DIR}/Peridigm_AbstractDiscretization.cpp
               ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
               ./utPeridigm_PdQuickGridDiscretization.cpp)
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
add_test (utPeridigm_PdQuickGridDiscretization python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_PdQuickGridDiscretization)

add_executable(utPeridigm_PdQuickGridDiscretization_MPI_np2
               ${DISCRETIZATION_DIR}/Peridigm_AbstractDiscretization.cpp
               ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
               ./utPeridigm_PdQuickGridDiscretization_MPI_np2.cpp)
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
add_test (utPeridigm_PdQuickGridDiscretization_MPI_np2  python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_PdQuickGridDiscretization_MPI_np2)

