
# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS} -d discretization/unit_test)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${IO_DIR} -b ${CMAKE_BINARY_DIR}/src/io ${UT_DIRS} )
endif (NOT ${buildInSource})

add_executable(utPeridigm_PdQuickGridDiscretization
               ${DISCRETIZATION_DIR}/Peridigm_Discretization.cpp
               ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
               ./utPeridigm_PdQuickGridDiscretization.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization
  ${Peridigm_LIBRARY}
  ${PDNEIGH_LIBS}
  ${Trilinos_LIBRARIES}
  ${MESH_INPUT_LIBS}
  ${REQUIRED_LIBS}
)
add_test (utPeridigm_PdQuickGridDiscretization python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_PdQuickGridDiscretization)

add_executable(utPeridigm_PdQuickGridDiscretization_MPI_np2
               ${DISCRETIZATION_DIR}/Peridigm_Discretization.cpp
               ${DISCRETIZATION_DIR}/Peridigm_PdQuickGridDiscretization.cpp
               ./utPeridigm_PdQuickGridDiscretization_MPI_np2.cpp)
target_link_libraries(utPeridigm_PdQuickGridDiscretization_MPI_np2
  ${Peridigm_LIBRARY}
  ${PDNEIGH_LIBS}
  ${Trilinos_LIBRARIES}
  ${MESH_INPUT_LIBS}
  ${REQUIRED_LIBS}
)
add_test (utPeridigm_PdQuickGridDiscretization_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_PdQuickGridDiscretization_MPI_np2)

add_executable(utPeridigm_ExodusDiscretization
               ${DISCRETIZATION_DIR}/Peridigm_Discretization.cpp
               ${DISCRETIZATION_DIR}/Peridigm_ExodusDiscretization.cpp
               ${IO_DIR}/Peridigm_ProximitySearch.cpp
               ./utPeridigm_ExodusDiscretization.cpp)
target_link_libraries(utPeridigm_ExodusDiscretization
  ${Peridigm_LIBRARY}
  ${PDNEIGH_LIBS}
  ${MESH_INPUT_LIBS}
  ${UTILITIES_LIBS}
  ${Trilinos_LIBRARIES}
  ${REQUIRED_LIBS}
)
add_test (utPeridigm_ExodusDiscretization python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ExodusDiscretization)
add_test (utPeridigm_ExodusDiscretization_MPI_np2 python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 2 ./utPeridigm_ExodusDiscretization)

add_executable(utPeridigm_GeometryUtils
               ${DISCRETIZATION_DIR}/Peridigm_GeometryUtils.cpp
               ./utPeridigm_GeometryUtils.cpp)
target_link_libraries(utPeridigm_GeometryUtils
  ${Peridigm_LIBRARY}
  ${Trilinos_LIBRARIES}
)
add_test (utPeridigm_GeometryUtils python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_GeometryUtils)

IF(PERIDIGM_PV)
  add_executable(utPeridigm_PartialVolumeCalculator
                 ${DISCRETIZATION_DIR}/Peridigm_PartialVolumeCalculator.cpp
                 ./utPeridigm_PartialVolumeCalculator.cpp)
  target_link_libraries(utPeridigm_PartialVolumeCalculator
    ${Peridigm_LIBRARY}
    ${Trilinos_LIBRARIES}
    ${MESH_INPUT_LIBS}
  )
  add_test (utPeridigm_PartialVolumeCalculator python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_PartialVolumeCalculator)

  add_executable(utPeridigm_PartialVolumeConvergence
                 ${DISCRETIZATION_DIR}/Peridigm_PartialVolumeCalculator.cpp
                 ./utPeridigm_PartialVolumeConvergence.cpp)
  target_link_libraries(utPeridigm_PartialVolumeConvergence
    ${Peridigm_LIBRARY}
    ${Trilinos_LIBRARIES}
    ${MESH_INPUT_LIBS}
  )
  add_test (utPeridigm_PartialVolumeConvergence python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_PartialVolumeConvergence)
ENDIF()
