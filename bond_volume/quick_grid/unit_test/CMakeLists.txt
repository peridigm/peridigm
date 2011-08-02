#
# QuickGrid Unit Tests
#                                                            
include_directories(${CMAKE_SOURCE_DIR}/material_evaluation ${CMAKE_SOURCE_DIR}/io ${CMAKE_SOURCE_DIR} ${TRILINOS_DIR}/include)

add_executable(ut_ringCellVolumeFraction ut_ringCellVolumeFraction.cxx)
target_link_libraries(ut_ringCellVolumeFraction  QuickGrid  PdNeigh PdVTK Utilities QuickGridVolFrac ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (ut_ringCellVolumeFraction python run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

add_executable(ut_VolumeFraction ut_VolumeFraction.cxx)
target_link_libraries(ut_VolumeFraction  QuickGrid  PdNeigh PdVTK Utilities QuickGridVolFrac ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (ut_ringCellVolumeFraction python run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

# Add test
add_executable(ut_bondVolumeConvergenceStudy ut_bondVolumeConvergenceStudy.cxx)
target_link_libraries(ut_bondVolumeConvergenceStudy QuickGridVolFrac QuickGrid PdNeigh Utilities PdVTK PdMaterialUtilities ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
add_test (ut_bondVolumeConvergenceStudy python run_unit_test.py ./ut_bondVolumeConvergenceStudy)

