
add_executable(ut_ringCellVolumeFraction ut_ringCellVolumeFraction.cxx)
target_link_libraries(ut_ringCellVolumeFraction  QuickGrid  PdNeigh PdVTK Utilities QuickGridVolFrac ${VTK_LIBS} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
#add_test (ut_ringCellVolumeFraction python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

add_executable(ut_VolumeFraction ut_VolumeFraction.cxx)
target_link_libraries(ut_VolumeFraction  QuickGrid  PdNeigh PdVTK Utilities QuickGridVolFrac ${VTK_LIBS} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
#add_test (ut_ringCellVolumeFraction python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

# Add test
add_executable(ut_bondVolumeConvergenceStudy ut_bondVolumeConvergenceStudy.cxx)
target_link_libraries(ut_bondVolumeConvergenceStudy QuickGridVolFrac QuickGrid PdNeigh Utilities PdVTK PdMaterialUtilities ${VTK_LIBS} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
add_test (ut_bondVolumeConvergenceStudy python ut_bondVolumeConvergenceStudy.py ./ut_bondVolumeConvergenceStudy)

# Add test
add_executable(ut_naiveQuadratureConvergenceStudy ut_naiveQuadratureConvergenceStudy.cxx)
target_link_libraries(ut_naiveQuadratureConvergenceStudy QuickGridVolFrac
QuickGrid PdNeigh Utilities PdVTK PdMaterialUtilities ${VTK_LIBS} ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
add_test (ut_naiveQuadratureConvergenceStudy python ut_naiveQuadratureConvergenceStudy.py ./ut_naiveQuadratureConvergenceStudy)

