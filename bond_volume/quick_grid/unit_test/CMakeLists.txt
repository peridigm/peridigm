#
# QuickGrid Unit Tests
#                                                            
include_directories(${CMAKE_SOURCE_DIR}/io ${CMAKE_SOURCE_DIR}/bond_volume ${TRILINOS_DIR}/include)

add_executable(ut_ringCellVolumeFraction ut_ringCellVolumeFraction.cxx)
target_link_libraries(ut_ringCellVolumeFraction  QuickGrid  PdNeigh PdVTK Utilities QuickGridVolFrac ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (ut_ringCellVolumeFraction python run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

add_executable(ut_VolumeFraction ut_VolumeFraction.cxx)
target_link_libraries(ut_VolumeFraction  QuickGrid  PdNeigh PdVTK Utilities QuickGridVolFrac ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#add_test (ut_ringCellVolumeFraction python run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)


