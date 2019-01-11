# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d bond_volume/quick_grid/unit_test)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${IO_DIR} -b ${CMAKE_BINARY_DIR}/src/io ${UT_DIRS} )
endif (NOT ${buildInSource})

if(USE_YAML)
  add_executable(ut_ringCellVolumeFraction ut_ringCellVolumeFraction.cxx)
  target_link_libraries(ut_ringCellVolumeFraction QuickGrid PdNeigh PdField Utilities QuickGridVolFrac ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
  #add_test (ut_ringCellVolumeFraction python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

  add_executable(ut_VolumeFraction ut_VolumeFraction.cxx)
  target_link_libraries(ut_VolumeFraction QuickGrid PdNeigh PdField Utilities QuickGridVolFrac ${Trilinos_LIBRARIES} ${UT_REQUIRED_LIBS})
  #add_test (ut_ringCellVolumeFraction python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py mpiexec -np 1 ./ut_ringCellVolumeFraction)

  # Add test
  add_executable(ut_bondVolumeConvergenceStudy ut_bondVolumeConvergenceStudy.cxx)
  target_link_libraries(ut_bondVolumeConvergenceStudy QuickGridVolFrac QuickGrid PdNeigh PdField Utilities PdMaterialUtilities ${Trilinos_LIBRARIES} ${PARSER_LIBS} ${UT_REQUIRED_LIBS})
  #add_test (ut_bondVolumeConvergenceStudy python ut_bondVolumeConvergenceStudy.py ./ut_bondVolumeConvergenceStudy)

  # Add test
  add_executable(ut_naiveQuadratureConvergenceStudy ut_naiveQuadratureConvergenceStudy.cxx)
  target_link_libraries(ut_naiveQuadratureConvergenceStudy QuickGridVolFrac
  QuickGrid PdNeigh PdField Utilities PdMaterialUtilities ${Trilinos_LIBRARIES} ${PARSER_LIBS} ${UT_REQUIRED_LIBS})
  #add_test (ut_naiveQuadratureConvergenceStudy python ut_naiveQuadratureConvergenceStudy.py ./ut_naiveQuadratureConvergenceStudy)
endif()
