add_subdirectory (utPeridigm_IsotropicElasticPlastic/)
add_subdirectory (twoPointMaxwellRelaxation/)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unitTest/utPeridigm_IsotropicElasticPlastic)
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unitTest/twoPointMaxwellRelaxation)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})
