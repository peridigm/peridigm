add_subdirectory (twoPoint_SLS_Relaxation/)
add_subdirectory (utPeridigm_IsotropicElasticPlastic/)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unitTest/twoPoint_SLS_Relaxation)
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unitTest/utPeridigm_IsotropicElasticPlastic)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})
