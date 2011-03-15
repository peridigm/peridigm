#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#
SET(PDITI_ROOT ${CMAKE_SOURCE_DIR}/pditi)

SET(UTILITIES_DIR ${PDITI_ROOT}/utilities)
SET(MESH_OUTPUT_DIR ${PDITI_ROOT}/mesh_output)
SET(PDNEIGH_DIR ${PDITI_ROOT}/pdneigh)
SET(MESH_INPUT_DIR ${PDITI_ROOT}/mesh_input)

# Project directories
add_subdirectory (aliti)
add_subdirectory (intrepid)
add_subdirectory (pdneigh)
add_subdirectory (utilities)
add_subdirectory (mesh_input)
add_subdirectory (mesh_output)
#add_subdirectory (operators)



# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d pditi/intrepid)
   set(UT_DIRS ${UT_DIRS}  -d pditi/aliti)
   set(UT_DIRS ${UT_DIRS}  -d pditi/pdneigh/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d pditi/utilities/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d pditi/mesh_input/quick_grid/unit_test)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
   set(SCRIPT_DIRS ${SCRIPT_DIRS}  -d scripts)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${SCRIPT_DIRS} )
endif (NOT ${buildInSource})

