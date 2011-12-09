#
# CMakeLists.txt file for io directory
#

# Project directories
add_subdirectory (bond_volume)
add_subdirectory (discretization)
add_subdirectory (mesh_input)
add_subdirectory (mesh_output)
add_subdirectory (pdneigh)
add_subdirectory (utilities)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d bond_volume/quick_grid/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d pdneigh/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d utilities/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d mesh_input/quick_grid/unit_test)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${IO_DIR} -b ${CMAKE_BINARY_DIR}/src/io ${UT_DIRS} )
endif (NOT ${buildInSource})
