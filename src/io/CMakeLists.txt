##
# CMakeLists.txt file for io directory
#

SET(IO_ROOT ${CMAKE_SOURCE_DIR}/src/io)
SET(UTILITIES_DIR ${IO_ROOT}/utilities)
SET(MESH_OUTPUT_DIR ${IO_ROOT}/mesh_output)
SET(PDNEIGH_DIR ${IO_ROOT}/pdneigh)
SET(MESH_INPUT_DIR ${IO_ROOT}/mesh_input)
SET(TRILINOS_LIBS ${Trilinos_LIBRARIES})


# Project directories
add_subdirectory (discretization)
add_subdirectory (mesh_input)
add_subdirectory (mesh_output)
add_subdirectory (pdneigh)
add_subdirectory (utilities)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d pdneigh/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d utilities/unit_test)
   set(UT_DIRS ${UT_DIRS}  -d mesh_input/quick_grid/unit_test)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${IO_ROOT} -b ${CMAKE_BINARY_DIR}/src/io ${UT_DIRS} )
endif (NOT ${buildInSource})



