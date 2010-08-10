Cmake_minimum_required (VERSION 2.6)

enable_testing ()

# Basic settings for unit testing PdStates project
project (utPdStates)

# MPI configuration
# use the FindMPI module to find a version of MPI
# expose CMAKE_CXX_COMPILER so thaet the user can change it to the MPI compiler
mark_as_advanced (CLEAR CMAKE_CXX_COMPILER)
mark_as_advanced (CLEAR CMAKE_CXX_FLAGS)
if (NOT CMAKE_CXX_FLAGS)
  set (CMAKE_CXX_FLAGS "-DMPICH_IGNORE_CXX_SEEK" CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif (NOT CMAKE_CXX_FLAGS)

# Boost configuration 
find_package (Boost 1.37.0)
if (Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif ()

set(libPdStates_DIR /home/jamitch/peridigm_projects/eclipse-vtk/states.build)
set(Boost_UNIT_TEST_FRAMEWORK_LIB boost_unit_test_framework-gcc41-mt)
include_directories(../ ${Boost_INCLUDE_DIRS})
link_directories(${libPdStates_DIR} ${Boost_LIB_DIRS})

# tests
add_executable(utFieldData utFieldSpec.cxx)
target_link_libraries(utFieldData ${Boost_UNIT_TEST_FRAMEWORK_LIB} PdStates)