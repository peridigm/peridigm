#
# CMakeLists.txt file for Peridigm "states" Directory
#
# To run cmake on this file user must define following cmake variables:
#  TRILINOS_DIR -- this is user installation of trilinos
#  CMAKE_CXX_COMPILER -- path to mpi compiler
#  PTHREAD_LIB -- FILEPATH to library eg /usr/lib64/libpthread_nonshared.a

#cmake \
# -D CMAKE_CXX_COMPILER=/usr/bin/mpic++ \
# -D vtk-5.4-InstallPath=/usr/local \
# -G"Eclipse CDT4 - Unix Makefiles" \
# -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
# -D ECLIPSE_CDT4_GENERATE_SOURCE_PROJECT=TRUE \
# ../states \

Cmake_minimum_required (VERSION 2.6)

enable_testing ()

# Basic settings for PdStates project
project (PdStates)

add_library(PdStates SHARED PdStates.cxx)

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

add_subdirectory(test)