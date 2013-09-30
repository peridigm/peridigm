#
# CMakeLists.txt file for Peridigm
#

Cmake_minimum_required (VERSION 2.6)

enable_testing ()
project(Peridigm)

# MPI configuration
# use the FindMPI module to find a version of MPI
# expose CMAKE_CXX_COMPILER so that the user can change it to the MPI compiler

mark_as_advanced (CLEAR CMAKE_CXX_COMPILER)

mark_as_advanced (CLEAR CMAKE_CXX_FLAGS)
if (NOT CMAKE_CXX_FLAGS)
  set (CMAKE_CXX_FLAGS "-DMPICH_IGNORE_CXX_SEEK" CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif (NOT CMAKE_CXX_FLAGS)

#
# Boost configuration
# Questions?  See http://www.cmake.org/cmake/help/cmake2.6docs.html#module:FindBoost
#
find_package (Boost 1.37.0 COMPONENTS unit_test_framework system filesystem thread regex)
if (Boost_FOUND)
  message ("\nBoost installation:")
  message ("   Version: ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
  message ("   Include Directory: ${Boost_INCLUDE_DIRS}")
  message ("   Library Directory: ${Boost_LIBRARY_DIRS}")
  message ("   The following libraries are included in Boost_LIBRARIES:")
  foreach (lib ${Boost_LIBRARIES})
    message ("      ${lib}")
  endforeach (lib)
else ()
  message ("\nBoost NOT FOUND!  (Try defining BOOST_ROOT)\n")
endif ()

#
# Blas and Lapack Libraries  
#
find_library(Blas_LIBRARY  
   NAMES blas  
   PATHS ${BLAS_LIBRARY_DIRS}  
)  
find_library(Lapack_LIBRARY  
   NAMES lapack  
   PATHS ${LAPACK_LIBRARY_DIRS}  
)  
set(BlasLapack_Libraries ${blas} ${lapack})  

#
# DAKOTA configuration
#
option (USE_DAKOTA 
   "Enable tests and functionality that depend on a DAKOTA installation." 
   OFF
)

#
# Trilinos configuration
#
set(TRILINOS_DIR CACHE PATH "Path to Trilinos directory")
SET(CMAKE_PREFIX_PATH ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)
MESSAGE("Trilinos installation:")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
# Hack for now until Trilinos exports cmake variable for binaries directory
SET(Trilinos_BINARY_PATH ${Trilinos_LIBRARY_DIRS}/../bin)
# verbose output
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("")

# Check for netcdf
LIST(FIND Trilinos_TPL_LIST Netcdf Netcdf_Package_Index)
IF(Netcdf_Package_Index LESS 0)
  MESSAGE("\n\nError:  Netcdf NOT found.  Netcdf is a REQUIRED Trilinos TPL.\n\n")
ENDIF()

#
# Enable performance testing
#
IF(TEST_PERFORMANCE)
  MESSAGE("-- Performance testing is enabled.\n")
ELSE()
  MESSAGE("-- Performance testing is NOT enabled.\n")
ENDIF()

#
# Optional dependency on Laboratory for Computational Mechanics (LCM)
#
IF(USE_LCM)
  SET(CMAKE_PREFIX_PATH ${ALBANY_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Albany REQUIRED)
  MESSAGE("Found Albany!  Here are the details: ")
  MESSAGE("   Albany_DIR = ${Albany_DIR}")
  MESSAGE("   Albany_LIBRARIES = ${Albany_LIBRARIES}")
  MESSAGE("   Albany_INCLUDE_DIRS = ${Albany_INCLUDE_DIRS}")
  MESSAGE("   Albany_LIBRARY_DIRS = ${Albany_LIBRARY_DIRS}")
  MESSAGE("End of Albany details\n")
  MESSAGE("-- LCM is enabled, compiling with -DPERIDIGM_LCM")
  MESSAGE("")
  ADD_DEFINITIONS(-DPERIDIGM_LCM)
  SET(PERIDIGM_LCM TRUE)
  include_directories(${Albany_INCLUDE_DIRS})
  link_directories(${Albany_LIBRARY_DIRS})
  set(LCM_LIBRARY ${Albany_LIBRARIES})
ELSE()
  MESSAGE("-- LCM is NOT enabled.\n")
  SET(PERIDIGM_LCM FALSE)
  SET(LCM_LIBRARY)
ENDIF()

#
# Enable PALS development
#
IF(USE_PALS)
  MESSAGE("-- PALS is enabled, compiling with -DPERIDIGM_PALS.\n")
  ADD_DEFINITIONS(-DPERIDIGM_PALS)
  SET(PERIDIGM_PALS TRUE)
ELSE()
  MESSAGE("-- PALS is NOT enabled.\n")
  SET(PERIDIGM_PALS FALSE)
ENDIF()

#
# Enable CJL development features
#
IF(USE_CJL)
  MESSAGE("-- CJL is enabled, compiling with -DPERIDIGM_CJL.\n")
  ADD_DEFINITIONS(-DPERIDIGM_CJL)
  SET(PERIDIGM_CJL TRUE)
ELSE()
  MESSAGE("-- CJL is NOT enabled.\n")
  SET(PERIDIGM_CJL FALSE)
  SET(CJL_LIBRARY)
ENDIF()

#
# Enable Sandia internal development
#
IF(USE_SANDIA_INTERNAL)
  MESSAGE("-- SANDIA_INTERNAL is enabled, compiling with -DPERIDIGM_SANDIA_INTERNAL.\n")
  ADD_DEFINITIONS(-DPERIDIGM_SANDIA_INTERNAL)
  SET(PERIDIGM_SANDIA_INTERNAL TRUE)
ELSE()
  MESSAGE("-- SANDIA_INTERNAL is NOT enabled.\n")
  SET(PERIDIGM_SANDIA_INTERNAL FALSE)
  SET(SANDIA_INTERNAL_LIBRARY)
ENDIF()

# Optional Installation helpers
# Note that some of this functionality depends on CMAKE > 2.8.8
SET(INSTALL_PERIDIGM FALSE)
IF(ENABLE_INSTALL) 
  IF(${CMAKE_VERSION} VERSION_GREATER 2.8.8 OR ${CMAKE_VERSION} VERSION_EQUAL 2.8.8)
    SET(INSTALL_PERIDIGM TRUE)
    MESSAGE("")
    MESSAGE("Installation files for Peridigm will be created")
    MESSAGE("-- install location = ${CMAKE_INSTALL_PREFIX}")
    MESSAGE("")
  ELSE()
    MESSAGE("")
    MESSAGE("Installation has been requested, but cannot continue")
    MESSAGE("-- Please upgrade CMAKE to version >= 2.8.8")
    MESSAGE("")
  ENDIF()
ENDIF()

IF(INSTALL_PERIDIGM)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR lib)
  include(CMakePackageConfigHelpers)
ENDIF()

# For some unit tests, we want to link only to specfic libraries
set(Epetra_LIBRARY epetra)
set(Zoltan_LIBRARY zoltan)

# other required libraries
set(MESH_INPUT_LIBS QuickGrid)
set(PARSER_LIBS muParser)
set(PDNEIGH_LIBS PdNeigh)
set(UTILITIES_LIBS Utilities)
set(PdMaterialUtilitiesLib PdMaterialUtilities)

# include for all directories
# Collection of include and lib directories
set(INCLUDE_DIRS 
${Trilinos_INCLUDE_DIRS} 
${Trilinos_TPL_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS} 
${CMAKE_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/src
${CMAKE_SOURCE_DIR}/src/core
${CMAKE_SOURCE_DIR}/src/io
${CMAKE_SOURCE_DIR}/src/io/discretization
${CMAKE_SOURCE_DIR}/src/io/mesh_input
${CMAKE_SOURCE_DIR}/src/io/mesh_input/quick_grid
${CMAKE_SOURCE_DIR}/src/io/mesh_output
${CMAKE_SOURCE_DIR}/src/io/bond_volume
${CMAKE_SOURCE_DIR}/src/io/utilities
${CMAKE_SOURCE_DIR}/src/io/pdneigh
${CMAKE_SOURCE_DIR}/src/contact
${CMAKE_SOURCE_DIR}/src/evaluators
${CMAKE_SOURCE_DIR}/src/materials
${CMAKE_SOURCE_DIR}/src/damage
)

include_directories(${INCLUDE_DIRS})

#
# Some convenience definitions
#
set(COMPUTE_DIR        ${CMAKE_SOURCE_DIR}/src/compute)
set(CONTACT_DIR        ${CMAKE_SOURCE_DIR}/src/contact)
set(CORE_DIR           ${CMAKE_SOURCE_DIR}/src/core)
set(DISCRETIZATION_DIR ${CMAKE_SOURCE_DIR}/src/io/discretization)
set(EVALUATORS_DIR     ${CMAKE_SOURCE_DIR}/src/evaluators)
set(IO_DIR             ${CMAKE_SOURCE_DIR}/src/io)
set(MATERIALS_DIR      ${CMAKE_SOURCE_DIR}/src/materials)
set(DAMAGE_DIR         ${CMAKE_SOURCE_DIR}/src/damage)
set(MESH_OUTPUT_DIR    ${CMAKE_SOURCE_DIR}/src/io/mesh_output)
set(QUICK_GRID_DIR     ${CMAKE_SOURCE_DIR}/src/io/mesh_input/quick_grid)
set(PDNEIGH_DIR        ${CMAKE_SOURCE_DIR}/src/io/pdneigh)
set(UTILITIES_DIR      ${CMAKE_SOURCE_DIR}/src/io/utilities)

# link directories 
set(LIB_DIRS 
  ${Boost_LIBRARY_DIRS})
link_directories(${LIB_DIRS})

set (REQUIRED_LIBS
  ${BlasLapack_Libraries}  
)

set (UT_REQUIRED_LIBS
  ${BlasLapack_Libraries}  
  ${Boost_LIBRARIES}
)

add_subdirectory (src)
add_subdirectory (test/regression)
add_subdirectory (test/verification)
add_subdirectory (test/performance)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d test/regression)
   set(UT_DIRS ${UT_DIRS}  -d test/verification)
   set(UT_DIRS ${UT_DIRS}  -d test/performance)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
   set(SCRIPT_DIRS ${SCRIPT_DIRS}  -d scripts)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${SCRIPT_DIRS} )
endif (NOT ${buildInSource})
# Set symbolic link to needed seacas tools in scripts dir
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/exodiff ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/epu ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/conjoin ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/exotxt ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/decomp ${CMAKE_BINARY_DIR}/scripts)
