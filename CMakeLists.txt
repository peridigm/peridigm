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
find_package (Boost 1.37.0 COMPONENTS unit_test_framework)
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
# VTK configuration
# Questions? See http://www.cmake.org/cmake/help/cmake2.6docs.html#module:FindVTK
#
find_package (VTK ${VTK_VERSION} NAMES VTK PATHS ${VTK_ROOT})
if (VTK_FOUND)
  message ("\nVTK installation:")
  message ("   Version: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  message ("   Include Directories: ${VTK_INCLUDE_DIRS}")
  message ("   Library Directory: ${VTK_LIBRARY_DIRS}")
  message ("   VTK Kits: ${VTK_KITS}")
  message ("   VTK_USE_FILE: ${VTK_USE_FILE}\n")
  include (${VTK_USE_FILE})
else ()
  message ("\nVTK NOT FOUND!  (Try defining VTK_DIR)\n")
endif ()
set( VTK_LIBS vtkIO vtkGraphics vtkParallel )

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
set(BlasLapack_Libraries blas lapack)  

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
# MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
# MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
# MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("")

# Check for netcdf
LIST(FIND Trilinos_TPL_LIST Netcdf Netcdf_Package_Index)
IF(Netcdf_Package_Index GREATER -1)
  MESSAGE("-- Found Netcdf")
ELSE()
  MESSAGE("Error:  Netcdf NOT found")
ENDIF()

# For some unit tests, we want to link only to specfic libraries
set(Epetra_LIBRARY epetra )
set(Zoltan_LIBRARY zoltan )

# other required libraries
set(MESH_INPUT_LIBS QuickGrid)

# other required libraries
set(MESH_OUTPUT_LIBS PdVTK)

# other required libraries
set(PARSER_LIBS muParser)

# load balancing and neighborhood creation
set(PDNEIGH_LIBS PdNeigh)
set(UTILITIES_LIBS Utilities)

# other required libraries
set(PdMaterialUtilitiesLib PdMaterialUtilities)

# include for all directories
# Collection of include and lib directories
set(INCLUDE_DIRS 
${Trilinos_INCLUDE_DIRS} 
${Trilinos_TPL_INCLUDE_DIRS}
${VTK_INCLUDE_DIRS} 
${Boost_INCLUDE_DIRS} 
${CMAKE_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/src
${CMAKE_SOURCE_DIR}/src/core
${CMAKE_SOURCE_DIR}/src/io
${CMAKE_SOURCE_DIR}/src/io/discretization
${CMAKE_SOURCE_DIR}/src/io/mesh_input
${CMAKE_SOURCE_DIR}/src/io/mesh_output
${CMAKE_SOURCE_DIR}/src/io/bond_volume
${CMAKE_SOURCE_DIR}/src/io/utilities
${CMAKE_SOURCE_DIR}/src/io/pdneigh
${CMAKE_SOURCE_DIR}/src/contact
${CMAKE_SOURCE_DIR}/src/evaluators
${CMAKE_SOURCE_DIR}/src/materials
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
set(MESH_OUTPUT_DIR    ${CMAKE_SOURCE_DIR}/src/io/mesh_output)

# link directories 
set(LIB_DIRS 
  ${VTK_LIBRARY_DIRS} 
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

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d test/regression)
   set(UT_DIRS ${UT_DIRS}  -d test/verification)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
   set(SCRIPT_DIRS ${SCRIPT_DIRS}  -d scripts)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${SCRIPT_DIRS} )
endif (NOT ${buildInSource})
# Set symbolic link to needed seacas tools in scripts dir
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/exodiff ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/epu ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/conjoin ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/exotxt ${CMAKE_BINARY_DIR}/scripts)
