#
# CMakeLists.txt file for Peridigm
#

Cmake_minimum_required (VERSION 2.6)

enable_testing ()

# Basic settings for Peridigm project
project (Peridigm.votd)
set (EXECUTABLE_NAME Peridigm)

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
find_package (VTK)
if (VTK_FOUND)
  message ("\nVTK installation:")
  message ("   Version: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  message ("   Include Directories: ${VTK_INCLUDE_DIRS}")
  message ("   Library Directory: ${VTK_LIBRARY_DIRS}")
  message ("   VTK Kits: ${VTK_KITS}\n")
else ()
  message ("\nVTK NOT FOUND!  (Try defining VTK_DIR)\n")
endif ()


# Trilinos configuration
set (TRILINOS_DIR CACHE PATH "Path to Trilinos directory")

find_library(Phalanx_LIBRARY
  NAMES phalanx
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Stratimikos_LIBRARY
  NAMES stratimikos
  PATHS ${TRILINOS_DIR}/lib
)


find_library(StratimikosAmesos_LIBRARY
  NAMES stratimikosamesos
  PATHS ${TRILINOS_DIR}/lib
)

find_library(StratimikosIFPack_LIBRARY
  NAMES stratimikosifpack
  PATHS ${TRILINOS_DIR}/lib
)

find_library(StratimikosML_LIBRARY
  NAMES stratimikosml
  PATHS ${TRILINOS_DIR}/lib
)

find_library(StratimikosAztecOO_LIBRARY
  NAMES stratimikosaztecoo
  PATHS ${TRILINOS_DIR}/lib
)

find_library(ML_LIBRARY
  NAMES ml
  PATHS ${TRILINOS_DIR}/lib
)

find_library(IFPack_LIBRARY
  NAMES ifpack
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Amesos_LIBRARY
  NAMES amesos
  PATHS ${TRILINOS_DIR}/lib
)

find_library(AztecOO_LIBRARY
  NAMES aztecoo
  PATHS ${TRILINOS_DIR}/lib
)

find_library(ThyraEpetraExt_LIBRARY
  NAMES thyraepetraext
  PATHS ${TRILINOS_DIR}/lib
)

find_library(ThyraEpetra_LIBRARY
  NAMES thyraepetra
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Thyra_LIBRARY
  NAMES thyra
  PATHS ${TRILINOS_DIR}/lib
)

find_library(EpetraExt_LIBRARY
  NAMES epetraext
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Zoltan_LIBRARY
  NAMES zoltan
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Shards_LIBRARY
  NAMES shards
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Epetra_LIBRARY
  NAMES epetra
  PATHS ${TRILINOS_DIR}/lib
)

find_library(RTOp_LIBRARY
  NAMES rtop
  PATHS ${TRILINOS_DIR}/lib
)

find_library(Teuchos_LIBRARY
  NAMES teuchos
  PATHS ${TRILINOS_DIR}/lib
)


set(TRILINOS_LIBS
phalanx
stratimikos
stratimikosamesos
stratimikosifpack
stratimikosml
stratimikosaztecoo
ml
ifpack
amesos
aztecoo
thyraepetraext
thyraepetra
thyra
epetraext
zoltan
shards
epetra
rtop
teuchos
)

# other required libraries
# VTK libs
set(VTK_LIBS vtkIO vtkGraphics)

set(PdQuickGridLibs PdQuickGrid PdNeighborhood)

# other required libraries
set(PdVTKLibs PdVTK)

# other required libraries
set(PdZoltanLibs PdZoltan)

# other required libraries
set(PdMaterialUtilitiesLib PdMaterialUtilities)

# include for all directories
# Collection of include and lib directories
set(INCLUDE_DIRS 
$ENV{TRILINOS_DIR}/include  
${VTK_INCLUDE_DIRS} 
${Boost_INCLUDE_DIRS} 
${CMAKE_SOURCE_DIR}/io
${CMAKE_SOURCE_DIR}/PdMaterialUtilities
${CMAKE_SOURCE_DIR}/src
${CMAKE_SOURCE_DIR}/src/contact
${CMAKE_SOURCE_DIR}/src/evaluators
${CMAKE_SOURCE_DIR}/src/materials
)

include_directories(${INCLUDE_DIRS})

# link directories 
set(LIB_DIRS ${TRILINOS_DIR}/lib ${VTK_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})
link_directories(${LIB_DIRS})


set (REQUIRED_LIBS
liblapack.a
libblas.a
libmpich.a
${PTHREAD_LIB}
libgfortran.a
libm.a
)

add_subdirectory (io)
add_subdirectory (io/unitTest)
add_subdirectory (src)
add_subdirectory (PdMaterialUtilities)
add_subdirectory (test/unit)
add_subdirectory (test/regression)
add_subdirectory (test/verification)

# Create symbolic links to files in directories specified.  Does not recurse folders -- yet!
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d test/unit)
   set(UT_DIRS ${UT_DIRS}  -d test/regression)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/cmake_symlinks.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})

