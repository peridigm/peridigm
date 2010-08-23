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
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Phalanx_LIBRARY})

find_library(Stratimikos_LIBRARY
  NAMES stratimikos
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Stratimikos_LIBRARY})

find_library(StratimikosAmesos_LIBRARY
  NAMES stratimikosamesos
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${StratimikosAmesos_LIBRARY})

find_library(StratimikosIFPack_LIBRARY
  NAMES stratimikosifpack
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${StratimikosIFPack_LIBRARY})

find_library(StratimikosML_LIBRARY
  NAMES stratimikosml
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${StratimikosML_LIBRARY})

find_library(StratimikosAztecOO_LIBRARY
  NAMES stratimikosaztecoo
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${StratimikosAztecOO_LIBRARY})

find_library(ML_LIBRARY
  NAMES ml
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${ML_LIBRARY})

find_library(IFPack_LIBRARY
  NAMES ifpack
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${IFPack_LIBRARY})

find_library(Amesos_LIBRARY
  NAMES amesos
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Amesos_LIBRARY})

find_library(AztecOO_LIBRARY
  NAMES aztecoo
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${AztecOO_LIBRARY})

find_library(ThyraEpetraExt_LIBRARY
  NAMES thyraepetraext
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${ThyraEpetraExt_LIBRARY})

find_library(ThyraEpetra_LIBRARY
  NAMES thyraepetra
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${ThyraEpetra_LIBRARY})

find_library(Thyra_LIBRARY
  NAMES thyra
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Thyra_LIBRARY})

find_library(EpetraExt_LIBRARY
  NAMES epetraext
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${EpetraExt_LIBRARY})

find_library(Zoltan_LIBRARY
  NAMES zoltan
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Zoltan_LIBRARY})

find_library(Shards_LIBRARY
  NAMES shards
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Shards_LIBRARY})

find_library(Epetra_LIBRARY
  NAMES epetra
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Epetra_LIBRARY})

find_library(RTOp_LIBRARY
  NAMES rtop
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${RTOp_LIBRARY})

find_library(Teuchos_LIBRARY
  NAMES teuchos
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS ${TRILINOS_LIBS} ${Teuchos_LIBRARY})

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
