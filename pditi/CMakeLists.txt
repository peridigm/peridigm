#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

Cmake_minimum_required (VERSION 2.8)

enable_testing ()

# Basic settings for PdITI project
project (PdITI)
set (EXECUTABLE_NAME PdITI)

set(HOME_BUILD 0)
set(HOME /home/awesome)
message("\nPdITI build from $ENV{HOME}\n")
if(EXISTS ${HOME})
    set(HOME_BUILD 1)
endif()


# MPI configuration
# use the FindMPI module to find a version of MPI
# expose CMAKE_CXX_COMPILER so that the user can change it to the MPI compiler
mark_as_advanced (CLEAR CMAKE_CXX_COMPILER)

mark_as_advanced (CLEAR CMAKE_CXX_FLAGS)
if (NOT CMAKE_CXX_FLAGS)
  set (CMAKE_CXX_FLAGS "-DMPICH_IGNORE_CXX_SEEK" CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif (NOT CMAKE_CXX_FLAGS)

if(HOME_BUILD)
# -D BoostUnitTest_Library:FILEPATH=/usr/local/boost-1.44/lib/libboost_unit_test_framework.so \
    set(BOOST_ROOT ${BOOST_ROOT})
    set(BOOST_VERSION 1.44)
    
    # Blas and Lapack Libraries
    find_library(ACML_LIBRARY
      NAMES acml
      PATHS ${BLAS_LIBRARY_DIRS}
    )
    set(BlasLapack_Libraries acml)
    
    # version of vtk
    set(VTK_VERSION ${VTK_VERSION})
    
else()
    # Initialize BOOST 
    set(BOOST_ROOT ${BOOST_ROOT})
    set(BOOST_VERSION 1.37)
    
    # Blas and Lapack Libraries
    find_library(Blas_LIBRARY
      NAMES blas
      PATHS ${BLAS_LIBRARY_DIRS}
    )
    find_library(Lapack_LIBRARY
      NAMES lapack
      PATHS ${LAPACK_LIBRARY_DIRS}
    )
    set(BlasLapack_Libraries blas lapack)
    
endif()

# find boost
find_package (Boost ${BOOST_VERSION} COMPONENTS unit_test_framework)
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

# VTK configuration
set(VTK_LIBS vtkIO vtkpng vtksys vtksqlite vtkNetCDF vtktiff vtkCommon vtkFiltering vtkalglib vtkDICOMParser vtkmetaio vtkGraphics)

find_package (VTK ${VTK_VERSION} COMPONENTS ${VTK_LIBS})
if (VTK_FOUND)
  message ("\nVTK installation:")
  message ("   Version: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  message ("   Include Directories: ${VTK_INCLUDE_DIRS}")
  message ("   Library Directory: ${VTK_LIBRARY_DIRS}")
  message ("   VTK Kits: ${VTK_KITS}\n")
else ()
  message ("\nVTK NOT FOUND!  (Try defining VTK_DIR)\n")
endif ()

#set(VTK_LIBS 
#vtkIO vtkzlib vtkpng vtkjpeg vtkexpat vtksys vtksqlite 
#vtkNetCDF vtktiff vtkCommon vtkFiltering vtkalglib 
#vtkDICOMParser vtkmetaio vtkGraphics)


# Trilinos configuration
set (TRILINOS_DIR CACHE PATH "Path to Trilinos installation directory")
set (PERIDIGM_IO_LIB CACHE PATH "Path to Peridigm io directory")
set (PERIDIGM_MAT_UTILS_LIB CACHE PATH "Path to Peridigm constitutive model directory")

find_library(Epetra_LIBRARY
  NAMES epetra
  PATHS ${TRILINOS_DIR}/lib
)
find_library(AztecOO_LIBRARY
  NAMES aztecoo
  PATHS ${TRILINOS_DIR}/lib
)
find_library(Teuchos_LIBRARY
  NAMES teuchos
  PATHS ${TRILINOS_DIR}/lib
)
find_library(Zoltan_LIBRARY
  NAMES zoltan
  PATHS ${TRILINOS_DIR}/lib
)
find_library(Intrepid_LIBRARY
  NAMES intrepid
  PATHS ${TRILINOS_DIR}/lib
)
find_library(Shards_LIBRARY
  NAMES shards
  PATHS ${TRILINOS_DIR}/lib
)
set(TRILINOS_LIBS epetra aztecoo teuchos zoltan intrepid shards)


# Required Peridigm libraries
find_library(PdQuickGrid_LIBRARY
  NAMES PdQuickGrid
  PATHS ${PERIDIGM_IO_LIB}
)
# Peridigm Material Utilities Library
find_library(PdMaterialUtilities_LIBRARY
  NAMES PdMaterialUtilities
  PATHS ${PERIDIGM_MAT_UTILS_LIB}
)
# PdITI IO VTK Library
find_library(PdVTK_LIBRARY
  NAMES PdVTK
  PATHS ${PERIDIGM_IO_LIB}
)
# Peridigm Neighborhood Library
find_library(PdNeighborhood_LIBRARY
  NAMES PdNeighborhood
  PATHS ${PERIDIGM_IO_LIB}
)
# Peridynamics Zoltan Library
find_library(PdZoltan_LIBRARY
  NAMES PdZoltan
  PATHS ${PERIDIGM_IO_LIB}
)
set(PERIDIGM_LIBS PdZoltan PdQuickGrid  PdNeighborhood PdVTK PdMaterialUtilities)

set (REQUIRED_LIBS
libmpich.a
${PTHREAD_LIB}
libgfortran.a
${BlasLapack_Libraries}
)

set (UT_REQUIRED_LIBS
libmpich.a
${PTHREAD_LIB}
libgfortran.a
${Boost_LIBRARIES}
${BlasLapack_Libraries}
)

# Collection of include and lib directories
set(INCLUDE_DIRS ${TRILINOS_DIR}/include  ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${PERIDIGM_IO_SRC} ../ ../.. ${PERIDIGM_MAT_UTILS_SRC})
include_directories(${INCLUDE_DIRS})

set(LIB_DIRS ${BLAS_LIBRARY_DIRS} ${LAPACK_LIBRARY_DIRS} ${TRILINOS_DIR}/lib ${VTK_LIBRARY_DIRS} ${PERIDIGM_IO_LIB} ../ ../.. ${PERIDIGM_MAT_UTILS_LIB})
link_directories(${LIB_DIRS})

# Project directories
add_subdirectory (aliti)
add_subdirectory (intrepid)
add_subdirectory (intrepid/unitTest)
add_subdirectory (operators)
add_subdirectory (operators/unitTest)
add_subdirectory (operators/unitTest/implicitDynamicsDemo_bar)
add_subdirectory (operators/unitTest/eigen1D)
add_subdirectory (operators/unitTest/twoPointPlasticity)


string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   message ("Creating symbolic links for out-of-source build\n")
   set(UT_DIRS ${UT_DIRS}  -d intrepid)
   set(UT_DIRS ${UT_DIRS}  -d aliti)
   set(UT_DIRS ${UT_DIRS}  -d operators/unitTest)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})

