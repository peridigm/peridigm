#
# CMakeLists.txt file for PdITI -- Peridynamics::Implicit
#

set(HOME_BUILD 0)
set(HOME /home/wow)
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
endif()

set(PERIDIGM_LIBS PdZoltan PdQuickGrid  PdNeighborhood PdVTK PdMaterialUtilities)

# Project directories
add_subdirectory (aliti)
add_subdirectory (intrepid)
add_subdirectory (intrepid/unitTest)
add_subdirectory (operators)
add_subdirectory (pdneigh)
add_subdirectory (pdneigh/unit_test)
add_subdirectory (operators/unitTest)
add_subdirectory (operators/unitTest/implicitDynamicsDemo_bar)
add_subdirectory (operators/unitTest/crackOpeningDemo)
add_subdirectory (operators/unitTest/eigen1D)
add_subdirectory (operators/unitTest/twoPointPlasticity)


#string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
#if (NOT ${buildInSource})
#   message ("Creating symbolic links for out-of-source build\n")
#   set(UT_DIRS ${UT_DIRS}  -d intrepid)
#   set(UT_DIRS ${UT_DIRS}  -d aliti)
#   set(UT_DIRS ${UT_DIRS}  -d operators/unitTest)
#   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
#endif (NOT ${buildInSource})

