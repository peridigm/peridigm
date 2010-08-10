#
# Discretization I/O and load balancing for Peridigm
#
#cmake \
# -D CMAKE_BUILD_TYPE:STRING=Debug \
# -D CMAKE_CXX_COMPILER:PATH=/usr/bin/mpic++ \
# -D VTK_DIR=/usr/local/vtk-install/lib/vtk-5.4 \
# -D TRILINOS_DIR=/home/jamitch/trilinos_votd/installed \
# -D BOOST_ROOT=/usr/include/boost \
# -D FIND_LIBRARY_USE_LIB64_PATHS:BOOL=ON \
# -D TPL_ENABLE_MPI:BOOL=ON  \
#     -D MPI_BASE_DIR:PATH=$MPICH2_HOME \
# -G"Eclipse CDT4 - Unix Makefiles" \
# -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
# -D ECLIPSE_CDT4_GENERATE_SOURCE_PROJECT=TRUE \
# ../peridigm 

# Collection of libs that we need to link with
include_directories(${TRILINOS_DIR}/include ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ../)
link_directories(${TRILINOS_DIR}/lib ${VTK_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS} ../)

# tests
add_executable(utPdQuickGridHorizon utPdQuickGridHorizon.cxx)
target_link_libraries(utPdQuickGridHorizon PdNeighborhood PdQuickGrid PdVTK ${VTK_LIBS} ${Boost_LIBRARIES} ${Epetra_LIBRARY} ${REQUIRED_LIBS} )
#
add_executable(utPdQ2CylinderRingHorizon utPdQ2CylinderRingHorizon.cxx)
target_link_libraries(utPdQ2CylinderRingHorizon PdNeighborhood PdQuickGrid PdVTK ${VTK_LIBS} ${Boost_LIBRARIES} ${Epetra_LIBRARY} ${REQUIRED_LIBS})
#
add_executable(utPdSmallMeshCylinder utPdSmallMeshCylinder.cxx)
target_link_libraries(utPdSmallMeshCylinder PdNeighborhood PdQuickGrid PdVTK ${VTK_LIBS} ${Boost_LIBRARIES} ${Epetra_LIBRARY} ${REQUIRED_LIBS})
#
add_executable(utPdSmallMeshCylinder_np4 utPdSmallMeshCylinder_np4.cxx)
target_link_libraries(utPdSmallMeshCylinder_np4 PdNeighborhood PdQuickGrid PdVTK ${VTK_LIBS} ${Boost_LIBRARIES} ${Epetra_LIBRARY} ${REQUIRED_LIBS})
#
add_executable(utPdQuickGrid_loadBal_np2_3x1x1 utPdQuickGrid_loadBal_np2_3x1x1.cxx)
target_link_libraries(utPdQuickGrid_loadBal_np2_3x1x1  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#
add_executable(utPdQuickGrid_loadBal_np2_4x4x4 utPdQuickGrid_loadBal_np2_4x4x4.cxx)
target_link_libraries(utPdQuickGrid_loadBal_np2_4x4x4  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#
add_executable(utPdQuickGrid_loadBal_np8_4x4x4 utPdQuickGrid_loadBal_np8_4x4x4.cxx)
target_link_libraries(utPdQuickGrid_loadBal_np8_4x4x4  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#
add_executable( utPdQuickGridMPI_np3  utPdQuickGridMPI_np3.cxx)
target_link_libraries( utPdQuickGridMPI_np3 PdZoltan PdQuickGrid PdNeighborhood  PdVTK  ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#
add_executable(utPdQuickGridMPI_np2_3x3x2 utPdQuickGridMPI_np2_3x3x2.cxx)
target_link_libraries(utPdQuickGridMPI_np2_3x3x2  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#
add_executable(utZoltanQuery_pointSizeInBytes utZoltanQuery_pointSizeInBytes.cxx)
target_link_libraries(utZoltanQuery_pointSizeInBytes  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})
#
add_executable(utZoltanQuery_packUnpackPointsMultiFunction utZoltanQuery_packUnpackPointsMultiFunction.cxx)
target_link_libraries(utZoltanQuery_packUnpackPointsMultiFunction  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPdVtkKdTree utPdVtkKdTree.cxx)
target_link_libraries(utPdVtkKdTree PdQuickGrid PdVTK PdZoltan ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPdNeighborhoodSort utPdNeighborhoodSort.cxx)
target_link_libraries(utPdNeighborhoodSort  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPdNeighborhoodProcessorFrame utPdNeighborhoodProcessorFrame.cxx)
target_link_libraries(utPdNeighborhoodProcessorFrame  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPd_frameSet_2x2x1_np4 utPd_frameSet_2x2x1_np4.cxx)
target_link_libraries(utPd_frameSet_2x2x1_np4  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPdNeighborhoodList utPdNeighborhoodList.cxx)
target_link_libraries(utPdNeighborhoodList  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPd_axisAlignedBoundaryPoints utPd_axisAlignedBoundaryPoints.cxx)
target_link_libraries(utPd_axisAlignedBoundaryPoints  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPdQuickGrid_solidCylinder_np1 utPdQuickGrid_solidCylinder_np1.cxx)
target_link_libraries(utPdQuickGrid_solidCylinder_np1  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utPdQuickGrid_solidCylinder_np2 utPdQuickGrid_solidCylinder_np2.cxx)
target_link_libraries(utPdQuickGrid_solidCylinder_np2  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

#
add_executable(utReloadBalance_np4 utReloadBalance_np4.cxx)
target_link_libraries(utReloadBalance_np4  PdZoltan PdQuickGrid PdNeighborhood  PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${Boost_LIBRARIES} ${REQUIRED_LIBS})

