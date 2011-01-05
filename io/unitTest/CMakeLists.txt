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
# ../peridigm.votd 


# convenience for this directory
set(IO_LIBS PdZoltan PdQuickGrid PdNeighborhood  PdVTK)

# tests
add_executable(utPdQuickGridHorizon utPdQuickGridHorizon.cxx)
target_link_libraries(utPdQuickGridHorizon ${IO_LIBS} ${VTK_LIBS} ${Epetra_LIBRARY} ${UT_REQUIRED_LIBS} )
#
add_executable(utPdQ2CylinderRingHorizon utPdQ2CylinderRingHorizon.cxx)
target_link_libraries(utPdQ2CylinderRingHorizon ${IO_LIBS} ${VTK_LIBS} ${Epetra_LIBRARY} ${UT_REQUIRED_LIBS})
#
add_executable(utPdSmallMeshCylinder utPdSmallMeshCylinder.cxx)
target_link_libraries(utPdSmallMeshCylinder ${IO_LIBS} ${VTK_LIBS} ${Epetra_LIBRARY} ${UT_REQUIRED_LIBS})
#
add_executable(utPdSmallMeshCylinder_np4 utPdSmallMeshCylinder_np4.cxx)
target_link_libraries(utPdSmallMeshCylinder_np4 ${IO_LIBS} ${VTK_LIBS}  ${Epetra_LIBRARY} ${UT_REQUIRED_LIBS})
#
add_executable(utPdQuickGrid_loadBal_np2_3x1x1 utPdQuickGrid_loadBal_np2_3x1x1.cxx)
target_link_libraries(utPdQuickGrid_loadBal_np2_3x1x1  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#
add_executable(utPdQuickGrid_loadBal_np2_4x4x4 utPdQuickGrid_loadBal_np2_4x4x4.cxx)
target_link_libraries(utPdQuickGrid_loadBal_np2_4x4x4  ${IO_LIBS}  ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#
add_executable(utPdQuickGrid_loadBal_np8_4x4x4 utPdQuickGrid_loadBal_np8_4x4x4.cxx)
target_link_libraries(utPdQuickGrid_loadBal_np8_4x4x4  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#
add_executable( utPdQuickGridMPI_np3  utPdQuickGridMPI_np3.cxx)
target_link_libraries( utPdQuickGridMPI_np3 ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#
add_executable(utPdQuickGridMPI_np2_3x3x2 utPdQuickGridMPI_np2_3x3x2.cxx)
target_link_libraries(utPdQuickGridMPI_np2_3x3x2  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#
add_executable(utZoltanQuery_pointSizeInBytes utZoltanQuery_pointSizeInBytes.cxx)
target_link_libraries(utZoltanQuery_pointSizeInBytes  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})
#
add_executable(utZoltanQuery_packUnpackPointsMultiFunction utZoltanQuery_packUnpackPointsMultiFunction.cxx)
target_link_libraries(utZoltanQuery_packUnpackPointsMultiFunction  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPdVtkKdTree utPdVtkKdTree.cxx)
target_link_libraries(utPdVtkKdTree ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPdNeighborhoodSort utPdNeighborhoodSort.cxx)
target_link_libraries(utPdNeighborhoodSort  ${IO_LIBS} ${VTK_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPdNeighborhoodProcessorFrame utPdNeighborhoodProcessorFrame.cxx)
target_link_libraries(utPdNeighborhoodProcessorFrame  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPd_frameSet_2x2x1_np4 utPd_frameSet_2x2x1_np4.cxx)
target_link_libraries(utPd_frameSet_2x2x1_np4  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPdNeighborhoodList utPdNeighborhoodList.cxx)
target_link_libraries(utPdNeighborhoodList  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPd_axisAlignedBoundaryPoints utPd_axisAlignedBoundaryPoints.cxx)
target_link_libraries(utPd_axisAlignedBoundaryPoints  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPdQuickGrid_solidCylinder_np1 utPdQuickGrid_solidCylinder_np1.cxx)
target_link_libraries(utPdQuickGrid_solidCylinder_np1  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utPdQuickGrid_solidCylinder_np2 utPdQuickGrid_solidCylinder_np2.cxx)
target_link_libraries(utPdQuickGrid_solidCylinder_np2  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

#
add_executable(utReloadBalance_np4 utReloadBalance_np4.cxx)
target_link_libraries(utReloadBalance_np4  ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

add_executable(utVTK_ioExample_np4 utVTK_ioExample_np4.cxx)
target_link_libraries(utVTK_ioExample_np4  ${PdMaterialUtilitiesLib} ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

add_executable(utPdVtkPlane utPdVtkPlane.cxx)
target_link_libraries(utPdVtkPlane  ${PdMaterialUtilitiesLib} ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

add_executable(utFinitePlaneFilter utFinitePlaneFilter.cxx)
target_link_libraries(utFinitePlaneFilter ${PdMaterialUtilitiesLib} ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

add_executable(utPd_Y-Z_crack utPd_Y-Z_crack.cxx)
target_link_libraries(utPd_Y-Z_crack ${PdMaterialUtilitiesLib} ${IO_LIBS} ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})


