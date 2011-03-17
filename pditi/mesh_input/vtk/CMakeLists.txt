#
# vtk demo executables
#

include_directories(${UTILITIES_DIR} ${PDNEIGH_DIR} ${MESH_INPUT_DIR} ${MESH_OUTPUT_DIR} ${TRILINOS_DIR}/include)

add_executable(vtkExodusReaderDemo vtkExodusReaderDemo.cxx)
target_link_libraries(vtkExodusReaderDemo QuickGrid PdVTK ${VTK_LIBS} ${UT_REQUIRED_LIBS})

add_executable(utVTK_ioExample_np4 utVTK_ioExample_np4.cxx)
target_link_libraries(utVTK_ioExample_np4 QuickGrid PdNeigh PdVTK ${VTK_LIBS} ${TRILINOS_LIBS} ${UT_REQUIRED_LIBS})

add_executable(utVTK_ParallelReaderDemo utVTK_ParallelReaderDemo.cxx)
target_link_libraries(utVTK_ParallelReaderDemo  ${VTK_LIBS} )

