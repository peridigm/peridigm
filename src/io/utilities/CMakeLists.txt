#
# private utilities used in 'pdneigh', 'mesh_input', 'mesh_output', and 'pditi'
#

add_subdirectory(unit_test)
#
# Utilities library
#
add_library(Utilities Sortable.cxx Peridigm_Memstat.cpp)

IF (INSTALL_PERIDIGM)
   install(TARGETS Utilities EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()