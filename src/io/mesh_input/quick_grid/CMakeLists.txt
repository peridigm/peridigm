#
# QuickGrid library
#

add_subdirectory(unit_test)
include_directories(${UTILITIES_DIR})
add_library(QuickGrid QuickGrid.cxx)

IF (INSTALL_PERIDIGM)
   install(TARGETS QuickGrid EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()
