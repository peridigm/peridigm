#add_subdirectory (unit_test)

# Collection of libs that we need to link with
set(PD_RKPM_SOURCES
    rkpmshapefunction.cxx
)

add_library(PdRKPMUtilities ${PD_RKPM_SOURCES})
IF (INSTALL_PERIDIGM)
   install(TARGETS PdRKPMUtilities EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()
