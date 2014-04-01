add_subdirectory (unit_test)

# Collection of libs that we need to link with
set(PD_MATERIAL_SOURCES
    ../core/Peridigm_InfluenceFunction.cpp
    elastic.cxx
    elastic_plastic.cxx
    elastic_plastic_hardening.cxx
    viscoelastic.cxx
    correspondence.cxx
    elastic_correspondence.cxx
    elastic_plastic_correspondence.cxx
    isotropic_hardening_correspondence.cxx
    viscoplastic_needleman_correspondence.cxx
    material_utilities.cxx
)

# Optional source files for Sandia internal development
IF(PERIDIGM_SANDIA_INTERNAL)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} partial_stress.cxx)
ENDIF()
IF(PERIDIGM_PALS)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} pals.cxx)
ENDIF()
IF(PERIDIGM_PV)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} elastic_pv.cxx correspondence_pv.cxx)
ENDIF()

# Optional source files for CJL development
IF(PERIDIGM_CJL)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} lammi_concrete_model.cxx)
ENDIF()

# PdMaterialUtilities library
add_library(PdMaterialUtilities ${PD_MATERIAL_SOURCES})

IF (INSTALL_PERIDIGM)
   install(TARGETS PdMaterialUtilities EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()
