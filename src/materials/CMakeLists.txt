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
    material_utilities.cxx
)

# Optional source files for Sandia internal development
IF(PERIDIGM_SANDIA_INTERNAL)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} elastic_partial_stress.cxx)
ENDIF()

# Optional source files for CJL development
IF(PERIDIGM_CJL)
  set(PD_MATERIAL_SOURCES ${PD_MATERIAL_SOURCES} lammi_concrete_model.cxx)
ENDIF()

# PdMaterialUtilities library
add_library(PdMaterialUtilities ${PD_MATERIAL_SOURCES})
