add_subdirectory (unit_test)

# Collection of libs that we need to link with

# PdMaterialUtilities library
IF(NOT PERIDIGM_CJL)
  add_library(PdMaterialUtilities
    ../core/Peridigm_InfluenceFunction.cpp
    elastic.cxx
    elastic_plastic.cxx
    elastic_plastic_hardening.cxx
    viscoelastic.cxx
    material_utilities.cxx
  )
ELSE()
  add_library(PdMaterialUtilities
    ../core/Peridigm_InfluenceFunction.cpp
    elastic.cxx
    elastic_plastic.cxx
    elastic_plastic_hardening.cxx
    viscoelastic.cxx
    lammi_concrete_model.cxx
    material_utilities.cxx
  )
ENDIF()


