add_subdirectory (unit_test)

# Collection of libs that we need to link with

# PdMaterialUtilities library
add_library(PdMaterialUtilities
  ../core/Peridigm_InfluenceFunction.cpp
  elastic.cxx
  thermoelastic.cxx
  elastic_plastic.cxx
  elastic_plastic_hardening.cxx
  viscoelastic.cxx
  material_utilities.cxx
)


