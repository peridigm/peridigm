if(USE_IMPROVED_QUADRATURE)
  add_test (DiffusionImprovedQuadrature_np1 python ./np1/DiffusionImprovedQuadrature.py)
  add_test (DiffusionImprovedQuadrature_np4 python ./np4/DiffusionImprovedQuadrature.py)
endif()

