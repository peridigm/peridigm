# Library for Peridynamics neighborhood creation
add_subdirectory(unit_test)

# include this path
add_library(PdNeigh ../Peridigm_VTKSearchTree.cpp ../Peridigm_JAMSearchTree.cpp ../Peridigm_ZoltanSearchTree.cpp NeighborhoodList.cxx PdZoltan.cxx BondFilter.cxx OverlapDistributor.cxx)
