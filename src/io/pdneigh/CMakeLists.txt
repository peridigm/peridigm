# Library for Peridynamics neighborhood creation
add_subdirectory(unit_test)

# include this path
add_library(PdNeigh ../Peridigm_JAMSearchTree.cpp ../Peridigm_ZoltanSearchTree.cpp NeighborhoodList.cxx PdZoltan.cxx BondFilter.cxx OverlapDistributor.cxx)

IF (INSTALL_PERIDIGM)
   install(TARGETS PdNeigh EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()