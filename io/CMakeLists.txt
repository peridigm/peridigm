#
# Discretization I/O and load balancing for Peridigm
#

# VTK configuration
set(VTK_LIBS vtkIO vtkGraphics)

#PdQuickGrid library
add_library(PdQuickGrid PdQuickGrid.cxx PdQuickGridParallel.cxx Field.cxx)

# PdVTK library
add_library(PdVTK PdVTK.cxx Field.cxx)

# PdZoltan library
add_library(PdZoltan PdZoltan.cxx)

# Parallel Neighborhood Search library
add_library(PdNeighborhood PdNeighborhood.cxx PdBondFilter.cxx)



