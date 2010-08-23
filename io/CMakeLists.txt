#
# Discretization I/O and load balancing for Peridigm
#
#cmake \
# -D CMAKE_BUILD_TYPE:STRING=Debug \
# -D CMAKE_CXX_COMPILER:PATH=/usr/bin/mpic++ \
# -D VTK_DIR=/usr/local/vtk-install/lib/vtk-5.4 \
# -D TRILINOS_DIR=/home/jamitch/trilinos_votd/installed \
# -D BOOST_ROOT=/usr/include/boost \
# -D FIND_LIBRARY_USE_LIB64_PATHS:BOOL=ON \
# -D TPL_ENABLE_MPI:BOOL=ON  \
#     -D MPI_BASE_DIR:PATH=$MPICH2_HOME \
# -G"Eclipse CDT4 - Unix Makefiles" \
# -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
# -D ECLIPSE_CDT4_GENERATE_SOURCE_PROJECT=TRUE \
# ../peridigm 

# Collection of libs that we need to link with
include_directories(${TRILINOS_DIR}/include ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
link_directories(${TRILINOS_DIR}/lib ${VTK_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})

# VTK configuration
set(VTK_LIBS vtkIO vtkGraphics)

#PdQuickGrid library
add_library(PdQuickGrid PdQuickGrid.cxx PdQuickGridParallel.cxx Field.cxx)

# PdVTK library
add_library(PdVTK PdVTK.cxx Field.cxx)

# PdZoltan library
add_library(PdZoltan PdZoltan.cxx)

# Parallel Neighborhood Search library
add_library(PdNeighborhood PdNeighborhood.cxx)



