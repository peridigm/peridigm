#!/bin/bash



mkdir build
pushd build

TRIRIONS_PATH=/usr/local/trilinos/lib/cmake/Trilinos/
C_COMPILER_PATH=/usr/bin/mpicc
CXX_COMPILER_PATH=/usr/bin/mpicxx
BOOST_PATH=/usr/local/boost/

cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D Trilinos_DIR:PATH=${TRIRIONS_PATH} \
-D CMAKE_C_COMPILER:STRING=${C_COMPILER_PATH} \
-D CMAKE_CXX_COMPILER:STRING=${CXX_COMPILER_PATH} \
-D BOOST_ROOT=${BOOST_PATH} \
-D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++14 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" \
..

make -j 4

popd
