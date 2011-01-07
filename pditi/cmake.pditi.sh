#! /bin/bash

function run_cmake {


	local BUILD_TYPE=$1

	cmake \
	 -D CMAKE_BUILD_TYPE:STRING=$BUILD_TYPE \
	 -D FIND_LIBRARY_USE_LIB64_PATHS:BOOL=ON \
	 -D TPL_ENABLE_MPI:BOOL=ON  \
	     -D MPI_BASE_DIR:PATH=$MPICH2_HOME \
	 -D CMAKE_LIBRARY_PATH:PATH=/usr/lib64 \
	 -D CMAKE_CXX_COMPILER:PATH=$CXX \
	 -D CMAKE_C_COMPILER:FILEPATH=$CC \
	 -D CMAKE_Fortran_COMPILER:FILEPATH=$F77 \
	 -D TPL_ENABLE_BLAS:STRING='ON' \
	 -D TRILINOS_DIR=$JAM_TRILINOS_NIGHTLY/$BUILD_TYPE \
	 -D PERIDIGM_MAT_UTILS_SRC=$PERIDIGM_HOME/PdMaterialUtilities \
	 -D PERIDIGM_MAT_UTILS_LIB=$PERIDIGM_HOME_BUILD/PdMaterialUtilities \
	 -D PERIDIGM_IO_SRC=$PERIDIGM_HOME/io \
	 -D PERIDIGM_IO_LIB=$PERIDIGM_HOME_BUILD/io \
	 -G"Eclipse CDT4 - Unix Makefiles" \
	 -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
	 -D ECLIPSE_CDT4_GENERATE_SOURCE_PROJECT=TRUE \
	 ../pditi-trunk

}

function run_cmake_home {

	VTK_LIB=/usr/local/lib/vtk-5.6
	MPICH2_HOME=/usr/local/mpich2-1.2.1-gcc4.4.1
	F77=$MPICH2_HOME/bin/mpif77
	CC=$MPICH2_HOME/bin/mpicc
	CXX=$MPICH2_HOME/bin/mpicxx
	export VTK_LIB MPICH2_HOME F77 CC CXX 

	cmake \
	 -D TPL_ENABLE_MPI:BOOL=ON  \
	     -D MPI_BASE_DIR:PATH=$MPICH2_HOME \
	 -D CMAKE_LIBRARY_PATH:PATH=/usr/local/gcc-4.4.1/lib64 \
	 -D CMAKE_CXX_COMPILER:PATH=$CXX \
	 -D CMAKE_C_COMPILER:FILEPATH=$CC \
	 -D CMAKE_Fortran_COMPILER:FILEPATH=$F77 \
	 -D TPL_ENABLE_BLAS:STRING='ON' \
	 -D BLAS_LIBRARY_DIRS:PATH="/usr/local/acml4.4.0/gfortran64/lib;/usr/local/acml4.4.0/gfortran64/include" \
	 -D BLAS_LIBRARY_NAMES:STRING=acml \
	 -D LAPACK_LIBRARY_DIRS:PATH="/usr/local/acml4.4.0/gfortran64/lib;/usr/local/acml4.4.0/gfortran64/include" \
	 -D LAPACK_LIBRARY_NAMES:STRING=acml \
	 -D TRILINOS_DIR=$HOME/trilinos-10.0.4/trilinos-10.0.4-installed \
	 -D PERIDIGM_MAT_UTILS_SRC=/home/awesome/c++/eclipseProjects/peridigm/PdMaterialUtilities \
	 -D PERIDIGM_MAT_UTILS_LIB=/home/awesome/c++/eclipseProjects/peridigm/PdMaterialUtilities.build \
	 -D PERIDIGM_IO_SRC=/home/awesome/c++/eclipseProjects/peridigm/io \
	 -D PERIDIGM_IO_LIB=/home/awesome/c++/eclipseProjects/peridigm/io.build \
	 -G"Eclipse CDT4 - Unix Makefiles" \
	 -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
	 -D ECLIPSE_CDT4_GENERATE_SOURCE_PROJECT=TRUE \
	 ../pditi

}


# This causes build type to default to DEBUG
# Otherwise, build type is set my first argument
INPUT_BUILD_TYPE=$1
BUILD_TYPE=${INPUT_BUILD_TYPE:-DEBUG}
BUILD_OPTIONS="DEBUG RELEASE"

# invoke cmake
if [ -f CMakeCache.txt ]; then
	echo Removing \'CMakeCache.txt\'
	rm CMakeCache.txt
fi
run_cmake $BUILD_TYPE

