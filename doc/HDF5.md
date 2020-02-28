## HDF5

[HDF5](https://www.hdfgroup.org/solutions/hdf5) is required by NetCDF and the SEACAS Trilinos package. HDF5 should be configured with the `--enable-parallel` option.

An example build process for [HDF5](https://www.hdfgroup.org/solutions/hdf5):

````
# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77

# Configure HDF5
./configure --prefix=/usr/local/HDF5 --enable-parallel

# Make and install HDF5
make -j 4
make install
````
