[NetCDF](https://github.com/Unidata/netcdf-c/releases) is required by the Trilinos SEACAS package. 
Prior to compiling NetCDF, it is recommended that you modify the file netcdf.h to better support 
large-scale Peridigm simulations, as described below.

An example build process for NetCDF:


````
# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77
````

Modify the following `#define` statements in the `netcdf.h` file.  Change the values to match what is given below.

````
#define NC_MAX_DIMS 65536                                                                                                    
#define NC_MAX_ATTRS 8192                                                                                      
#define NC_MAX_VARS 524288                                                                                                    
#define NC_MAX_NAME 256                                                                                                      
#define NC_MAX_VAR_DIMS 8   
````

Alternatively, you can install [this fork](https://github.com/johntfoster/netcdf-c) of NetCDF where the changes have already
been made.  This version is known to work with Trilinos and Peridigm.

````
# Configure NetCDF
./configure --prefix=/usr/local/netcdf --enable-parallel--enable-netcdf-4 --disable-v2 --disable-fsync --disable-dap
````

````
# Make, test, and install NetCDF
make -j 8
make check
make install
````
