#Boost

Peridigm requires [Boost C++ Libraries](http://www.boost.org/), version 1.37 or later, including the `regex` and `unit_test`
compiled libraries. Boost installations on many systems include header files only.  This is not sufficient, the required libraries must be compiled and installed. To ensure proper execution of Peridigm and its unit tests, add the boost directory `<install_dir>/lib` to your `LD_LIBRARY_PATH` and/or `DYLD_LIBRARY_PATH` environment variables.

An example build process for Boost:

````
# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77
````

````
# Run the Boost bootstrap script
./bootstrap.sh
````

````
# Compile and install Boost using the Boost's bjam build system
./b2 install --prefix=/user/local/boost/
````
