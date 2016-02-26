# Building Peridigm

Peridigm utilizes the CMake build system. It is recommended that Makefiles be created by running `cmake` from the command line, as opposed to using the `ccmake` GUI.

Below is an example CMake configuration script for Peridigm. Note that the option `-std=c++11` within the `CMAKE_CXX_FLAGS` list is specific to compilers that support C++11 features. A compiler that is C++11 compliant (e.g., GCC 4.7.2 or later) is required for recent versions of Trilinos.

````
rm -f CMakeCache.txt
 
cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D Trilinos_DIR:PATH=/usr/local/trilinos/lib/cmake/Trilinos/ \
-D CMAKE_C_COMPILER:STRING=/usr/local/bin/mpicc \
-D CMAKE_CXX_COMPILER:STRING=/usr/local/bin/mpicxx \
-D BOOST_ROOT=/usr/local/boost/ \
-D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++11 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" \
<path to Peridigm source directory>
````

Once Peridigm has been successfully configured, it can be compiled as follows:

````
make -j 4
````

The Peridigm test suite is run as follows:

````
make test
````
