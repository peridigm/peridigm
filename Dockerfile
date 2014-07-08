FROM johntfoster/trilinos

MAINTAINER John Foster <johntfosterjr@gmail.com>

ENV HOME /root

RUN apt-get -yq install libboost1.49-all-dev
RUN apt-get -yq install openmpi-bin

RUN mkdir -p /peridigm
ADD CMakeLists.txt /peridigm/CMakeLists.txt
ADD src /peridigm/src
ADD test /peridigm/test
ADD examples /peridigm/examples
ADD misc /peridigm/misc
ADD scripts /peridigm/scripts

#Build Peridigm
RUN mkdir -p /peridigm/build
WORKDIR /peridigm/build/
RUN cmake \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D CMAKE_INSTALL_PREFIX:PATH=/usr/local/Peridigm \
    -D CMAKE_EXE_LINKER_FLAGS:STRING="" \
    -D CMAKE_C_FLAGS:STRING="-O3" \
    -D CMAKE_CXX_FLAGS:STRING="-O3" \
    -D TRILINOS_DIR:PATH=/usr/local/trilinos \
    -D CMAKE_C_COMPILER:STRING=/usr/bin/mpicc \
    -D CMAKE_CXX_COMPILER:STRING=/usr/bin/mpicxx \
    -D BOOST_ROOT=/usr/include/boost \
    -D USE_DAKOTA:BOOL=OFF \
    ..

WORKDIR /peridigm/build/
RUN make && make install
WORKDIR /
RUN mv /peridigm/scripts /usr/local/Peridigm/scripts
RUN rm -rf /peridigm

ENV LD_LIBRARY_PATH /usr/local/netcdf/lib

