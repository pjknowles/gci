FROM ubuntu:jammy-20220531
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake git g++ doxygen graphviz bash rsync curl
RUN DEBIAN_FRONTEND=noninteractive apt install -y mpich libhdf5-mpich-dev
RUN DEBIAN_FRONTEND=noninteractive apt install -y libopenblas-dev
RUN DEBIAN_FRONTEND=noninteractive apt install -y liblapack-dev liblapacke-dev
RUN DEBIAN_FRONTEND=noninteractive apt install -y libeigen3-dev ccache
RUN DEBIAN_FRONTEND=noninteractive apt install -y wget
RUN export VERSION=5.8.1 && wget https://github.com/GlobalArrays/ga/releases/download/v${VERSION}/ga-${VERSION}.tar.gz && tar xzf ga-${VERSION}.tar.gz && cd ga-${VERSION} && ./configure --disable-f77 --with-mpi3 --without-blas --without-lapack && make install && cd .. && rm -rf ga-${VERSION}*
