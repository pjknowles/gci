# parent image
FROM ubuntu:20.04

# Install any needed packages
RUN apt-get update 
RUN apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git g++ gfortran doxygen graphviz bash rsync curl
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y mpich libhdf5-mpich-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libblas-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y liblapack-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libeigen3-dev

