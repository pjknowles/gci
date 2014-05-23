cmake_minimum_required(VERSION 2.8)

file (GLOB SRCS *.cpp *.h)
add_executable(gci.exe ${SRCS})

configure_file(FCIDUMP FCIDUMP COPYONLY)
