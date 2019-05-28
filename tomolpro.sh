#!/bin/sh
if [ -z "$1" ]; then
    echo "Please enter path to Molpro tree"
    read tree
    set $tree
fi
echo "Synchronize to $1"

rm -f $1/src/global/{SMat,SMatMat,Operator,SMatfunction}.cpp
rm -f $1/src/gci/gciOperator.*
rsync src/gci.cpp $1/src/gci
rsync lib/*.{cpp,h} $1/src/gci
rsync dependencies/memory/lib/memory.F90 $1/src/global
rsync dependencies/memory/lib/memory.h $1/src/global
rsync dependencies/memory/lib/bytestreamC* $1/src/global
rsync dependencies/symmetry_matrix/lib/{SMat,Operator,SMatMat}.h dependencies/symmetry_matrix/lib/SymmetryMatrixF.F90 $1/src/global
rsync dependencies/symmetry_matrix/lib/{SMat,Operator,SMatMat}-{double.cpp,implementation.h} $1/src/blas
rsync dependencies/Profiler/lib/Profiler* $1/src/global
rsync dependencies/FCIdump/lib/FCIdump* $1/src/global
rsync -av dependencies/IterativeSolver/{lib,example,CMakeLists.txt,Doxyfile.in,test} $1/src/IterativeSolver
rsync gci.registry $1/lib
