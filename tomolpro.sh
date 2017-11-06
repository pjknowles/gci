#!/bin/sh
if [ -z "$1" ]; then
    echo "Please enter path to Molpro tree"
    read tree
    set $tree
fi
echo "Synchronize to $1"

rsync src/*.{h,cpp} $1/src/gci
rsync submodules/memory/memory.h submodules/memory/bytestreamC* $1/src/global
rsync submodules/symmetry_matrix/*.{cpp,h} submodules/symmetry_matrix/SymmetryMatrixF.F90 $1/src/global
rsync submodules/Profiler/Profiler* $1/src/global
rsync submodules/FCIdump/FCIdump* $1/src/global
rsync submodules/IterativeSolver/LinearAlgebra.h $1/src/global
rsync submodules/IterativeSolver/IterativeSolver/I*.{cpp,h} $1/src/IterativeSolver
