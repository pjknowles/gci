cmake_minimum_required(VERSION 2.8)
add_executable(gci.exe
 gci.cpp 
    gciHamiltonian.cpp 
    gciState.cpp 
    gciDeterminant.cpp 
    gciString.cpp 
    gciWavefunction.cpp 
    FCIdump.cpp 
    gciStringSet.cpp 
    gciExcitationSet.cpp 
    gciOperator.cpp 
    gciPrintable.cpp 
    gciSymmetrySpace.cpp 
    gciOrbitalSpace.cpp 
    gciTransitionDensity.cpp 
    gciFile.cpp 
    gciMolpro.cpp
 gci.h 
    gciHamiltonian.h 
    gciState.h 
    gciDeterminant.h 
    gciString.h 
    gciWavefunction.h 
    FCIdump.h 
    gciStringSet.h 
    gciExcitationSet.h 
    gciOperator.h 
    gciPrintable.h 
    gciSymmetrySpace.h 
    gciOrbitalSpace.h 
    gciTransitionDensity.h 
    gciFile.h 
    gciMolpro.h
)

configure_file(FCIDUMP FCIDUMP COPYONLY)
