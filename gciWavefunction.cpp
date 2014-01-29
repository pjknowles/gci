#include "gciWavefunction.h"

Wavefunction::Wavefunction(unsigned int Nelec, int MS2)
{
    nelec = Nelec;
    ms2 = MS2;
}

void Wavefunction::setMS2(int MS2)
{
    ms2 = MS2;
}

void Wavefunction::setNelec(unsigned int Nelec)
{
    nelec = Nelec;
}

unsigned int Wavefunction::Nelec()
{
    return nelec;
}

int Wavefunction::MS2()
{
    return ms2;
}
