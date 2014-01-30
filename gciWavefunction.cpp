#include "gciWavefunction.h"

void Wavefunction::buildStrings()
{
    String string;
    string.first((nelec+ms2)/2);
    String::buildStrings(string,&alphaStrings);
    string.first((nelec-ms2)/2);
    String::buildStrings(string,&betaStrings);
}
