#include "gci.h"
#include "gciHamiltonian.h"
#include "gciDeterminant.h"
#include "gciWavefunction.h"
#include "gciStringSet.h"
#include "gciExcitationSet.h"
#include "gciOrbitalSpace.h"
#include "FCIdump.h"
#include <iostream>
using namespace gci;
#ifndef MOLPRO
//int main(int argc, char *argv[])
int main()
{
    try {
        FCIdump dump("FCIDUMP");
//        OrbitalSpace os("FCIDUMP");
//        xout <<"Orbital space:" << os << std::endl;
//       xout << "before Hamiltonian constructore"<<std::endl;
    Hamiltonian hh(&dump);
    xout << "Hamiltonian: " <<hh.str()<<std::endl;
    OrbitalSpace os = hh;
    xout << "Orbital space: " << os.str(1) <<std::endl;
//    exit(0);
    State ss(&dump);
    ss.orbitalSpace=&os;

//    Determinant d1(&ss);

//    d1.create(3);
//    d1.create(-1);


//    String s1(&ss);
//    int phase;
//    unsigned int orbital;
//    orbital=1;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
//    orbital=3;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
//    orbital=2;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
//    s1.first(3);
////    orbital=5;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
//    int n=0;
//    for (bool i=true; i && n<20; n++,i=s1.next()) {
//        xout << "Advance string=" <<s1.printable() <<std::endl;
//    }
//    xout <<"Total number of string="<<n<<std::endl;

//    xout << "Scan through constructing determinants:" << std::endl;
//    d1.first(); while(d1.next()) {
//        xout << " Determinant " <<d1.printable() <<std::endl;
//    }
//    xout <<"done scanning through determinants"<<std::endl;

    Wavefunction w(&dump);
    xout << "Wavefunction after constructor:"<<w.str(2)<<std::endl
    <<"...end of Wavefunction after constructor."<<std::endl<<std::endl;
    w.set((double)0.12345);
    xout << "Wavefunction after assign:"<<w.str(2)<<std::endl
    <<"...end of Wavefunction after assign."<<std::endl<<std::endl;
//    w.buildStrings();
//    xout << "Wavefunction after buildStrings:"<<w.str(1)<<std::endl;
    Wavefunction w2=w;
    xout << "Copied wavefunction:"<<w2.str(2)<<std::endl
    <<"...end of copied wavefunction."<<std::endl<<std::endl;
    xout << "Original wavefunction after copy:"<<w.str(2)<<std::endl
    <<"...end of original wavefunction."<<std::endl<<std::endl;
    w.set((double)1);
    xout << "Original wavefunction after original changed:"<<w.str(2)<<std::endl
    <<"...end of original wavefunction."<<std::endl<<std::endl;
    xout << "Copied wavefunction after original changed:"<<w2.str(2)<<std::endl
    <<"...end of copied wavefunction."<<std::endl<<std::endl;

    xout << "w.w=" << w*w << std::endl;
    xout << "w2.w2=" << w2*w2 << std::endl;
    xout << "Original wavefunction after w2.w2:"<<w.str(2)<<std::endl
    <<"...end of original wavefunction."<<std::endl<<std::endl;

    //Wavefunction w3;
    xout << "w:"<<w.str(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    w2 = w;
    xout << "after w2=w, w:"<<w.str(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    xout << "after w2=w, w2:"<<w2.str(2)<<std::endl <<"...end w2."<<std::endl<<std::endl;
    Wavefunction w3 = w2+(double)98*w-(w*((double)99));
    xout << "back from w3=..." <<std::endl;
    xout << "w:"<<w.str(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    xout << "w2:"<<w2.str(2)<<std::endl <<"...end w2."<<std::endl<<std::endl;
    xout << "w3:"<<w3.str(2)<<std::endl <<"...end w3."<<std::endl<<std::endl;


    w.diagonalHamiltonian(hh);
    xout << "Diagonal elements: " << w.str(2) << std::endl;
    size_t i = w.minloc();
    Determinant d = w.determinantAt(i);
    xout << "Lowest determinant " << d <<" with energy "<<w.at(i)<<std::endl;

    for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symb = syma ^ w.symmetry;
        std::vector<ExcitationSet> seta;
        seta = w.alphaStrings[syma].allExcitations(w.alphaStrings[syma],1,1);
        xout << "Excitations from alpha strings of symmetry " << syma+1 <<std::endl;
        for (std::vector<ExcitationSet>::iterator a=seta.begin(); a!=seta.end(); a++)
            xout <<"ExcitationSet: " <<a->str()<<std::endl;
        xout << "Alpha occupation numbers"<<std::endl;
        std::vector<double> on = w.alphaStrings[syma].occupationNumbers();
        for (size_t i=0; i < w.alphaStrings[syma].size(); i++) {
            for (size_t j=0; j < w.orbitalSpace->total(); j++)
                xout << " " << on[i+j*w.alphaStrings[syma].size()];
            xout << std::endl;
        }
    }


  return 0;
    }
    catch (const std::string& ex) {
        xout << "uncaught exception: " << ex << std::endl;
       throw "Error";
    }
    catch (char const* ex) {
        xout << "uncaught exception: " << ex << std::endl;
       throw "Error";
    }
//    catch (...) {
//        xout << "uncaught exception: "  << std::endl;
//       throw "Error";
//    }
}
#endif

