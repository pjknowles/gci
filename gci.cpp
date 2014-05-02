#include "gci.h"
#include "gciHamiltonian.h"
#include "gciDeterminant.h"
#include "gciWavefunction.h"
#include "gciStringSet.h"
#include "gciExcitationSet.h"
#include "FCIdump.h"
#include <iostream>
using namespace gci;
#ifndef MOLPRO
//int main(int argc, char *argv[])
int main()
{
    try {
        FCIdump dump("FCIDUMP");
    Hamiltonian hh(&dump);
    xout << "Hamiltonian: " <<hh.toString(3)<<std::endl;
    State ss(&dump);
    ss.hamiltonian=&hh;

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
    xout << "Wavefunction after constructor:"<<w.toString(2)<<std::endl
    <<"...end of Wavefunction after constructor."<<std::endl<<std::endl;
    w.set((double)0.12345);
    xout << "Wavefunction after assign:"<<w.toString(2)<<std::endl
    <<"...end of Wavefunction after assign."<<std::endl<<std::endl;
//    w.buildStrings();
//    xout << "Wavefunction after buildStrings:"<<w.toString(1)<<std::endl;
    Wavefunction w2=w;
    xout << "Copied wavefunction:"<<w2.toString(2)<<std::endl
    <<"...end of copied wavefunction."<<std::endl<<std::endl;
    xout << "Original wavefunction after copy:"<<w.toString(2)<<std::endl
    <<"...end of original wavefunction."<<std::endl<<std::endl;
    w.set((double)1);
    xout << "Original wavefunction after original changed:"<<w.toString(2)<<std::endl
    <<"...end of original wavefunction."<<std::endl<<std::endl;
    xout << "Copied wavefunction after original changed:"<<w2.toString(2)<<std::endl
    <<"...end of copied wavefunction."<<std::endl<<std::endl;

    xout << "w.w=" << w*w << std::endl;
    xout << "w2.w2=" << w2*w2 << std::endl;
    xout << "Original wavefunction after w2.w2:"<<w.toString(2)<<std::endl
    <<"...end of original wavefunction."<<std::endl<<std::endl;

    //Wavefunction w3;
    xout << "w:"<<w.toString(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    w2 = w;
    xout << "after w2=w, w:"<<w.toString(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    xout << "after w2=w, w2:"<<w2.toString(2)<<std::endl <<"...end w2."<<std::endl<<std::endl;
    Wavefunction w3 = w2+(double)98*w-(w*((double)99));
    xout << "back from w3=..." <<std::endl;
    xout << "w:"<<w.toString(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    xout << "w2:"<<w2.toString(2)<<std::endl <<"...end w2."<<std::endl<<std::endl;
    xout << "w3:"<<w3.toString(2)<<std::endl <<"...end w3."<<std::endl<<std::endl;


    for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symb = syma ^ w.symmetry;
        std::vector<ExcitationSet> seta;
        seta = w.alphaStrings[syma].allExcitations(w.alphaStrings[syma],1,1);
        xout << "Excitations from alpha strings of symmetry " << syma+1 <<std::endl;
        for (std::vector<ExcitationSet>::iterator a=seta.begin(); a!=seta.end(); a++)
            xout <<"ExcitationSet: " <<a->toString()<<std::endl;
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

