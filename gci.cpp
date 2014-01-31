#include "gci.h"
#include "gciHamiltonian.h"
#include "gciDeterminant.h"
#include "gciWavefunction.h"
#include "gciString.h"
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
    xout << "Hamiltonian: " <<hh.printable(3)<<std::endl;
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
    xout << "Wavefunction after constructor:"<<w.printable(1)<<std::endl;
//    w.buildStrings();
//    xout << "Wavefunction after buildStrings:"<<w.printable(1)<<std::endl;

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

