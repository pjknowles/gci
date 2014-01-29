#include "gci.h"
#include "gciHamiltonian.h"
#include "gciDeterminant.h"
#include "gciString.h"
#include <iostream>
using namespace gci;
#ifndef MOLPRO
//int main(int argc, char *argv[])
int main()
{
    Hamiltonian hh("FCIDUMP");

    gci::Determinant d1(&hh);

    d1.create(3);
    d1.create(-1);

    String s1(&hh);
    int phase;
    unsigned int orbital;
    orbital=1;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<endl;
    orbital=2;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<endl;
    orbital=3;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<endl;
    s1.first(3);
//    orbital=5;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<endl;
    int n=0;
    for (bool i=true; i && n<20; n++,i=s1.next()) {
        xout << "Advance string=" <<s1.printable() <<endl;
    }
    xout <<"Total number of string="<<n<<std::endl;

  return 0;
}
#endif

