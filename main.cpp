#include "TreeCIHamiltonian.h"
#include "TreeCIDeterminant.h"
using namespace TreeCI;
int main(int argc, char *argv[])
{

    Hamiltonian hh("FCIDUMP");

    TreeCI::Determinant d1(&hh);
    d1.create(3);
    d1.create(-1);

  return 0;
}
