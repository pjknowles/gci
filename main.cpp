#include "TreeCIHamiltonian.h"
#include "TreeCIDeterminant.h"

int main(int argc, char *argv[])
{

    TreeCIHamiltonian hh("FCIDUMP");

    TreeCIDeterminant d1(&hh);
    d1.create(3);
    d1.create(-1);

  return 0;
}
