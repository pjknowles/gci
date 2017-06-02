#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H
#include "Operator.h"
#include "FCIdump.h"


namespace gci {

  SymmetryMatrix::Operator makeOperator(const FCIdump& dump);
}

#endif // GCIOPERATOR_H
