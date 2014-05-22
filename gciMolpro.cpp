#include "gciMolpro.h"

void MxmDrvGen( double *Out, uint nRowStOut, uint nColStOut,
    double *A, uint nRowStA, uint nColStA,
    double *B, uint nRowStB, uint nColStB,
    uint nRows, uint nLink, uint nCols, bool AddToDest)
{
  if (! AddToDest)
    for (uint s=0; s<nCols; s++)
      for (uint r=0; r<nRows; r++)
        Out[r*nRowStOut+s*nColStOut]=(double)0;
  for (uint s=0; s<nCols; s++)
    for (uint r=0; r<nRows; r++)
      for (uint t=0; t<nLink; t++)
        Out[r*nRowStOut+s*nColStOut] +=
            A[r*nRowStA+t*nColStA] * B[t*nRowStB+s*nColStB];
}

