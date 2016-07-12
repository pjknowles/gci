#ifndef GCIMOLPRO_H
#define GCIMOLPRO_H
#include <cstdlib>
#ifdef MOLPRO
#include <mpp/CxMpp.h>
#include <cic/ItfMpp.h>
#include <cic/ItfCommon.h>
#include <cic/ItfFortranInt.h>
using namespace itf;
#else
typedef unsigned int uint;
// matrix multiplication routines. All of them call mxma/mxmb as appropriate
//
// C = A^T B matrix multiplication driver routine proxy:
// C_rs = A_tr B_ts.  dim(r)=nRows, dim(s)=nCols, dim(t)=nLink
// Matrices are col-major (i.e., t is the fast dimension) with
// col-stride==nrows. A^T B is supposed to be the fastest
// MxM-driver on most machine architectures.
void MxmDrvTN(double *Out, const double *A,
    double *B, uint nRows, uint nLink, uint nStrideLink, uint nCols, bool AddToDest = false);
// as above, but for C = A B
void MxmDrvNN( double *Out, double *A,
    double *B, uint nRows, uint nLink, uint nCols, bool AddToDest = false);

// as MOLPRO's mxma/mxmb, except for 'Out' occuring before A/B.
void MxmDrvGen( double *Out, uint nRowStOut, uint nColStOut,
    double *A, uint nRowStA, uint nColStA,
    double *B, uint nRowStB, uint nColStB,
    uint nRows, uint nLink, uint nCols, bool AddToDest = false );
// adds support for a prefactor factor, but resorts to BLAS for
// actual computation. Therefore for each matrix one of the strides
// must be 1! Also, this will likely be very slow for small matrices.
void MxmDrvGenF( double *Out, uint nRowStOut, uint nColStOut,
    double *A, uint nRowStA, uint nColStA,
    double *B, uint nRowStB, uint nColStB,
    uint nRows, uint nLink, uint nCols, bool AddToDest, double Factor );

// also a generic driver for mxva/b
void MxvDrvGen(double *A, uint nRowStA, uint nColStA,
    double *V, uint nStrideV, double *R, uint nStrideR,
    uint nRowsR, uint nLink, bool AddToDest = false);

// diagonalize a symmetric matrix in-place, storing nDim eigenvalues at pEigValues.
void Diagonalize( double *pMatrix, double *pEigValues, uint nDim, uint nColStride );


#endif
#endif // GCIMOLPRO_H
