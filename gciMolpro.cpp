#include "gciMolpro.h"
#ifndef MOLPRO

void MxmDrvNN(double *Out, double *A, double *B, uint nRows, uint nLink, uint nCols, bool AddToDest)
{
  if (! AddToDest)
    for (uint s=0; s<nCols; s++)
      for (uint r=0; r<nRows; r++)
        Out[r+s*nRows]=(double)0;
  for (uint s=0; s<nCols; s++)
    for (uint r=0; r<nRows; r++)
      for (uint t=0; t<nLink; t++)
        Out[r+s*nRows] +=
            A[r+t*nRows] * B[t+s*nLink];
}

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

#include <cmath>
#include <assert.h>
void Diagonalize( double *x, double *d, unsigned int m, unsigned int nColStride ) {
  assert (nColStride == m);
  unsigned int n=m;
//       subroutine diag2(m,n,d,x)
// c
// c      computation of all eigenvalues and eigenvectors of a real
// c      symmetric matrix by the method of qr transformations.
// c      if the euclidean norm of the rows varies   s t r o n g l y
// c      most accurate results may be obtained by permuting rows and
// c      columns to give an arrangement with increasing norms of rows.
// c
// c      two machine constants must be adjusted appropriately,
// c      eps = minimum of all x such that 1+x is greater than 1 on the
// c            computer,
// c      tol = inf / eps  with inf = minimum of all positive x represen-
// c            table within the computer.
// c      a dimension statement e(160) may also be changed appropriately.
// c
// c      input
// c
// c      (m)   not larger than 160,  corresponding value of the actual
// c            dimension statement a(m,m), d(m), x(m,m),
// c      (n)   not larger than (m), order of the matrix,
// c      (a)   the matrix to be diagonalized, its lower triangle has to
// c            be given as  ((a(i,j), j=1,i), i=1,n),
// c.....
// c.....the matrix #a# has been removed from the procedure
// c.....the matrix #x# has to be put up by a predecessor routine c.....
// c
// c      output
// c
// c      (d)   components d(1), ..., d(n) hold the computed eigenvalues
// c            in ascending sequence. the remaining components of (d) are
// c            unchanged,
// c      (x)   the computed eigenvector corresponding to the j-th eigen-
// c            value is stored as column (x(i,j), i=1,n). the eigenvectors
// c            are normalized and orthogonal to working accuracy. the
// c            remaining entries of (x) are unchanged.
// c
// c      array (a) is unaltered. however, the actual parameters
// c      corresponding to (a) and (x)  may be identical, ''overwriting''
// c      the eigenvectors on (a).
// c
// c      leibniz-rechenzentrum, munich 1965
// c
  const int maxdim=500;
  const double eps=2.5e-16,dinf=2.3e-308,tol=dinf/eps;
  assert (m <= maxdim);
  double e[maxdim];
  if(n == 1)  {
    //     special treatment of case n = 1
    d[0]=x[0];
    x[0]=(double) 1;
    return;
  }
  for (int i=0; i<n; i++) {
    d[i]=(double) 0;
    e[i]=(double) 0;
  }

  //     householder's reduction

  for (int i=n-1;i>0;i--) {
    int l=i-1;
    double h=(double)0;
    double g=x[i + i*m];
    if (l>0) {
      for (int k=0; k<l; k++)
	h+=x[i+k*m]*x[i+k*m];
      double s=h+g*g;
      if(s < tol)
        h=(double)0;
      else if(h > 0){
	l++;
        double f=g;
        g=std::sqrt(s);
        if(f >0) g=-g;
        h=s-f*g;
        x[i+(i-1)*m]=f-g;
        f=(double)0;

	for (int j=0; j<l; j++) {
	  x[j+m*i]=x[i+m*j]/h;
	  s=(double)0;
	  for (int k=0; k<=j; k++)
	    s+=x[j+m*k]*x[i+m*k];
	  //        j1=j+1
	  //if(j1.gt.l) go to 100
	  for (int k=j+1; k<l; k++)
	    s+=x[k+m*j]*x[i+m*k];
	  e[j]=s/h;
	  f+=s*x[j+m*i];
	}

	f=f/(2*h);

	for (int j=0; j<l; j++)
	  e[j]-=f*x[i+m*j];

	for (int j=0; j<l; j++) {
	  f=x[i+m*j];
	  s=e[j];
	  for (int k=0; k<=j; k++)
	    x[j+m*k]-=(f*e[k]+x[i+m*k]*s);
	}

      }
    }
    d[i]=h;
    e[i-1]=g;
  }

  //     accumulation of transformation matrices

  d[0]=x[0];
  x[0]=(double)1;
  for (int i=1; i<n; i++) {
    if (d[i] > (double)0) {
      for (int j=0; j<i; j++) {
	double s=(double)0;
	for (int k=0; k<i; k++)
	  s+=x[i+m*k]*x[k+m*j];
	for (int k=0; k<i; k++)
	  x[k+m*j]-=s*x[k+m*i];
      }
    }
    d[i]=x[i+m*i];
    x[i+m*i]=(double)1;
      for (int j=0; j<i; j++) {
	x[i+m*j]=(double)0;
	x[j+m*i]=(double)0;
      }
  }

  //     diagonalization of the tridiagonal matrix

  double b=(double)0;
  double f=(double)0;
  e[n-1]=(double)0;

  for (int l=0; l<n; l++) {
    double h=eps*(std::abs(d[l])+std::abs(e[l]));
    if (h > b) b=h;

    //     test for splitting

    int j;
    for (j=l; j<n; j++)
      if (std::abs(e[j]) <= b) break;

    //     test for convergence

    if(j != l) {
      while (std::abs(e[l]) > b) {

	//     shift from upper 2*2 minor

	double p=(d[l+1]-d[l])*(double)0.5;
	double r=std::sqrt(p*p+e[l]*e[l]);
	if (p<0)
	  p+=r;
	else
	  p-=r;
	h=d[l]+p;
	for (int i=l;i<n;i++)
	  d[i]-=h;
	f+=h;

	//     qr transformation

	p=d[j];
	double c=(double)1;
	double s=(double)0;

	for (int i=j-1; i>=l; i--) {
	  double g=c*e[i];
	  h=c*p;

	  //     protection against underflow of exponents

	  if (std::abs(p)>=std::abs(e[i])) {
	    c=e[i]/p;
	    r=std::sqrt(c*c+(double)1);
	    e[i+1]=s*p*r;
	    s=c/r;
	    c=(double)1/r;
	} else {
	    c=p/e[i];
	    r=std::sqrt(c*c+(double)1);
	    e[i+1]=s*e[i]*r;
	    s=(double)1/r;
	    c=c/r;
	}
	  p=c*d[i]-s*g;
	  d[i+1]=h+s*(c*g+s*d[i]);
	  for (int k=0; k<n; k++) {
	    h=x[k+m*(i+1)];
	    x[k+m*(i+1)]=x[k+m*i]*s+h*c;
	    x[k+m*i]=x[k+m*i]*c-h*s;
	  }
	}

	e[l]=s*p;
	d[l]=c*p;
      }
    }

    //     convergence

    d[l]=d[l]+f;
  }

  //     ordering of eigenvalues

  for (int i=0; i<n-1; i++) {
    int k=i;
    double p=d[i];
    for (int j=i+1; j<n; j++) {
      if(d[j] <p) {
	k=j;
	p=d[j];
      }
    }
    if (k != i) {
      d[k] =d[i];
      d[i]=p;
      for (int j=0; j<n; j++) {
	p=x[j+m*i];
	x[j+m*i]=x[j+m*k];
	x[j+m*k]=p;
      }
    }
  }

  //     fixing of sign

  for (int i=0; i<n-1; i++) {
    double pm=(double)0;
    int k;
    for (int j=0; j<n-1; j++) {
      if(pm <= std::abs(x[j+m*i]))  {
	pm =std::abs(x[j+m*i]);
	k=j;
      }
    }
    if(x[k+m*i] < (double)0) {
    for (int j=0; j<n-1; j++)
	   x[j+m*i]=-x[j+m*i];
    }
  }
}

#endif
