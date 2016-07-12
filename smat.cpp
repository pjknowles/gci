#include "smat.h"

#ifdef MOLPRO
#include "cic/ItfFortranInt.h"
#include "cic/ItfBlasInt.h"
#define dsyev_y FD(dsyev_x)
#define dgemm_y FD(dgemm_x)
#define dger_y FD(dger_x)
#define dgemv_y FD(dgemv_x)
#define dgeev_y FD(dgeev_x)
#else
#define xout std::cout
#define FD(x) x ## _
typedef int64_t FORTINT;
typedef double FORTDBL;
typedef FORTINT const & FINTARG;
typedef FORTDBL const & FDBLARG;
#include <mkl_lapack.h>
#include <mkl_blas.h>

extern "C" {
void dgemm_y(char const &TransA, char const &TransB, FINTARG M, FINTARG N, FINTARG K,
             FORTDBL const &Alpha, FORTDBL const *A, FINTARG lda, FORTDBL const *B, FINTARG ldb,
             FORTDBL const &Beta, FORTDBL *C, FINTARG ldc)
{
  MKL_INT M_=M;
  MKL_INT N_=N;
  MKL_INT K_=K;
  MKL_INT lda_=lda;
  MKL_INT ldb_=ldb;
  MKL_INT ldc_=ldc;
  dgemm(&TransA, &TransB, &M_, &N_, &K_,
        &Alpha, A, &lda_, B, &ldb_,
        &Beta, C, &ldc_);
}

void dgemv_y(char const &Trans, FINTARG M, FINTARG N,
             FORTDBL const &Alpha, FORTDBL const *A, FINTARG lda, FORTDBL const *X, FINTARG incx,
             FORTDBL const &Beta, FORTDBL *Y, FINTARG incy);

void dger_y(FINTARG M, FINTARG N, FORTDBL const &Alpha, FORTDBL const *X, FINTARG incx, FORTDBL const *Y, FINTARG incy, FORTDBL *A, FINTARG lda);

void dsyev_y(char const &jobz, char const &uplo, FINTARG N, FORTDBL *A, FINTARG lda, FORTDBL* w, FORTDBL* work, FINTARG lwork, FORTINT& info)
{
  MKL_INT N_=N;
  MKL_INT lda_=lda;
  MKL_INT lwork_=lwork;
  MKL_INT info_=info;
  dsyev(&jobz, &uplo, &N_, A, &lda_, w, work, &lwork_, &info_);
  info=info_;
}

void dgeev_y(char const &jobvl, char const &jobvr, FINTARG N, FORTDBL *A, FINTARG lda, FORTDBL* wr, FORTDBL* wi,
             FORTDBL* vl, FINTARG ldvl,
             FORTDBL* vr, FINTARG ldvr,
             FORTDBL* work, FINTARG lwork, FORTINT& info)
{
  MKL_INT N_=N;
  MKL_INT lda_=lda;
  MKL_INT ldvl_=ldvl;
  MKL_INT ldvr_=ldvr;
  MKL_INT lwork_=lwork;
  MKL_INT info_=info;
  dgeev(&jobvl, &jobvr, &N_, A, &lda_,  wr, wi,
        vl, &ldvl_,
        vr, &ldvr_,
        work, &lwork_, &info_);
}
}
#include <iostream>
#endif

 smat::smat(std::vector<std::vector<size_t> > dimensions, int parity, int symmetry, unsigned int rank, int create, double *buffer, std::string description)
{
     initialise(dimensions,symmetry,parity,rank,create,buffer,description);
 }

#ifdef MOLPRO
 smat::smat(std::string space, int parity, int symmetry, int create, double *buffer, std::string description)
 {
     std::vector<std::vector<size_t> > spaces;
     std::vector<size_t> nt1(8);
     smat_get_orbital_space(space[0],&(nt1[0]));
     spaces.push_back(nt1);
     smat_get_orbital_space(space[(space.size()>1 ? 1:0)],&(nt1[0]));
     spaces.push_back(nt1);
     initialise(spaces,symmetry,parity,space.size(),create,buffer,description);
 }
#endif

smat::smat(smat const & source, int parity, int symmetry, unsigned int rank, int create, double* buffer, std::string description)
{
  int symmetry_=(symmetry==9?source.symmetry:symmetry);
  int parity_=(parity==9?source.parity:parity);
  int rank_=(rank==0?source.rank:rank);
  initialise(source.dimensions,symmetry_,parity_,rank_,create,buffer,description);
}

smat::smat()
{
  std::vector<std::vector<size_t> > spaces;
  std::vector<size_t> nt1(8,0);
  spaces.push_back(nt1);
  spaces.push_back(nt1);
  initialise(spaces);
}

void smat::initialise(std::vector<std::vector<size_t> > dimensions, int symmetry, int parity, unsigned int rank, int create, double *buffer, std::string description)
{
  this->dimensions.resize(2);
  for (size_t i=0; i < dimensions.size(); i++) (this->dimensions)[i] = dimensions[i];
  this->symmetry = symmetry;
  this->parity = parity;
  this->rank = rank;
  this->buffer = buffer;
  if (create > 0 && buffer != NULL) throw "Illegal specification of create and buffer in smat::smat";
  if (create > 0 || (create < 0 && buffer == NULL)) ensure_buffer();
}

#include <stdlib.h>
#include "memory.h"
void smat::ensure_buffer()
{
  size_t desired_size=size();
  if (buffer != NULL) {
    // no smart resize
    if (buffer_size >= desired_size) return;
    memory_release(buffer);
  }
  buffer = (double *) memory_allocate(desired_size*sizeof(double));
  buffer_size = desired_size;
  block.resize(0);
  for (unsigned int block_symmetry=0; block_symmetry < 8; block_symmetry++) {
    block.push_back( &(buffer[block_offset(block_symmetry)]) );
  }
}


smat::~smat()
{
  if (buffer != NULL) memory_release(buffer);
}

size_t smat::size() const
{
  size_t result=0;
  if (rank == 1)
    if (symmetry >=0)
      result = dimensions[0][symmetry];
    else
      for (int i=0; i<8; i++)
	result += dimensions[0][i];
  else
    result = block_offset(8);
  return result;
}

size_t smat::block_offset(unsigned int block_symmetry) const
{
  size_t result=0;
  if (rank == 2) {
    unsigned int bs=block_symmetry, bs2=symmetry^block_symmetry;
    if (bs2 > bs) {bs2=block_symmetry; bs=symmetry^bs2;}
    for (unsigned int ks=0; ks < bs; ks++) {
      unsigned int ls=ks^symmetry;
      if (parity != 0 && symmetry == 0)
	result += dimensions[0][ks]*(dimensions[0][ks]+1)/2;
      else if (ks > ls || parity == 0)
	result += dimensions[0][ks] * dimensions[1][ls];
    }
  } else if (rank == 1) {
    if (symmetry < 0)
      for (unsigned int ks=0; ks < block_symmetry; ks++)
	result += dimensions[0][ks];
    else if (block_symmetry > (unsigned int)symmetry)
      result = dimensions[0][symmetry]; // usual use: block_symmetry=9 for total size
  }
  return result;
}

std::vector<size_t> smat::block_dimensions(unsigned int block_symmetry) const
{
  std::vector<size_t> result;
  if (rank == 2) {
    size_t transpose = ((block_symmetry^(unsigned int)symmetry) > block_symmetry && parity != 0) ? 1 : 0;
    if (transpose) {
      result.push_back(dimensions[0][block_symmetry^symmetry]);
      result.push_back(dimensions[0][block_symmetry]);
      result.push_back(1);
    } else {
      result.push_back(dimensions[0][block_symmetry]);
      result.push_back(dimensions[1][block_symmetry^symmetry]);
      result.push_back(0);
    }
  } else if (rank == 1) {
    result.resize(1,0);
    if (symmetry < 0)
      for (unsigned int ks=0; ks < block_symmetry; ks++)
	result[0] += dimensions[0][ks];
    else if (block_symmetry == (unsigned int)symmetry)
      result[0] = dimensions[0][symmetry];
    result.push_back(1);
    result.push_back(0);
  }
  return result;
}

double smat::trace() const
{
  if ((rank != 2) || (symmetry != 0) || (parity < 0) ) return (double)0;
  double result=(double)0;
  for (unsigned int ks=0; ks<8; ks++) {
    double* buff=&buffer[block_offset(ks)];
    for (size_t i=0; i<dimensions[0][ks]; i++) {
      result += *buff;
      buff += (parity>0) ? i+2 : dimensions[0][ks]+1;
    }
  }
  return result;
}

void smat::zero()
{
  if (buffer == NULL) return;
  for (size_t i=0; i < size(); i++) buffer[i]=(double)0;
}

#include <sstream>
std::string smat::str(std::string title, int number) const
{
  std::stringstream s;
  s << "Matrix "<<title;
  if(number!=999999) s <<"("<<number<<") " ;
  s <<"; symmetry="<<symmetry+1<<" parity="<<parity<<std::endl;
  // s << "buffer:"; for (size_t k=0; k<size(); k++) s<<" "<<buffer[k]; s<<std::endl;
  for (unsigned int ks=0; ks<8; ks++) {
    size_t nr=dimensions[0][ks], nc=(rank > 1 ? dimensions[1][ks^symmetry] : 1);
    if (nr < 1 || nc < 1) continue;
    if (rank == 1 && symmetry >=0 && (unsigned int) symmetry != ks) continue;
    if (rank == 2 && parity != 0 && ks > (ks^symmetry)) continue;
    if (rank == 1)
      s << "Block ("<<ks+1<<"), dimensions ("<<nr<<")"<<std::endl;
    else
      s << "Block ("<<ks+1<<","<<(ks^symmetry)+1<<"), dimensions ("<<nr<<","<<nc<<")"<<std::endl;
    if (buffer==NULL) continue;
    if (rank == 1) {
      for (size_t k=0; k<nr; k++) s <<" "<<block[ks][k];
      s << std::endl;
    }
    else {
      for (size_t k=0; k<nr; k++) {
	if (parity != 0 && symmetry==0)
	  for (size_t l=0; l<=k; l++)
	    s <<" "<<block[ks][(k)*(k+1)/2+l];
	else
	  for (size_t l=0; l<nc; l++)
	    s <<" "<<block[ks][k+l*nr];
	s << std::endl;
      }
    }
  }

  return s.str();
}

void smat::copy(smat const & source, int parity, std::string description)
{
  if (parity != 999999) this->parity=parity;
  if (buffer!=NULL) { // buffer exists: honour existing dimensions
    for (size_t k=0; k<2; k++)
      if (dimdiff(dimensions[k],source.dimensions[k])) {
	xout << (k ? "Rows:" : "Columns:"); for (size_t l=0; l<8; l++) xout <<"  source="<<source.dimensions[k][l]<<", this="<<dimensions[k][l]; xout <<std::endl;
	throw "Incorrect receiving dimensions, smat::copy";
      }
    if (buffer_size != size()) {
      xout << "buffer_size="<<buffer_size<<", size()="<<size()<<std::endl;
      throw "Wrong receiving buffer size, smat::copy";
    }
  } else { // buffer does not exist: force dimensions from source
    for (size_t k=0; k<2; k++)
      dimensions[k] = source.dimensions[k];
    this->ensure_buffer();
  }
  if (this->description=="") this->description=source.description;
  if (description != "") this->description=description;
  if (rank==2 && (source.parity || this->parity) && dimdiff(dimensions[0],dimensions[1]))
    throw "smat::copy: matrix must be square";

  if (rank == 1 || this->parity == source.parity)
    for (size_t i=0; i<size(); i++)
      buffer[i] = source.buffer[i];
  else if (this->parity == 0) {  // copy triangle to square
    for (unsigned int ks=0; ks<8; ks++) {
      double* to=this->block[ks];
      unsigned int ls=ks^symmetry;
      size_t nk=dimensions[0][ks];
      size_t nl=dimensions[1][ls];
      if (ks == ls) {
	double* from=source.block[ks];
	if (source.parity > 0)
	  for (size_t k=0; k<nk; k++)
	    for (size_t l=0; l<=k; l++) {
	       to[k*nl+l] = to[l*nk+k] = from[k*(k+1)/2+l];
	    }
	else
	  for (size_t k=0; k<nk; k++)
	    for (size_t l=0; l<=k; l++)
	      to[k*nl+l] = -(to[l*nk+k] = from[k*(k+1)/2+l]);
      } else if (ks > ls) {
	double* from=source.block[ks];
	for (size_t kl=0; kl<nk*nl; kl++)
	  to[kl] = from[kl];
      } else {
	double* from=source.block[ls];
	if (this->parity > 0)
	  for (size_t k=0; k<nk; k++)
	    for (size_t l=0; l<=k; l++)
	      to[l*nk+k] = from[k*nl+l];
	else
	  for (size_t k=0; k<nk; k++)
	    for (size_t l=0; l<=k; l++)
	      to[l*nk+k] = -from[k*nl+l];
      }
    }
  } else { // copy square to triangle
    for (unsigned int ks=0; ks<8; ks++) {
      double* from=source.block[ks];
      unsigned int ls=ks^symmetry;
      size_t nk=dimensions[0][ks];
      size_t nl=dimensions[1][ls];
      if (ks == ls) {
	double* to=this->block[ks];
	for (size_t k=0; k<nk; k++)
	  for (size_t l=0; l<=k; l++)
	    to[k*(k+1)/2+l] = from[k*nk+l];
      } else if (ks > ls) {
	double* to=this->block[ks];
	for (size_t kl=0; kl<nk*nl; kl++)
	  to[kl] = from[kl];
      }
    }
  }
}

void smat::scal(double a)
{
  for (size_t i=0; i<size(); i++)
    buffer[i] *= a;
}

void smat::axpy(double a, smat& x)
{
  for (size_t i=0; i<size(); i++)
    buffer[i] += a*x.buffer[i];
}

void smat::ger(double alpha, smat &x, smat &y)
{
  throw "not implemented";
}

void smat::gemm(smat const & a, smat const & b, char transa, char transb, double alpha, double beta)
{
  int transposea=(transa!='N'&&transa!='n')?1:0;
  int transposeb=(transb!='N'&&transb!='n')?1:0;
  if (a.rank !=2) throw "a not a matrix in smat::gemm";
  if (b.rank !=2) throw "b not a matrix in smat::gemm";
  if (b.symmetry < 0 || b.symmetry > 7) throw "b.symmetry problem";
  unsigned int bs=b.symmetry;
  for (unsigned int ks=0; ks<8; ks++)
    if (a.dimensions[1-transposea][ks] != b.dimensions[transposeb][ks]) throw "Mismatch in a,b dimensions in smat::gemm";
  symmetry=a.symmetry^b.symmetry;
  for (unsigned int ks=0; ks<8; ks++)
    dimensions[0][ks] = a.dimensions[transposea][ks];
  for (unsigned int ks=0; ks<8; ks++)
    dimensions[1][ks] = b.dimensions[1-transposea][ks];
  rank = 2;
  parity = 0;
  ensure_buffer();
  if (size()==0) return;
  if (a.parity!=0 || b.parity!=0) throw "Symmetric matrix support not yet complete in smat::gemm";
  for (unsigned int ks=0; ks<8; ks++) {
    unsigned int ls=ks^a.symmetry;
    double* cblock=block[ks];
    size_t cr=dimensions[0][ks];
    size_t cc=dimensions[1][ks^symmetry];
    if (cr*cc == 0) continue;
    double* ablock=a.block[transposea ? ls : ks];
    size_t lda=a.dimensions[0][transposea?ls:ks];
    size_t k=a.dimensions[1-transposea][ls];
    unsigned int lst=ls; if (transposeb) lst=ls^bs;
    double* bblock=b.block[lst];
    size_t ldb=b.dimensions[0][lst];
    dgemm_y(transa,transb,cr,cc,k,alpha,ablock,lda,bblock,ldb,beta,cblock,cr);
  }
}

void smat::gemv(smat const &a, smat const &x, char transa, double alpha, double beta)
{
  throw "not implemented";
}


#include <cmath>
bool smat::dimdiff(std::vector<size_t> const v1, std::vector<size_t> const v2) const
{
  int diff=0; for (size_t ks=0; ks<8; ks++) diff+=std::fabs(v1[ks]-v2[ks]);
  return diff!=0;
}

#include <limits>
int smat::ev(smat & val, smat* vec, smat* vali, smat* vecl, const std::string algorithm, const std::string sort) const
{
  FORTINT info=0;
  size_t base = memory_save();
  if (this->rank!=2 || dimdiff(dimensions[0],dimensions[1])) throw "Eigenvalues/vectors only for square matrix, smat::ev";
  if (this->symmetry!=0) throw "Eigenvalues/vectors only for matrix of symmetry 1, smat::ev";
  if (val.rank != 1 || val.buffer==NULL || val.symmetry !=-1 || dimdiff(dimensions[0],val.dimensions[0])) throw "Invalid val, smat::ev";
  val.description="Eigenvalues";
  if (vec != NULL && (vec->rank != 2 || vec->buffer==NULL || vec->symmetry !=1 || dimdiff(dimensions[0],vec->dimensions[0])) ) throw "Invalid vec, smat::ev";
  if (vali != NULL && (vali->rank != 1 || vali->buffer==NULL || vali->symmetry !=-1 || dimdiff(dimensions[0],vali->dimensions[0])) ) throw "Invalid vali, smat::ev";
  if (vecl != NULL && (vecl->rank != 2 || vecl->buffer==NULL || vecl->symmetry !=1 || dimdiff(dimensions[0],vecl->dimensions[0])) ) throw "Invalid vecl, smat::ev";
  smat* vvpt=vec; if (vec == NULL) vvpt = new smat(*this,0);
  vvpt->copy(*this,0);
  vvpt->description="Eigenvectors";
  char jobz = (vec==NULL ? 'N' : 'V');
  char jobvl = (vecl==NULL ? 'N' : 'V');
  for (unsigned int k=0; k<8; k++) {
    size_t n=this->dimensions[0][k];
    if (n < 1) continue;
    size_t bases = memory_save();
    double* vblock=vvpt->block[k];
    double* valblock=val.block[k];
    if (this->parity>0) {
      double lwork;
      dsyev_y(jobz,'L',n,vblock,n,valblock,&lwork,-1,info);
      double* work = memory_allocate(sizeof(double)*(size_t)(lwork+1));
      dsyev_y(jobz,'L',n,vblock,n,valblock,work,(size_t)lwork,info);
      if (info!=0) { memory_release_saved(base); return info; }
    } else if (this->parity == 0) {
      double* mpt = memory_allocate(sizeof(double)*n*n);
      double* mpt2 = this->block[k];
      for (size_t l=0; l<n*n; l++) mpt[l]=mpt2[l];
      double* vpt = memory_allocate(sizeof(double)*n);
      mpt2 = NULL; if (vecl!=NULL) mpt2=vecl->block[k];
      double* mpt3 = NULL; if (vec!=NULL) mpt3=vec->block[k];
      double lwork;
      dgeev_y(jobvl,jobz,n,mpt,n,valblock,vpt,mpt2,n,mpt3,n,&lwork,-1,info);
      double* work = memory_allocate(sizeof(double)*(size_t)(lwork+1));
      dgeev_y(jobvl,jobz,n,mpt,n,valblock,vpt,mpt2,n,mpt3,n,work,(size_t)lwork,info);
      memory_release(work);
      for (size_t m=0; m<n; m++)
	if (vpt[m]!=0.0) {
	  if (vali==NULL) throw "Matrix has complex eigenvalues, but no array to receive imaginary part was passed, smat::ev";
	  for (size_t l=0; l<n; l++)
	    vali->block[m][l]=vpt[l];
	  goto validone;
	}
    validone: ;
    } else
      throw " Eigenvalues/vectors cannot be computed for antisymmetric matrix, smat::ev";
    memory_release_saved(bases);
    // sort eigensolutions
    bases = memory_save();
    double* vpt = memory_allocate(sizeof(double)*n);
    if (sort[0]=='A' || sort[0]=='a')
      for (size_t l=0; l<n; l++)
	vpt[l] = valblock[l];
    else if (sort[0]=='D' || sort[0]=='d')
      for (size_t l=0; l<n; l++)
	vpt[l] = -valblock[l];
    else if (sort[0]=='O' || sort[0]=='o') {
      if (vec == NULL) throw "Sorting of eigensolutions by overlap requested, but eigenvectors not calculated, smat::ev";
      for (size_t l=0; l<n; l++) {
	double tester=0;
	for (size_t m=0; l<n; l++)
	  if (std::fabs(vec->block[k][m+l*n]) > tester) {
	    tester=std::fabs(vec->block[k][m+l*n]);
	    vpt[l] = (double)m;
	  }
      }
    } else
      throw "Unknown sorting algorithm, smat::ev";

    // naive sort
    std::vector<long> map(n,-1);
    {
      std::vector<bool> chosen(n,false);
      for (size_t l=0; l<n; l++) {
	double tester=-(std::numeric_limits<double>::max())/4;
	for (size_t m=0; m<n; m++) {
	  if (chosen[m] || vpt[m] < tester) continue;
	  tester = vpt[m];
	  map[l]=m;
	}
	chosen[map[l]]=true;
      }
    }
    memory_release(vpt);

    // sort eigenvalues
    vpt = memory_allocate(sizeof(double)*n);
    double* mpt = valblock;
    for (size_t m=0; m<n; m++)
      vpt[m] = mpt[map[m]];
    for (size_t m=0; m<n; m++)
      mpt[m] = vpt[m];
    if (vali != NULL) {
      mpt = valblock;
      for (size_t m=0; m<n; m++)
	vpt[m] = mpt[map[m]];
     for (size_t m=0; m<n; m++)
       mpt[m] = vpt[m];
   }
   memory_release(vpt);
   // sort eigenvectors
  if (vec != NULL) {
    double* mpt = vec->block[k];
    double* mpt2 = memory_allocate(sizeof(double)*n*n);
    for (size_t l=0; l<n; l++) {
      size_t iphase=0;
      for (size_t m=0; m<n; m++)
	if (std::fabs(mpt[m+map[l]*n]) > std::fabs(mpt[iphase+map[l]*n])) iphase=m;
      if (mpt[iphase+l*n]>0)
	for (size_t m=0; m<n; m++)
	  mpt2[m+l*n] = mpt[m+map[l]];
      else
	for (size_t m=0; m<n; m++)
	  mpt2[m+l*n] = -mpt[m+map[l]];
    }
    for (size_t m=0; m<n*n; m++)
      mpt[m] = mpt2[m];
    memory_release(mpt2);
  }
  if (vecl != NULL) {
    double* mpt = vecl->block[k];
    double* mpt2 = memory_allocate(sizeof(double)*n*n);
    for (size_t l=0; l<n; l++) {
      size_t iphase=0;
      for (size_t m=0; m<n; m++)
	if (std::fabs(mpt[m+map[l]*n]) > std::fabs(mpt[iphase+map[l]*n])) iphase=m;
      if (mpt[iphase+l*n]>0)
	for (size_t m=0; m<n; m++)
	  mpt2[m+l*n] = mpt[m+map[l]];
      else
	for (size_t m=0; m<n; m++)
	  mpt2[m+l*n] = -mpt[m+map[l]];
    }
    for (size_t m=0; m<n*n; m++)
      mpt[m] = mpt2[m];
    memory_release(mpt2);
  }
  memory_release_saved(bases);
  }
  if (vec==NULL) delete vvpt;
  return info;
}
smat smat::desymmetrise() const
{
  std::vector<std::vector<size_t> > newdim;
  std::vector<size_t> dim(8,0);
  for (int i=0; i<2; i++) {
    newdim.push_back(dim);
    for (int j=0; j<8; j++) newdim[i][0] += this->dimensions[i][j];
  }
  smat newmat(newdim,this->parity,0,this->rank,1,NULL,this->description);
  newmat.ensure_buffer();
  newmat.zero();
  double* to=newmat.block[0];
  size_t nrnew=newdim[0][0];
  if (rank == 1 )
    if (symmetry==0)
      for (size_t i=0; i<newmat.dimensions[0][0]; i++)
	to[i]=this->block[0][i];
    else {
    to += block_offset(symmetry);
    for (size_t i=0; i<dimensions[0][symmetry]; i++)
      to[i]=this->block[symmetry][i];
    }
  else { // rank==2
    if (parity==0 || symmetry!=0) {
      for (int sr=0; sr<8; sr++) {
	size_t nr = dimensions[0][sr];
	int sc=sr^symmetry;
	if (sc>sr) continue;
	size_t nc = dimensions[0][sc];
	double* from=this->block[sr];
	for (size_t c=0; c<nc; c++)
	  for (size_t r=0; r<nr; r++)
	    to[r+c*nrnew] = from[r+c*nr];
	to += nrnew*nc + nr;
      }
    }
    else { // symmetry==0 && parity !=0
      size_t cabs=0;
      for (int sr=0; sr<8; sr++) {
	size_t nr = dimensions[0][sr];
	double* from=this->block[sr];
	for (size_t c=0; c<nr; c++) {
	  for (size_t r=0; r<=c; r++)
	    to[r] = from[r];
	  from += c+1;
	  cabs++;
	  to += cabs;
	}
	to += nr;
      }
    }
  }
  return newmat;
}

#include <cmath>
extern "C" {
  void symmetry_matrix_module_test_c(int printlevel) {
#ifdef MOLPRO
    memory_module_test();
    try {
    smat om("YY",1,0);
    double t=0;
    for (unsigned int ks=0; ks<8; ks++)
      for (size_t k=0; k<om.dimensions[0][ks]; k++) {
	for (size_t l=0; l<k; l++)
	  om.block[ks][k*(k+1)/2+l]=(k+l+2)*0.001;
	om.block[ks][(k+1)*(k+2)/2-1]=(k+1)%3+(k+1)*0.4;
	t+= om.block[ks][(k+1)*(k+2)/2-1];
      }
    if (printlevel > 1) xout <<om.str("om")<<std::endl;
    if (printlevel > 1) xout <<"trace: "<<t<<" "<<om.trace()<<std::endl;
    if (std::fabs(t-om.trace()) > 1e-48) throw "trace error, symmetry_matrix_module_test_c";

    for (size_t k=0; k<om.size(); k++) om.buffer[k]=(k+1)*0.01;
    if (printlevel > 1) xout <<om.str("om")<<std::endl;
    smat oms(om,0);
    oms.copy(om,0);
    if (printlevel > 1) xout <<oms.str("oms")<<std::endl;

    for (size_t k=0; k<oms.size(); k++) oms.buffer[k]=0;
    for (unsigned int ks=0; ks<8; ks++) {
      size_t n=oms.dimensions[0][ks];
      for (size_t k=0; k<n; k++) oms.block[ks][k*(n+1)]=1;
      if (n<2) continue;
      double theta=.01*(ks+1);
      oms.block[ks][n+1]=oms.block[ks][0]=std::cos(theta);
      oms.block[ks][1]=std::sin(theta);
      oms.block[ks][n]=-oms.block[ks][1];
    }
    if (printlevel > 1) xout <<oms.str("oms")<<std::endl;


    smat om2(oms);
    if (printlevel > 1) xout <<om2.str("om2 after creation")<<std::endl;
    om2.gemm(oms,oms,'n','t',(double)1,(double)0);
    if (printlevel > 1) xout <<om2.str("om2")<<std::endl;

    smat ov("Y",0,-1,1,NULL,"Eigenvalues");
    if (printlevel > 1) xout << om.str("Matrix to be diagonalised");
    om.ev(ov);
    if (printlevel > 1) xout <<ov.str("eigenvalues")<<std::endl;

    smat us = oms.desymmetrise();
    if (printlevel > -1) xout <<oms.str("oms")<<std::endl;
    if (printlevel > -1) xout <<us.str("oms desymmetrised")<<std::endl;

    smat symmetric("  ",1,0); symmetric.copy(oms,1);
    if (printlevel > -1) xout <<symmetric.str("symmetric")<<std::endl;
    if (printlevel > -1) xout <<symmetric.desymmetrise().str("symmetric desymmetrised")<<std::endl;


  }
    catch(char const* msg) {
      xout <<msg<<std::endl;
    }
#endif
  }
}
