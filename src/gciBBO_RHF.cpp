#include "gciRun.h"
#include "gciBBO_RHF.h"
#include "gciOperatorBBO.h"

using namespace gci;

void Run::BBO_RHF(const State &prototype) {
  int maxIterations;
  double thresh;
  maxIterations = options.parameter("HF_MAX_ITER", std::vector<int>{50}).at(0);
  thresh = options.parameter("HF_E_THRESH", std::vector<double>{1.0e-8}).at(0);
  OperatorBBO molHam(options); // molecular Hamiltonian
//  xout << molHam << std::endl;
  //time to do HF
  // For now there is no need to be very smart about it.
  // Just copy the RHF looping structure and split the body into an electronic
  // and bosonic components.
  // The MO vectors will be handled on the local level.
  profiler->start("BBO_RHF");
  std::vector<int> symmetries;
  for (const auto &s : prototype.orbitalSpace->orbital_symmetries) {
    symmetries.push_back(s + 1);
  }
  dim_t dim;
  for (const auto s: *prototype.orbitalSpace) dim.push_back(s);
  Operator C(dim, symmetries, 1, false, prototype.symmetry, false, false, "MO");
  dim_t occ(8,0);
  std::vector<int> occInt(8, 0);
  occInt = options.parameter("OCC", std::vector<int>{8});
  SMat Csplice({dim, occ}, parityNone);
  SMat CspliceT(&Csplice, parityNone, 0, 2);
  SMat Cmat = C.O1(true);
  Cmat.setIdentity();
  Operator P(dim, symmetries, 1, false, prototype.symmetry, false, true, "density");
  vector<SMat> U; // vibrational modals
  for (int iMod = 0; iMod < molHam.m_nMode; ++iMod){
    dim.assign(8, 0);
    dim[molHam.m_symMode[iMod] - 1] = (unsigned int) molHam.m_nModal;
    U.push_back(SMat({dim, dim}, parityNone, molHam.m_symMode[iMod], "Modals for mode " + std::__cxx11::to_string(iMod + 1)));
  }
  double energy, energyPrev=0.0, err;
  std::vector<double> energy, energyPrev=0.0, err;
  std::cout << "Iter " << " E " << " dE " << std::endl;
  // Energy and convergence check first.
  for (int iIter = 0; iIter < maxIterations; ++iIter) {
    Csplice.splice(Cmat);
    CspliceT = Csplice;
    CspliceT.transpose();
    P.O1(true) = 2.0 * (Csplice * CspliceT);
// Electronic contribution
    Operator F = molHam.m_Hel.fock(P, true, "Fock operator");
    energy = molHam.m_Hel.m_O0 + 0.5 * ( P.O1() & (molHam.m_Hel.O1() + F.O1()));
    err = energy - energyPrev;
    std::cout << iIter << "    " << energy << "    " << err << std::endl;
    if (abs(err) < thresh) break;
    energyPrev = energy;
    SMat eigVal({dim}, parityNone, -1, "Eigenvalues");
    F.O1().ev(eigVal, &Cmat);
  }
  profiler->stop("RHF");
}