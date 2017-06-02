#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H
#include "gci.h"
#include "FCIdump.h"
#include "gciSymmetrySpace.h"
#include "gciDeterminant.h"
#include "gciOrbitalSpace.h"
#include "SMat.h"
#include "Operator.h"
#include <string>
#include <vector>
#include <climits>

namespace gci {
/**
 * @brief
 * Class holds hamiltonian or other operator for FCI or other calculation
 *
 */
class OldOperator : public OrbitalSpace
{
public:

  /*!
   * \brief Construct from SymmetryMatrix::Operator
   * \param source
   */
  OldOperator(const SymmetryMatrix::Operator& source);

  /*!
 \brief construct Operator object

 \param filename : if present, call load
*/
  OldOperator(std::string filename="");
  /*!
     \brief construct Operator object

     \param dump : if present, call load
    */   OldOperator(FCIdump* dump);
  /*!
   * \brief copy constructor
   * \param source
   */
  OldOperator(const OldOperator& source);
  /*!
   * \brief general copy constructor
   * \param source
   * \param forceSpinUnrestricted whether to force conversion to a UHF object
   * \param oneElectron whether to copy the 1-electron part of source
   * \param twoElectron whether to copy the 2-electron part of source
   */
  explicit OldOperator(const OldOperator &source, const bool forceSpinUnrestricted, const bool oneElectron=true, const bool twoElectron=true);
  /*!
   * \brief Construct an operator templated on another, but with a special specification
   * \param source
   * \param special
   * \param forceSpinUnrestricted whether to force conversion to a UHF object
   * - "Q" projector onto space containing satellite orbital
   * - "P" 1-Q
   */
  explicit OldOperator(const std::string special, const OldOperator &source, const bool forceSpinUnrestricted=false);
  virtual ~OldOperator();
  void load(std::string filename="FCIDUMP", int verbosity=0); /**< \brief load integrals from FCIDUMP */
  void load(FCIdump* dump, int verbosity=0); /**< \brief load integrals from FCIDUMP */
  void unload(); /**< \brief destroy loaded integrals */
  /*!
     * \brief Construct a printable representation of the operator
     * \param verbosity how much information to include
     * \param columns Page width
     * \return printable representation of the operator
     */
  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;
  bool loaded;  /**< \brief whether the integrals are loaded */
  double coreEnergy; /**< \brief core energy */
  std::vector<double> *integrals_a;  /**< \brief point to aa integrals */
  std::vector<double> *integrals_b; /**< \brief point to bb integrals */
  std::vector<double> *integrals_aa; /**< \brief point to aaaa integrals */
  std::vector<double> *integrals_ab; /**< \brief point to aabb integrals */
  std::vector<double> *integrals_bb; /**< \brief point to bbbb integrals */
  std::vector<double> *bracket_integrals_a; /**< \brief point to aa integrals in bra-ket form, antisymmetrised */
  std::vector<double> *bracket_integrals_b; /**< \brief point to bb integrals in bra-ket form, antisymmetrised */
 private:
 public:
  std::vector<double> *bracket_integrals_aa; /**< \brief point to aaaa integrals in bra-ket form, antisymmetrised */
  std::vector<double> *bracket_integrals_ab; /**< \brief point to abab integrals in bra-ket form */
  std::vector<double> *bracket_integrals_bb; /**< \brief point to bbbb integrals in bra-ket form, antisymmetrised */
  size_t basisSize;///< \brief size of orbital basis set
  /*!
     * \brief calculate canonical index of a pair of orbitals.
     * \param i Absolute number (starting with 1) of first orbital.
     * \param j Absolute number (starting with 1) of second orbital.
     * \return Number of orbital pairs of the same symmetry canonically before ij.
     */
  size_t int1Index(unsigned int i, unsigned int j) const;
  /*!
     * \brief calculate canonical address of a 2-electron integral
     * \param i Absolute number (starting with 1) of first orbital.
     * \param j Absolute number (starting with 1) of second orbital.
     * \param k Absolute number (starting with 1) of third orbital.
     * \param l Absolute number (starting with 1) of fourth orbital.
     * \return
     */
  size_t int2Index (unsigned int i, unsigned int j, unsigned int k, unsigned int l) const;

  /*!
     * \brief int1 Generate array of diagonal one-electron integrals
     * \param spin positive for alpha, negative for beta
     * \return one-dimensional array with h(i,i) at i-1
     */
  std::vector<double> int1(int spin) const;


  /*!
     * \brief intJ Generate array of two-electron exchange integrals
     * \param spini positive for alpha, negative for beta, first index
     * \param spinj positive for alpha, negative for beta, second index
     * \return one-dimensional array with (ii|jj) at i-1 + (j-1)*basisSize
     */
  std::vector<double> intJ(int spini, int spinj) const;
  /*!
     * \brief intK Generate array of two-electron Coulomb integrals
     * \param spin positive for alpha, negative for beta
     * \return one-dimensional array with (ij|ji) at i-1 + (j-1)*basisSize
     */
  std::vector<double> intK(int spin) const;

  /*!
     * \brief Generate a new object containing the Fock operator corresponding to a given reference determinant
     * \param reference the reference Slater determinant
     * \return  the Fock operator
     */
  OldOperator FockOperator(const Determinant& reference) const;

  OldOperator& operator+=(const OldOperator &other); ///< add another Operator
  OldOperator& operator-=(const OldOperator &other); ///< subtract another Operator
  OldOperator& operator*=(const double factor); ///< scale by a factor
  OldOperator& operator=(const OldOperator &h1); ///< assignment overload
  /*!
   * \brief Construct the same-spin perturbed Operator
   * \param reference the zero-order Slater determinant
   * \return  the same-spin Operator
   */
  OldOperator sameSpinOperator(const Determinant &reference) const;

  /*!
   * \brief construct the BraKet form of the Operator
   * \param neleca if non-zero, incorporate 1-electron integrals into 2-electron for neleca electrons
   * \param nelecb if non-zero, incorporate 1-electron integrals into 2-electron for nelecb electrons
   */
  void constructBraKet(int neleca=0, int nelecb=0);
  /*!
   * \brief Delete the BraKet form of the Operator
   */
  void deconstructBraKet();

  void rotate1(std::vector<double>* integrals, const std::vector<double> *rot);
  void rotate2(std::vector<double>* integrals, const std::vector<double> *rot1, const std::vector<double> *rot2);
  void rotate(std::vector<double> const * rota, std::vector<double> const * rotb=NULL);
  void rotate(SMat const * rota, SMat const * rotb=NULL);

private:
  const SymmetryMatrix::Operator& m_Operator; ///< reference to the underlying operator
  size_t ijSize;
  size_t ijklSize;
  OldOperator& plusminusOperator(const OldOperator &other, const char operation='+');
  void plusminusEqualsHelper(std::vector<double>*& me, std::vector<double> * const &other, const char operation='+');
  void starEqualsHelper(std::vector<double> *&me, const double factor);
  /*!
   * \brief general copy of another object
   * \param source
   * \param forceSpinUnrestricted whether to force converstion to a UHF object
   * \param oneElectron whether to copy the 1-electron part of source
   * \param twoElectron whether to copy the 2-electron part of source
   */
  void _copy(const OldOperator &source, const bool forceSpinUnrestricted=false, const bool oneElectron=true, const bool twoElectron=true);
};
OldOperator operator+(const OldOperator &h1, const OldOperator &h2); ///< add two Operator objects
OldOperator operator-(const OldOperator &h1, const OldOperator &h2); ///< subtract two Operator objects
OldOperator operator*(const OldOperator &h1, const double factor); ///< return a factor times the Operator object

}

using namespace gci;

#endif // GCIOPERATOR_H
