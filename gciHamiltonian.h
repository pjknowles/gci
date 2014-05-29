#ifndef GCIHAMILTONIAN_H
#define GCIHAMILTONIAN_H
#include "gci.h"
#include "FCIdump.h"
#include "gciSymmetrySpace.h"
#include "gciDeterminant.h"
#include "gciOrbitalSpace.h"
#include <string>
#include <vector>

namespace gci {
/**
 * @brief
 * Class holds hamiltonian operator for FCI or other calculation
 *
 */
class Hamiltonian : public OrbitalSpace
{
public:

  /*!
 \brief construct Hamiltonian object

 \param filename : if present, call load
*/
  Hamiltonian(std::string filename="");
  /*!
     \brief construct Hamiltonian object

     \param dump : if present, call load
    */   Hamiltonian(FCIdump* dump);
  /*!
   * \brief copy constructor
   * \param source
   * \param forceSpinUnrestricted whether to force converstion to a UHF object
   */
  Hamiltonian(const Hamiltonian &source, const bool forceSpinUnrestricted=false);
  ~Hamiltonian();
  void load(std::string filename="FCIDUMP", int verbosity=0); /**< \brief load integrals from FCIDUMP */
  void load(FCIdump* dump, int verbosity=0); /**< \brief load integrals from FCIDUMP */
  void unload(); /**< \brief destroy loaded integrals */
  /*!
     * \brief Construct a printable representation of the hamiltonian
     * \param verbosity how much information to include
     * \return printable representation of the hamiltonian
     */
  std::string str(int verbosity=0) const;
  bool loaded;  /**< \brief whether the integrals are loaded */
  double coreEnergy; /**< \brief core energy */
  std::vector<double> *integrals_a;  /**< \brief point to aa integrals */
  std::vector<double> *integrals_b; /**< \brief point to bb integrals */
  std::vector<double> *integrals_aa; /**< \brief point to aaaa integrals */
  std::vector<double> *integrals_ab; /**< \brief point to aabb integrals */
  std::vector<double> *integrals_bb; /**< \brief point to bbbb integrals */
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
     * \brief Generate a new object containing the Fock hamiltonian corresponding to a given reference determinant
     * \param reference the reference Slater determinant
     * \return  the Fock hamiltonian
     */
  Hamiltonian FockHamiltonian(const Determinant& reference) const;

  Hamiltonian& operator+=(const Hamiltonian &other); ///< add another Hamiltonian
  Hamiltonian& operator-=(const Hamiltonian &other); ///< subtract another Hamiltonian
  Hamiltonian& operator*=(const double factor); ///< scale by a factor
  /*!
   * \brief Construct the same-spin perturbed Hamiltonian
   * \param reference the zero-order Slater determinant
   * \return  the same-spin Hamiltonian
   */
  Hamiltonian sameSpinHamiltonian(const Determinant &reference) const;

private:
  size_t ijSize;
  size_t ijklSize;
  Hamiltonian& plusminusOperator(const Hamiltonian &other, const char operation='+');
  void plusminusEqualsHelper(std::vector<double>*& me, std::vector<double> * const &other, const char operation='+');
  void starEqualsHelper(std::vector<double> *&me, const double factor);
};
Hamiltonian operator+(const Hamiltonian &h1, const Hamiltonian &h2); ///< add two Hamiltonian objects
Hamiltonian operator-(const Hamiltonian &h1, const Hamiltonian &h2); ///< subtract two Hamiltonian objects
Hamiltonian operator*(const Hamiltonian &h1, const double factor); ///< return a factor times the Hamiltonian object

}

using namespace gci;

#endif // GCIHAMILTONIAN_H
