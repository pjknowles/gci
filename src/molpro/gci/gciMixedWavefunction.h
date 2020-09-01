#ifndef GCI_GCIMIXEDWAVEFUNCTION_H
#define GCI_GCIMIXEDWAVEFUNCTION_H

#include <mpi.h>
#include <vector>

#include "molpro/gci/gciHProductSet.h"
#include "molpro/gci/gciMixedOperator.h"
#include "molpro/gci/gciMixedOperatorSecondQuant.h"
#include "molpro/gci/gciWavefunction.h"

namespace molpro {
namespace gci {
namespace array {
class Array;
}

/*!
 * @brief Mixed bosonic and fermionic wavefunction represented as configuration expansion in the direct product
 * space of Slater Determenants and Hartree Product of modals. It implements storage of the wavefunction
 * coefficients and application of the Hamiltonitan ( r = H c).
 *
 * FCIdump parameters:
 *      NMODE
 *          -- number of vibrational modes
 *          default = 0
 *      NMODAL
 *          -- number of modal basis functions per mode (for ground state NMODAL=1)
 *          default = 1
 *      VIB_EXC_LVL
 *          -- maximum number of modes simultaneously excited in the wavefunction
 *          default = 1
 *
 * Vibrational basis is specified by mode coupling level, CIS, CISD etc, and maximum excitation level.
 *
 * Nomenclature:
 *      modal - one particle vibrational basis function
 *
 * The full wavefunction is stored in the global array buffer. It is ordered by the vibrational basis, with each
 * vibrational basis having a corresponding full electronic CI vector.
 * Psi = {psi_{p,q}_0, psi_{p,q}_1A, psi_{p,q}_2A,..., psi_{p,q}_1B, ..., psi_{p,q}_1A_1B, ...}
 * where psi_{p,q} is the electronic CI wavefunction and indices IA imply direct product with HO basis function I of
 * mode A.
 *
 * Prallelism
 * ----------
 * The full wavefunction is stored in a global array. Relevant chunks are copied into the Wavefunction buffer
 * for computation. The buffer is copied (or accumulated etc) back into GA after which it is released.
 *
 *
 */
class MixedWavefunction : virtual public Printable {
public:
  MPI_Comm m_child_communicator; //!< Communicator for children Wavefunction objects
  using value_type = double;
  using mapped_type = double;

protected:
  VibSpace m_vibSpace;    //!< Parameters defining the vibrational space of current wavefunction
  HProductSet m_vibBasis; //!< Vibrational basis for the full space of current wavefunction
  size_t m_elDim;         //!< Dimension of the electronic (slater determinant) space

  /*!
   * @brief Prototype electronic wavefunction
   *
   * It's buffer is populated with relevant section from GA before computation
   */
  Wavefunction m_prototype;

public:
  std::unique_ptr<array::Array> m_array; //!< array storing the actual buffer and performing linear algebra

  explicit MixedWavefunction(const Options &options, const State &prototype, MPI_Comm head_commun = mpi_comm_compute);

  MixedWavefunction(const MixedWavefunction &source);
  MixedWavefunction(const MixedWavefunction &source, int option);
  MixedWavefunction &operator=(const MixedWavefunction &source);

  ~MixedWavefunction();

  /*!
   * @brief Returns a wavefunction corresponding to vibrational product under offset
   * @param iVib Index to the vibrational product
   * @return Reference to a wavefunction under offset
   */
  Wavefunction wavefunctionAt(size_t iVib, MPI_Comm commun) const;

  VibSpace vibSpace() const { return m_vibSpace; }      ///< copy of the vibrational space
  size_t elDim() const { return m_elDim; }              ///< size of electronic Fock space
  size_t vibDim() const { return m_vibBasis.vibDim(); } ///< size of vibrational Fock space

  /*!
   * @brief Gets boundaries in GA for a block corresponding to electronic wavefunction under vibrational
   * index ``indKet``
   * @param iVib vibrational index of the electronic Wavefunction
   * @param lo index of the start of the block
   * @param hi index of the end of the block (inclusive)
   */
  static void ga_wfn_block_bound(int iVib, int *lo, int *hi, int dimension);
  static void copy_to_local(const MixedWavefunction &w, int iVib, Wavefunction &wfn);
  void put(int iVib, Wavefunction &wfn);
  void accumulate(int iVib, Wavefunction &wfn);

  /*!
   * \brief Add to this object the action of an operator on another wavefunction
   * \param ham Fully second quantized mixed Hamiltonian operator
   * \param w Other mixed Wavefunction
   * \param parallel_stringset whether to use parallel algorithm in StringSet construction
   * \param with_sync whether to syncronise the processes at the end of operation.
   *                  Calculations without a sync have to DivideTasks beforhand,
   *                  since it requires a sync.
   */
  void operatorOnWavefunction(const MixedOperatorSecondQuant &ham, const MixedWavefunction &w,
                              bool parallel_stringset = false, bool with_sync = true);

  /*!
   * @brief Set this object to the diagonal elements of the hamiltonian
   * \param ham Fully second quantized mixed Hamiltonian operator
   */
  void diagonalOperator(const MixedOperatorSecondQuant &ham, bool parallel_stringset = false);

  /*!
   * @brief Calculates vibrational density matrix
   * @return
   */
  std::vector<double> vibDensity();

  /*!
   * @brief Checks that the two wavefunctions are of the same electronic State and of the same dimension.
   */
  bool compatible(const MixedWavefunction &other) const;
  //! Checks that wavefunctions are compatible during Array operations
  //    bool compatible(const Array &other) const override;

  void settilesize(int a, int b, int c){};

  std::string str(int v = 0, unsigned int c = 0) const override { return ""; }

  [[nodiscard]] double dot(const MixedWavefunction &w) const;
  [[nodiscard]] double dot(const std::map<unsigned long int, double> &w) const;
  void axpy(double a, const std::map<unsigned long int, double> &w) const;
  void axpy(double a, const MixedWavefunction &w);
  //! allocates the array buffer
  void allocate_buffer();
  void sync() const;
  void zero();
  size_t size() const;
  double at(unsigned long i) const;
  void set(unsigned long i, double v);
  void scal(double a);
  void divide(const MixedWavefunction *y, const MixedWavefunction *z, double shift = 0, bool append = false,
              bool negative = false);
  std::vector<size_t> minlocN(int n) const;

}; // class MixedWavefunction

} // namespace gci
} // namespace molpro
#endif // GCI_GCIMIXEDWAVEFUNCTION_H
