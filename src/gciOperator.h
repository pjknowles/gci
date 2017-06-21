#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H
#include "Operator.h"
#include "FCIdump.h"
#include "Eigen/Dense"
#include "gciDeterminant.h"
#include "gciOrbitalSpace.h"


namespace gci {

  class Operator : public SymmetryMatrix::Operator
  {
  public:
    using SymmetryMatrix::Operator::Operator;
    /*!
     * \brief Construct an operator from its symmetry-block dimensions.
     * \param dimensions Numbers of one-particle functions in each symmetry block
     * for each of the four two-particle indices (bra, bra, ket, ket).
     * \param rank The rank (0, 1 or 2) of the operator.
     * \param uhf Whether the underlying 1-particle spaces are different for alpha and beta spin.
     * \param hermiticity Conjugation symmetry, [ai] -> [ia], [bj] -> [jb]: 0=none, 1=symmetric, -1=antisymmetric, default 0.
     * \param exchange For each of bra, ket, whether there is an interchange symmetry between two
     * indices in the two-particle part of the operator: [abij] -> [baij]: 0=none, 1=symmetric, -1=antisymmetric, default -1.
     * \param symmetry Overall symmetry of operator (0-7).
     * \param description A string describing the object.
     */
    explicit Operator(const dims_t dimensions,
                      int rank=2,
                      bool uhf=false,
                      const std::vector<int> hermiticity={1,1},
                      const std::vector<int> exchange={-1,-1},
                      unsigned int symmetry=0,
                      std::string description="")
      : SymmetryMatrix::Operator(dimensions,rank,uhf,hermiticity,exchange,symmetry,description)
    {

    }
    /*!
     * \brief Construct a real hermitian fermionic operator
     * \param dimension Numbers of orbitals in each symmetry block.
     * \param rank The rank (0, 1 or 2) of the operator.
     * \param uhf Whether the underlying 1-particle spaces are different for alpha and beta spin.
     * \param symmetry Overall symmetry of operator (0-7).
     * \param description A string describing the object.
     */
    explicit Operator(const dim_t dimension,
                      int rank=2,
                      bool uhf=false,
                      unsigned int symmetry=0,
                      std::string description="")
      :  Operator(dims_t{dimension,dimension,dimension,dimension}, rank, uhf, std::vector<int>{1,1}, std::vector<int>{-1,-1}, symmetry, description) { }

    /*!
     * \brief Copy constructor. A complete (deep) copy is made.
     * \param source Object to be copied.
     */
    Operator(const Operator& source) : SymmetryMatrix::Operator(source) {}

    /*!
     * \brief Assigment operator
     * \param source Object to be copied.
     * \return A reference to this.
     */
    Operator& operator=(const Operator& source)
    {
      SymmetryMatrix::Operator::operator =(source);
      return *this;
    }

    /*!
     * \brief Construct an object from what is produced by bytestream(). If the bytestream
     * contains data, it will be loaded, otherwise the contents of the object are undefined,
     * and only the dimensions and parameters are loaded.
     * \param dump The raw buffer of a bytestream produced by bytestream()
     */
    static Operator construct(const char *dump);
    /*!
     * \brief Construct an object from what is produced by bytestream(). If the bytestream
     * contains data, it will be loaded, otherwise the contents of the object are undefined,
     * and only the dimensions and parameters are loaded.
     * \param bs The bytestream produced by bytestream()
     */
    static Operator construct(const class bytestream& bs) { return construct((const char*)&(bs.data()[0])); }
    /*!
     * \brief Obtain a reference to 1-particle matrix elements.
     * \param spinUp alpha or beta spin.
     * \return
     */
    /*!
     * \brief Construct an object from an FCIdump. If the FCIdump
     * contains data, it will be loaded, otherwise the contents of the object are undefined,
     * and only the dimensions and parameters are loaded.
     * \param dump The raw buffer of a FCIdump.
     */
    static Operator construct(const FCIdump& dump);

    /*!
     * \brief Construct an operator templated on this, but with a special specification
     * \param special
     * - "Q" projector onto space containing satellite orbital
     * - "P" 1-Q
     * \param forceSpinUnrestricted whether to force conversion to a UHF object
     */
    Operator* projector(const std::string special, const bool forceSpinUnrestricted) const;

    /*!
       * \brief int1 Generate array of diagonal one-electron integrals
       * \param spin positive for alpha, negative for beta
       * \return one-dimensional array with h(i,i) at i-1
       */
    Eigen::VectorXd int1(int spin) const;

    /*!
       * \brief intJ Generate array of two-electron exchange integrals
       * \param spini positive for alpha, negative for beta, first index
       * \param spinj positive for alpha, negative for beta, second index
       * \return array with (ii|jj)
       */
    Eigen::MatrixXd intJ(int spini, int spinj) const;
    /*!
       * \brief intK Generate array of two-electron Coulomb integrals
       * \param spin positive for alpha, negative for beta
       * \return array with (ij|ji)
       */
    Eigen::MatrixXd intK(int spin) const;

    std::vector<unsigned int> orbital_symmetries() const { return m_orbital_symmetries;}

    /*!
     * \brief Build a Fock operator from the density arising from a single Slater determinant
     * \param reference The Slater determinant
     * \param description Descriptive text
     */
    Operator fockOperator(const Determinant& reference, const std::string description="Fock") const;

    /*!
     * \brief Build a same-spin operator from the density arising from a single Slater determinant
     * \param reference The Slater determinant
     * \param description Descriptive text
     */
    Operator sameSpinOperator(const Determinant &reference, const std::string description="Same Spin Hamiltonian") const;

   std::vector<OrbitalSpace> m_orbitalSpaces;
  private:
  std::vector<unsigned int> m_orbital_symmetries;

    /*!
     * \brief offset Return the number of orbitals of the same symmetry before the given one.
     * \param i Absolute number (starting with 1) of the orbital.
     * \return
     */
    size_t offset(unsigned int i) const {
      return i < 1 ? 0 : std::count(m_orbital_symmetries.begin(),m_orbital_symmetries.begin()+i-1,m_orbital_symmetries[i-1]);
    }
  public:
    const FCIdump* m_fcidump; // temporary, to support migration
  };

}

#endif // GCIOPERATOR_H
