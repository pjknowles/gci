#ifndef GCI_GCIMIXEDOPERATORSECONDQUANT_H
#define GCI_GCIMIXEDOPERATORSECONDQUANT_H

#include <FCIdump.h>

#include "gciVibOperator.h"
#include "gciOptions.h"
#include "gciPersistentOperator.h"
#include <hdf5.h>


namespace gci {


/*!
 * @brief Mixed electron-vibration operator in a fully second quantized form.
 *
 * sqH_el = pd q h_pq + 2e- + pd q E_ij^A h_{pq,ij}^A
 *
 * @note Only 1MC operators are currently implemented
 *
 * @note Also includes polaritonic Hamiltonian for coupling of cavity photon to
 * the electrons via molecular dipole moment operator
 *   H = HO(cavity) + (g.D) (a_dg + a)  + (g.D)^2 / w
 *
 * Large electronic operators from the mixed-tensor are stored on an hdf5-file.
 * This is handled by ``PersistentOperator``, which requires that the hdf5-file be open
 * externally and remain open while ``PersistentOperator`` is in use.
 *
 *
 */
class MixedOperatorSecondQuant {
protected:
    std::string m_fcidump_f; //!< base name for dump files storing electronic operators of the mixed Tensor
    std::string hdf5_fname; //!< file name for hdf5 where large eletronic operators are stored
    bool restart; //!< operator is already stored on the hdf5, no overwriting is done
    hid_t hid_file; //!< id of the hdf5 file
    bool hdf5_file_owner; //!< whether it owns the hdf5 file
    std::string m_description; //!< desctiption of the Hamiltonian
public:
    using hel_t = PersistentOperator;
    const int nMode; //!< Number of vibrational modes
    const int nModal; //!< Number of modals per mode (for now assumed the same for each mode)
    VibOperator<double> Hvib;//!< Purely vibrational term
    std::map<std::string, hel_t> elHam; //!< Purely electronic terms
    //! Purely electronic Hamiltonian that has to be applied twice
    //! e.g. dipole moment operator
    std::map<std::string, hel_t> elHam2;
    std::map<std::string, VibOperator<hel_t >> mixedHam;//!< Mixed electronic-vibrational terms

    explicit MixedOperatorSecondQuant(const Options &options);
    MixedOperatorSecondQuant(const MixedOperatorSecondQuant &other, const std::string &include_ham);
    ~MixedOperatorSecondQuant();
    /*!
     * @brief Checks if bra and ket vibrational basis are connected by the mixed Hamiltonian
     */
    bool connected(const HProduct &bra, const HProduct &ket) const;

    std::string description() {return m_description;}

public:

    /*!
     * @brief From FCIdump file generates antisymmetric electronic operator
     */
    static SymmetryMatrix::Operator constructOperatorAntisymm1el(const FCIdump &dump, bool collective = true);

    static SymmetryMatrix::Operator constructK(const FCIdump &dump, bool collective = true);

    static SymmetryMatrix::Operator constructD(const FCIdump &dump, bool collective = true);

    void initializeHel(const FCIdump &fcidump, bool h0_only = false);

    void initializeLambda(const FCIdump &fcidump);

    void initializeK(const FCIdump &fcidump);

    void initializeD(const FCIdump &fcidump);
};

/*!
 * @brief Constructs purely vibrational operator
 * @note I'm assuming it's 1MC only
 * @param Hvib empty vibrational operator
 * @param fcidump_name Name of the main fcidump file
 * @param nmode number of modes
 * @param nmodal  number of modals
 */
void constructHvib(VibOperator<double> &Hvib, const std::string &fcidump_name, int nmode, int nmodal);


/*!
 * @brief Constructs harmonic oscillator operator
 * @param Hvib empty vibrational operator
 * @param freq harmonic frequencies for each mode
 * @param nmode number of modes
 * @param nmodal number of modal basis functions
 */
void constructHvib(VibOperator<double> &Hvib, int nmode, int nmodal, std::vector<double> freq);

/*!
 * @brief constructs coupling Hamiltonian of cavity photon to dipole moment of electronic field
 *        and dipole moment self-interaction
 *        mixedHam = (g.D) (a_dg + a)
 *        elHam2 = (g.D)^2 / w
 * @param elHam2
 * @param mixedHam
 * @param fcidump_f
 * @param gamma
 * @param freq
 */
void constructDMcoupling(std::map<std::string, MixedOperatorSecondQuant::hel_t> &elHam2,
                         std::map<std::string, VibOperator<MixedOperatorSecondQuant::hel_t >> &mixedHam,
                         const std::string &fcidump_f, const std::vector<double> &gamma,
                         const std::vector<double> &freq, int nmode, int nmodal);

} // namespace gci

#endif //GCI_GCIMIXEDOPERATORSECONDQUANT_H
