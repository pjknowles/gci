#ifndef GCI_GCIITERATIVESOLVERUTIL_H
#define GCI_GCIITERATIVESOLVERUTIL_H
#include <vector>

#include "gci.h"
namespace gci {
namespace run {

//TODO This is a bad design! Get rid of global variables
double _lastEnergy;
double _mu;
double _residual_q;
bool parallel_stringset;

/*!
 * @brief
 */
template<class t_Wavefunction>
using ParameterVectorSet = std::vector<t_Wavefunction>;

using Pvector = std::map<size_t, double>;


//TODO Why do we even need it? This can easily be a member function for each iterative solver. There are only about 10 lines that are actually relevant!
/*!
 * @todo Write documentation
 * @brief
 * @tparam t_Wavefunction
 * @tparam t_Operator
 * @tparam ParameterVectorSet
 */
template<class t_Wavefunction, class t_Operator>
struct residual {
protected:
    using value_type = typename t_Wavefunction::value_type;
    const t_Operator &m_hamiltonian;
    const bool m_subtract_Energy;
    const t_Operator *m_Q;
public:
    /*!@todo Write documentation
     * @brief
     * @param hamiltonian
     * @param subtract_Energy
     * @param Q
     * @param parallel_stringset
     */
    residual(const t_Operator &hamiltonian, bool subtract_Energy, t_Operator *Q = nullptr,
             bool parallel_stringset = false) :
            m_hamiltonian(hamiltonian), m_subtract_Energy(subtract_Energy), m_Q(Q) { }

    /*!@todo Improve documentation
     * @brief
     * @param psx Vectors from a parameter set???
     * @param outputs The residual?
     * @param active Flags which vectors are actively
     * @param append
     */
    void operator()(const ParameterVectorSet<t_Wavefunction> &psx, ParameterVectorSet<t_Wavefunction> &outputs,
                    std::vector<bool> active, bool append = false) const {
        for (size_t k = 0; k < psx.size(); k++) {
            const t_Wavefunction &x = psx[k];
            t_Wavefunction &g = outputs[k];
            if (not append)
                g.zero();
            if (active[k]) {
                auto prof = profiler->push("Hc");
                g.operatorOnWavefunction(m_hamiltonian, x, parallel_stringset);
            }
            if (m_subtract_Energy) {
                double cc = x.dot(x);
                double cg = x.dot(g);
                _lastEnergy = cg / cc;
                double epsilon = cg / cc;
                if (m_Q != nullptr) {
                    t_Wavefunction m(g);
                    m.zero();
                    m.operatorOnWavefunction(*m_Q, x);
                    double cm = x.dot(m);
                    double gm = g.dot(m);
                    _mu = cm == 0 ? 0 : (cg * cm - cc * gm) / (cm * cm - cm * cc);
                    epsilon = (cg - cm * _mu + cc * _mu * _residual_q) / (cc);
                    g.axpy(-_mu, m);
                    // FIXME idempotency constraint to follow
                    _lastEnergy = epsilon - _mu * _residual_q;
                }
                g.axpy(-_lastEnergy, x);
            }
        }
    }
};

template<class t_Wavefunction, class t_Operator>
struct Presidual {
private:
    const t_Operator &m_hamiltonian;
    const std::vector<Pvector> &m_P;
public:
    Presidual(const t_Operator &hamiltonian, const std::vector<Pvector> &P) : m_hamiltonian(hamiltonian), m_P(P) { }
    void
    operator()(const std::vector<std::vector<double> > &Pcoeff, ParameterVectorSet<t_Wavefunction> &outputs) const {
        for (size_t k = 0; k < Pcoeff.size(); k++) {
            assert(m_P.size() == Pcoeff[k].size());
            t_Wavefunction &g = outputs[k];
            t_Wavefunction w(g);
            w.m_sparse = true;
            for (size_t i = 0; i < m_P.size(); i++) {
                w.buffer_sparse.insert({m_P[i].begin()->first, Pcoeff[k][i]});
            }
            auto prof = profiler->push("HcP");
            g.operatorOnWavefunction(m_hamiltonian, w);
        }
    }
};

template<class t_Wavefunction>
class updater {
public:
    updater(const t_Wavefunction &diagonals, bool subtractDiagonal)
            : m_diagonals(diagonals), m_subtractDiagonal(subtractDiagonal) { }
protected:
    using value_type = typename t_Wavefunction::value_type;
    const t_Wavefunction &m_diagonals;
    const bool m_subtractDiagonal;
public:
    /*!@todo Add more detailed description
     * @brief For linear case (matrix diagonalisation) applies a preconditioner
     * @param psc Set of current solutions
     * @param psg Set of residuals for non-linear optimization or actions from matrix vector operation in linear diagonalisation
     * @param shift Per state shift to use in the preconditioner
     * @param append Whether to add  gw / (diagonals - shift) to cw or overwrite it
     */
    void operator()(ParameterVectorSet<t_Wavefunction> &psc,
                    const ParameterVectorSet<t_Wavefunction> &psg,
                    std::vector<value_type> shift = {},
                    bool append = false) const {
        std::vector<value_type> shifts = shift;
        for (size_t state = 0; state < psc.size(); state++) {
            if (m_subtractDiagonal)
                shifts[state] -= m_diagonals.at(m_diagonals.minloc(state + 1));
            t_Wavefunction &cw = psc[state];
            const t_Wavefunction &gw = psg[state];
            //FIXME change to a threshold
            if (shift[state] == 0) {
                //FIXME shouldn't this be division?
                cw.times(&gw, &m_diagonals);
            } else {
                shifts[state] += 2 * std::numeric_limits<value_type>::epsilon()
                                 * std::max<value_type>(1, std::abs(
                        m_diagonals.at(m_diagonals.minloc(state + 1)))); // to guard against zero
                cw.divide(&gw, &m_diagonals, shifts[state], append, true);
            }
        }
    }
};

}  // namespace run
}  // namespace gci
#endif //GCI_GCIITERATIVESOLVERUTIL_H
