#ifndef GCI_SRC_MOLPRO_GCI_PROBLEM_H_
#define GCI_SRC_MOLPRO_GCI_PROBLEM_H_
#include "gciWavefunction.h"
#include <molpro/Operator.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>

namespace molpro::gci {

class Problem : public molpro::linalg::itsolv::Problem<Wavefunction> {
  const molpro::Operator& m_hamiltonian;

public:
  using linalg::itsolv::Problem<Wavefunction>::container_t;
  using P = std::map<size_t, container_t::value_type>;
  using linalg::itsolv::Problem<container_t>::value_t;
  Problem(const Operator& hamiltonian);
  Problem() = delete;
  //  value_t residual(const R& parameters, R& residual) const override;
  void action(const CVecRef<container_t>& parameters, const VecRef<container_t>& actions) const override;
  bool diagonals(container_t& d) const override;
  std::vector<double> pp_action_matrix(const std::vector<P>& pparams) const override;
  void p_action(const std::vector<std::vector<value_t>>& p_coefficients, const CVecRef<P>& pparams,
                const VecRef<container_t>& actions) const override;
  void precondition(const VecRef<container_t>& action, const std::vector<double>& shift,
                    const container_t& diagonals) const override;
};
} // namespace molpro::gci

#endif // GCI_SRC_MOLPRO_GCI_PROBLEM_H_
