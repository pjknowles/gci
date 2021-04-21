#include "Problem.h"

molpro::gci::Problem::Problem(const molpro::Operator& hamiltonian) : m_hamiltonian(hamiltonian) {}

void molpro::gci::Problem::action(const CVecRef<container_t>& parameters, const VecRef<container_t>& actions) const {
  for (size_t k = 0; k < parameters.size(); k++) {
    const auto& v = parameters[k].get();
    auto& a = actions[k].get();
    a.fill(0);
    a.operatorOnWavefunction(m_hamiltonian, v);
  }
}

bool molpro::gci::Problem::diagonals(container_t& d) const {
  d.diagonalOperator(m_hamiltonian);
  return true;
}

void molpro::gci::Problem::p_action(const std::vector<std::vector<value_t>>& p_coefficients,
                                    const CVecRef<std::map<size_t, container_t::value_type>>& pparams,
                                    const VecRef<container_t>& actions) const {}

std::vector<double>
molpro::gci::Problem::pp_action_matrix(const std::vector<std::map<size_t, container_t::value_type>>& pparams) const {
  return std::vector<double>();
}

void molpro::gci::Problem::precondition(const VecRef<container_t>& action, const std::vector<double>& shift,
                                        const container_t& diagonals) const {
  auto dlb = diagonals.distr_buffer->local_buffer();
  auto offset = dlb->start();
  for (int k = 0; k < action.size(); k++) {
    auto& a = action[k].get();
    auto alb = a.distr_buffer->local_buffer();
    for (int i = 0; i < alb->size(); i++)
      (*alb)[i] /= ((*dlb)[i] - shift[k] + 1e-15);
  }
}
