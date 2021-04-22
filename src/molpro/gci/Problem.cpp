#include "Problem.h"
#include <molpro/mpi.h>

molpro::gci::Problem::Problem(const molpro::Operator& hamiltonian, const State& prototype)
    : m_hamiltonian(hamiltonian), m_prototype(prototype) {}

void molpro::gci::Problem::action(const CVecRef<container_t>& parameters, const VecRef<container_t>& actions) const {
  for (size_t k = 0; k < parameters.size(); k++) {
    const auto& v = parameters[k].get();
//#ifdef HAVE_MPI_H
//    auto distribution = v.distr_buffer->distribution();
//    for (int rank = 0; rank < mpi::size_global(); rank++)
//      MPI_Bcast((void*)(v.buffer.data() + distribution.range(rank).first),
//                distribution.range(rank).second - distribution.range(rank).first, MPI_DOUBLE, rank, mpi::comm_global());
//#endif
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
                                    const VecRef<container_t>& actions) const {
  for (size_t k = 0; k < p_coefficients.size(); k++) {
    Wavefunction& g = actions[k];
    Wavefunction w(g);
    w.m_sparse = true;
    for (size_t i = 0; i < pparams.size(); i++) {
      assert(pparams[i].get().size() == 1);
      w.buffer_sparse.insert({pparams[i].get().begin()->first, p_coefficients[k][i]});
    }
    auto prof = profiler->push("HcP");
    g.operatorOnWavefunction(m_hamiltonian, w);
  }
}

std::vector<double>
molpro::gci::Problem::pp_action_matrix(const std::vector<std::map<size_t, container_t::value_type>>& pparams) const {
  auto prof = profiler->push("HPP");
  constexpr size_t NP = 0;
  const auto newNP = pparams.size();
  std::vector<double> addHPP(newNP * (newNP - NP), (double)0);
  for (size_t p0 = NP; p0 < newNP; p0++) {
    Wavefunction wsparse(m_prototype);
    wsparse.m_sparse = true;
    Wavefunction gsparse(m_prototype);
    gsparse.m_sparse = true;
    assert(pparams[0].size() == 1);
    wsparse.set(pparams[p0].begin()->first, (double)1);
    gsparse.operatorOnWavefunction(m_hamiltonian, wsparse);
    for (size_t p1 = 0; p1 < newNP; p1++) {
      auto jdet1 = pparams[p1].begin()->first;
      if (gsparse.buffer_sparse.count(jdet1))
        addHPP[p1 + (p0 - NP) * newNP] = gsparse.buffer_sparse.at(jdet1);
    }
  }

  return addHPP;
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
#ifdef HAVE_MPI_H
    auto distribution = a.distr_buffer->distribution();
    for (int rank = 0; rank < mpi::size_global(); rank++)
      MPI_Bcast((void*)(a.buffer.data() + distribution.range(rank).first),
                distribution.range(rank).second - distribution.range(rank).first, MPI_DOUBLE, rank, mpi::comm_global());
#endif
  }
}
