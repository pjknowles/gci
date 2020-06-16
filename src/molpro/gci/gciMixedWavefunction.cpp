#include "gciMixedWavefunction.h"
#include "gciHProductSet.h"

#include <ga.h>
#include <mpi.h>

namespace molpro {
namespace gci {

MixedWavefunction::MixedWavefunction(const Options &options, const State &prototype, MPI_Comm head_commun)
    : Array(head_commun), m_child_communicator(_sub_communicator),
      m_vibSpace(options.parameter("NMODE", 0), options.parameter("NMODAL", 1), options.parameter("VIB_EXC_LVL", 1)),
      m_vibBasis(m_vibSpace), m_elDim(0), m_prototype(prototype, m_child_communicator) {
  m_elDim = m_prototype.size();
  m_vibBasis.generateFullSpace();
  m_dimension = m_vibBasis.vibDim() * m_elDim;
  allocate_buffer();
}

MixedWavefunction::MixedWavefunction(const MixedWavefunction &source, int option)
    : Array(source), m_child_communicator(source.m_child_communicator), m_vibSpace(source.m_vibSpace),
      m_vibBasis(source.m_vibBasis), m_elDim(source.m_elDim), m_prototype(source.m_prototype) {}

Wavefunction MixedWavefunction::wavefunctionAt(size_t iVib, MPI_Comm commun) const {
  auto wfn = Wavefunction{m_prototype, 0, commun};
  wfn.allocate_buffer();
  if (!empty())
    copy_to_local(m_ga_handle, iVib, wfn);
  return wfn;
}

void MixedWavefunction::ga_wfn_block_bound(int iVib, int *lo, int *hi, int dimension) {
  lo[0] = iVib * dimension;
  hi[0] = lo[0] + dimension - 1;
}

void MixedWavefunction::copy_to_local(int ga_handle, int iVib, Wavefunction &wfn) const {
  auto p = profiler->push("copy_to_local");
  double *buffer = wfn.buffer.data();
  auto dimension = wfn.dimension;
  int lo, hi, ld = dimension;
  ga_wfn_block_bound(iVib, &lo, &hi, dimension);
  NGA_Get(ga_handle, &lo, &hi, buffer, &ld);
}

void MixedWavefunction::put(int iVib, Wavefunction &wfn) {
  auto p = profiler->push("put");
  double *buffer = wfn.buffer.data();
  auto dimension = wfn.dimension;
  int lo, hi, ld = dimension;
  ga_wfn_block_bound(iVib, &lo, &hi, dimension);
  Array::put(lo, hi, buffer);
}

void MixedWavefunction::accumulate(int iVib, Wavefunction &wfn, double scaling_constant) {
  auto p = profiler->push("accumulate");
  double *buffer = wfn.buffer.data();
  auto dimension = wfn.dimension;
  int lo, hi, ld = dimension;
  ga_wfn_block_bound(iVib, &lo, &hi, dimension);
  NGA_Acc(m_ga_handle, &lo, &hi, buffer, &ld, &scaling_constant);
}

void MixedWavefunction::operatorOnWavefunction(const MixedOperatorSecondQuant &ham, const MixedWavefunction &w,
                                               bool parallel_stringset, bool with_sync) {
  if (with_sync)
    DivideTasks(1000000000, 1, 1, m_communicator);
  auto prof = profiler->push("MixedWavefunction::operatorOnWavefunction");
  auto res = Wavefunction{m_prototype, 0, m_child_communicator};
  std::unique_ptr<Wavefunction> res2;
  auto ketWfn = Wavefunction{m_prototype, 0, m_child_communicator};
  bool zeroed = false;
  res.allocate_buffer();
  if (!ham.elHam2.empty()) {
    res2 = std::make_unique<Wavefunction>(m_prototype, 0, m_child_communicator);
    res2->allocate_buffer();
  }
  ketWfn.allocate_buffer();
  for (const auto &bra : m_vibBasis) {
    auto iBra = m_vibBasis.index(bra);
    // Purely electronic operators
    if (NextTask(m_communicator) && !ham.elHam.empty()) {
      auto p = profiler->push("Hel");
      copy_to_local(w.m_ga_handle, iBra, ketWfn);
      if (!zeroed) {
        res.zero();
        zeroed = true;
      }
      for (const auto &hel : ham.elHam) {
        auto p = profiler->push(hel.first);
        auto op = hel.second.get();
        res.operatorOnWavefunction(*op, ketWfn, parallel_stringset);
      }
    }
    // Purely electronic operators applied twice
    if (NextTask(m_communicator) && !ham.elHam2.empty()) {
      auto p = profiler->push("Hel2");
      copy_to_local(w.m_ga_handle, iBra, ketWfn);
      if (!zeroed) {
        res.zero();
        zeroed = true;
      }
      for (const auto &hel : ham.elHam2) {
        auto p = profiler->push(hel.first);
        auto op = hel.second.get();
        res2->zero();
        res2->operatorOnWavefunction(*op, ketWfn, parallel_stringset);
        res.operatorOnWavefunction(*op, *res2, parallel_stringset);
      }
    }
    // Purely vibrational operators
    if (NextTask(m_communicator)) {
      auto p = profiler->push("Hvib");
      if (!zeroed) {
        res.zero();
        zeroed = true;
      }
      for (const auto &vibEl : ham.Hvib.tensor) {
        auto val = vibEl.second.oper;
        auto ket = bra.excite(vibEl.second.exc);
        if (!ket.withinSpace(m_vibSpace))
          continue;
        auto iKet = m_vibBasis.index(ket);
        copy_to_local(w.m_ga_handle, iKet, ketWfn);
        res.axpy(val, ketWfn);
      }
    }
    // Mixed operators
    for (const auto &ket : m_vibBasis) {
      auto p = profiler->push("Hmixed");
      auto iKet = m_vibBasis.index(ket);
      if (!ham.connected(bra, ket))
        continue;
      if (!NextTask(m_communicator))
        continue;
      copy_to_local(w.m_ga_handle, iKet, ketWfn);
      if (!zeroed) {
        res.zero();
        zeroed = true;
      }
      for (const auto &mixedTerm : ham.mixedHam) {
        const auto &vibTensor = mixedTerm.second;
        for (const auto &vibEl : vibTensor.tensor) {
          auto &p_op = vibEl.second.oper;
          auto connected_ket = bra.excite(vibEl.second.exc);
          if (connected_ket != ket)
            continue;
          auto op = p_op.get();
          res.operatorOnWavefunction(*op, ketWfn, parallel_stringset);
        }
      }
    }
    if (zeroed) {
      accumulate(iBra, res);
      zeroed = false;
    }
  }
  if (with_sync)
    sync();
}

void MixedWavefunction::diagonalOperator(const MixedOperatorSecondQuant &ham, bool parallel_stringset) {
  auto p = profiler->push("MixedWavefunction::diagonalOperator");
  zero();
  DivideTasks(1000000000, 1, 1, m_communicator);
  auto res = Wavefunction{m_prototype, 0, m_child_communicator};
  auto wfn = Wavefunction{m_prototype, 0, m_child_communicator};
  res.allocate_buffer();
  wfn.allocate_buffer();
  for (const auto &bra : m_vibBasis) {
    auto iBra = m_vibBasis.index(bra);
    // Purely electronic operators
    if (NextTask(m_communicator)) {
      res.zero();
      for (const auto &hel : ham.elHam) {
        auto op = hel.second.get();
        res.diagonalOperator(*op);
      }
      accumulate(iBra, res);
    }
    // Pure vibrational operator
    if (NextTask(m_communicator)) {
      res.zero();
      for (const auto &vibEl : ham.Hvib.tensor) {
        auto val = vibEl.second.oper;
        auto &vibExc = vibEl.second.exc;
        auto ket = bra.excite(vibExc);
        if (ket != bra)
          continue;
        res += val;
      }
      accumulate(iBra, res);
    }
    // all mixed vibrational - electronic operators
    for (const auto &mixedTerm : ham.mixedHam) {
      for (const auto &vibEl : mixedTerm.second.tensor) {
        auto p_op = vibEl.second.oper;
        auto &vibExc = vibEl.second.exc;
        auto ket = bra.excite(vibExc);
        if (ket != bra)
          continue;
        if (!NextTask(m_communicator))
          continue;
        res.zero();
        auto op = p_op.get();
        res.diagonalOperator(*op);
        accumulate(iBra, res);
      }
    }
  }
  sync();
}

bool MixedWavefunction::compatible(const MixedWavefunction &other) const {
  bool sameSize = (m_vibBasis.vibDim() == other.m_vibBasis.vibDim());
  if (!Array::compatible(other))
    return false;
  bool sameVibBasis = (m_vibSpace == other.m_vibSpace);
  bool sameElectronicWfn = m_prototype.compatible(other.m_prototype);
  return sameSize && sameElectronicWfn && sameVibBasis;
}

std::vector<double> MixedWavefunction::vibDensity() {
  auto nM = m_vibSpace.nModal;
  auto size = nM * nM;
  auto dm = std::vector<double>(size, 0.);
  auto bra = Wavefunction{m_prototype, 0, m_child_communicator};
  auto ket = Wavefunction{m_prototype, 0, m_child_communicator};
  bra.allocate_buffer();
  ket.allocate_buffer();
  int local_i = -1;
  int local_j = -1;
  auto norm = dot(this);
  if (std::abs(norm) < 1.0e-8)
    throw std::runtime_error("Norm of wavefunction is too small, " + std::to_string(norm));
  DivideTasks(size, 1, 1, m_communicator);
  for (size_t i = 0; i < nM; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      if (NextTask(m_communicator)) {
        if (local_i != i)
          copy_to_local(m_ga_handle, i, bra);
        local_i = i;
        if (local_j != i)
          copy_to_local(m_ga_handle, j, ket);
        local_j = j;
        dm[i * nM + j] = dm[j * nM + i] = bra.dot(ket) / norm;
      }
    }
  }
  // reduce the whole array together
  MPI_Allreduce(MPI_IN_PLACE, dm.data(), size, MPI_DOUBLE, MPI_SUM, m_communicator);
  return dm;
}

// bool MixedWavefunction::compatible(const Array &other) const {
//    auto other_wfn = dynamic_cast<const MixedWavefunction *>(&other);
//    if (other_wfn == nullptr) return Array::compatible(other);
//    else return compatible(*other_wfn);
//}

} // namespace gci
} // namespace molpro
