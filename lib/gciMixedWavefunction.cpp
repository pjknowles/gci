#include "gciMixedWavefunction.h"
#include "gciHProductSet.h"

#include <ga.h>
#include <mpi.h>
#include <stdexcept>

namespace gci {

MPI_Comm create_new_comm(MPI_Comm head_comm = MPI_COMM_COMPUTE) {
    MPI_Comm new_comm = MPI_COMM_NULL;
    MPI_Comm_split(head_comm, GA_Nodeid(), GA_Nodeid(), &new_comm);
    if (new_comm == MPI_COMM_NULL) GA_Error((char *) "Failed to create a new communicator", 0);
    return new_comm;
}

MixedWavefunction::MixedWavefunction(const Options &options, const State &prototype, MPI_Comm head_commun)
        : m_head_communicator(head_commun), m_child_communicator(create_new_comm(m_head_communicator)),
          m_vibSpace(options.parameter("NMODE", 0), options.parameter("NMODAL", 1),
                     options.parameter("VIB_EXC_LVL", 1)), m_vibBasis(m_vibSpace), m_elDim(0), m_dimension(0),
          m_prototype(prototype), m_ga_handle(-1), m_ga_chunk(1), m_ga_allocated(false) {
    m_elDim = m_prototype.size();
    m_vibBasis.generateFullSpace();
    m_dimension = m_vibBasis.vibDim() * m_elDim;
}

MixedWavefunction::MixedWavefunction(const MixedWavefunction &source, int option, MPI_Comm head_commun)
        : m_head_communicator(head_commun), m_child_communicator(create_new_comm(m_head_communicator)),
          m_vibSpace(source.m_vibSpace), m_vibBasis(source.m_vibBasis), m_elDim(source.m_elDim),
          m_dimension(source.m_dimension), m_prototype(source.m_prototype), m_ga_handle(-1),
          m_ga_chunk(source.m_ga_chunk), m_ga_allocated(false) {
    copy_buffer(source);
    allocate_buffer();
}

MixedWavefunction::~MixedWavefunction() {
    if (!empty()) GA_Destroy(m_ga_handle);
    MPI_Comm_free(&m_child_communicator);
}


bool MixedWavefunction::empty() const {
    return !m_ga_allocated;
}

void MixedWavefunction::allocate_buffer() {
    if (_ga_pgroups.find(m_head_communicator) == _ga_pgroups.end()){
//  global processor ranks
        int loc_size, glob_size, glob_rank;
        MPI_Comm_size(m_head_communicator, &loc_size);
        MPI_Comm_size(MPI_COMM_COMPUTE, &glob_size);
        MPI_Comm_rank(MPI_COMM_COMPUTE, &glob_rank);
        std::vector<int> glob_ranks{loc_size};
        MPI_Allgather(&glob_rank, 1, MPI_INT, &glob_ranks[0], loc_size, MPI_INT, m_head_communicator);
// create new GA processor group
        m_ga_pgroup = GA_Pgroup_create(&glob_ranks[0], loc_size);
        _ga_pgroups[m_head_communicator] = m_ga_pgroup;
    } else m_ga_pgroup = _ga_pgroups[m_head_communicator];
    m_ga_handle = NGA_Create_handle();
    NGA_Set_pgroup(m_ga_handle, m_ga_pgroup);
    auto dims = (int) m_dimension;
    NGA_Set_data(m_ga_handle, 1, &dims, C_DBL);
    NGA_Set_array_name(m_ga_handle, (char *) "MixedWavefunction");
    NGA_Set_chunk(m_ga_handle, &m_ga_chunk);
    auto succ = GA_Allocate(m_ga_handle);
    if (!succ) GA_Error((char *) "Failed to allocate", 0);
    m_ga_allocated = true;
}

//TODO does it still work if source has a different communicator?
void MixedWavefunction::copy_buffer(const MixedWavefunction &source) {
    if (source.empty()) return;
    if (empty()) {
        char name[] = "MixedWavefunction";
        m_ga_handle = GA_Duplicate(source.m_ga_handle, name);
        m_ga_chunk = source.m_ga_chunk;
        m_ga_allocated = true;
    } else GA_Copy(source.m_ga_handle, m_ga_handle);
}

/*
 * Each process searches for lowest n numbers in it's local GA buffer, if it has one.
 * Then, sends them to root who does the final sort
 */
std::vector<size_t> MixedWavefunction::minlocN(size_t n) const {
    if (empty()) return {};
    auto iproc = GA_Nodeid();
    auto nproc = GA_Nnodes();
    // get distribution of GA local to this process, and access the buffer
    int lo, hi;
    NGA_Distribution(m_ga_handle, iproc, &lo, &hi);
    auto length = hi - lo;
    double *buffer = nullptr;
    int ld;
    std::array<double,10> arr;
    for (int i =0; i < 10; ++i){
        arr[i] = at(i);
    }
    std::cout << arr[0] << std::endl;
    NGA_Access(m_ga_handle, &lo, &hi, &buffer, &ld);
    if (buffer == nullptr) GA_Error((char *) "Failed to access local GA buffer", 0);
    // search for n lowest numbers
    auto nmin = length > n ? n : length;
    std::vector<size_t> minInd(n, 0);
    std::vector<value_type> minVal(n, DBL_MAX);
    auto ptr_to_min = std::min_element(buffer, buffer + length);
    size_t min_ind = (ptr_to_min - buffer);
    value_type prev_min = *ptr_to_min;
    minInd[0] = lo + min_ind;
    minVal[0] = prev_min;
    for (int i = 1; i < nmin; ++i) {
        ptr_to_min = std::min_element(buffer, buffer + length,
                                      [prev_min](double el1, double el2) {return el1 < el2 && el1 > prev_min;});
        min_ind = (ptr_to_min - buffer);
        prev_min = *ptr_to_min;
        minInd[i] = lo + min_ind;
        minVal[i] = prev_min;
    }
    NGA_Release(m_ga_handle, &lo, &hi);
    // root collects values, does the final sort and sends the result back
    if (iproc == 0) {
        minInd.resize(n * nproc);
        minVal.resize(n * nproc);
    }
    if (iproc == 0) {
        MPI_Gather(MPI_IN_PLACE, n, MPI_SIZE_T, minInd.data(), n * nproc, MPI_SIZE_T, 0, m_head_communicator);
        MPI_Gather(MPI_IN_PLACE, n, MPI_DOUBLE, minVal.data(), n * nproc, MPI_DOUBLE, 0, m_head_communicator);
        std::vector<value_type> sort_permutation(minVal.size());
        std::iota(sort_permutation.begin(), sort_permutation.end(), 0);
        std::sort(sort_permutation.begin(), sort_permutation.end(),
                  [minVal](const auto &i1, const auto &i2) {return minVal[i1] < minVal[i2];});
        std::transform(sort_permutation.begin(), sort_permutation.begin() + n, minInd.begin(),
                       [minInd](auto i) {return minInd[i];});
        minInd.resize(n);
    }
    else{
        MPI_Gather(minInd.data(), n, MPI_SIZE_T, minInd.data(), n * nproc, MPI_SIZE_T, 0, m_head_communicator);
        MPI_Gather(minVal.data(), n, MPI_DOUBLE, minVal.data(), n * nproc, MPI_DOUBLE, 0, m_head_communicator);
    }
    MPI_Bcast(minInd.data(), n, MPI_SIZE_T, 0, m_head_communicator);
    return minInd;
}


double MixedWavefunction::at(size_t ind) const {
    if (ind >= m_dimension) throw std::logic_error("Out of bounds");
    if (empty()) throw std::logic_error("GA was not allocated");
    value_type buffer;
    int lo = ind, high = ind, ld = 1;
    NGA_Get(m_ga_handle, &lo, &high, &buffer, &ld);
    return buffer;
}

Wavefunction MixedWavefunction::wavefunctionAt(size_t iVib) const {
    auto wfn = Wavefunction{m_prototype};
    wfn.allocate_buffer();
    ga_copy_to_local(m_ga_handle, iVib, wfn, m_elDim);
    return wfn;
}

std::string MixedWavefunction::str() const {
    std::string out = "MixedWavefunction elements: ";
    for (size_t i = 0; i < m_dimension; ++i) {
        out += std::to_string(at(i)) + ", ";
    }
    return out;
}

void MixedWavefunction::ga_wfn_block_bound(int iVib, int *lo, int *hi, int dimension) {
    lo[0] = iVib * dimension;
    hi[0] = lo[0] + dimension - 1;
}

void MixedWavefunction::ga_copy_to_local(int ga_handle, int iVib, Wavefunction &wfn, int dimension) {
    auto p = profiler->push("ga_copy_to_local");
    value_type *buffer = wfn.buffer.data();
    int lo, hi, ld = dimension;
    ga_wfn_block_bound(iVib, &lo, &hi, dimension);
    NGA_Get(ga_handle, &lo, &hi, buffer, &ld);
}

void MixedWavefunction::ga_accumulate(int ga_handle, int iVib, Wavefunction &wfn, int dimension, value_type scaling_constant) {
    auto p = profiler->push("ga_accumulate");
    value_type *buffer = wfn.buffer.data();
    int lo, hi, ld = dimension;
    ga_wfn_block_bound(iVib, &lo, &hi, dimension);
    NGA_Acc(ga_handle, &lo, &hi, buffer, &ld, &scaling_constant);
}

void MixedWavefunction::operatorOnWavefunction(const MixedOperatorSecondQuant &ham, const MixedWavefunction &w,
                                               bool parallel_stringset) {
    DivideTasks(1000000000, 1, 1, m_head_communicator);
    auto prof = profiler->push("MixedWavefunction::operatorOnWavefunction");
    auto res = Wavefunction{m_prototype};
    auto ketWfn = Wavefunction{m_prototype};
    res.allocate_buffer();
    ketWfn.allocate_buffer();
    for (const auto &bra : m_vibBasis) {
        auto iBra = m_vibBasis.index(bra);
        // Purely electronic operators
        if (NextTask(m_head_communicator)) {
            auto p = profiler->push("Hel");
            ga_copy_to_local(m_ga_handle, iBra, ketWfn, m_elDim);
            res.zero();
            for (const auto &hel : ham.elHam) {
                auto p = profiler->push(hel.first);
                res.operatorOnWavefunction(*hel.second, ketWfn, parallel_stringset);
            }
            ga_accumulate(m_ga_handle, iBra, res, m_elDim);
        }
        // Purely vibrational operators
        if (NextTask(m_head_communicator)) {
            auto p = profiler->push("Hvib");
            res.zero();
            for (const auto &vibEl : ham.Hvib.tensor) {
                auto val = vibEl.second.oper;
                auto vibExc = VibExcitation{vibEl.second.exc};
                vibExc.conjugate();
                auto ket = bra.excite(vibExc);
                if (!ket.withinSpace(m_vibSpace)) continue;
                auto iKet = m_vibBasis.index(ket);
                ga_copy_to_local(m_ga_handle, iKet, ketWfn, m_elDim);
                res.axpy(val, ketWfn);
            }
            ga_accumulate(m_ga_handle, iBra, res, m_elDim);
        }
        // Mixed operators
        for (const auto &ket : m_vibBasis) {
            auto p = profiler->push("Hmixed");
            auto iKet = m_vibBasis.index(bra);
            if (!ham.connected(bra, ket)) continue;
            if (!NextTask(m_head_communicator)) continue;
            ga_copy_to_local(m_ga_handle, iKet, ketWfn, m_elDim);
            res.zero();
            for (const auto &mixedTerm : ham.mixedHam) {
                const auto &vibTensor = mixedTerm.second;
                for (const auto &vibEl : vibTensor.tensor) {
                    auto &op = vibEl.second.oper;
                    auto vibExc = VibExcitation{vibEl.second.exc};
                    vibExc.conjugate();
                    auto connected_ket = bra.excite(vibExc);
                    if (connected_ket != ket) continue;
                    res.operatorOnWavefunction(op, ketWfn, parallel_stringset);
                }
            }
            ga_accumulate(m_ga_handle, iBra, res, m_elDim);
        }
    }
}

void MixedWavefunction::diagonalOperator(const MixedOperatorSecondQuant &ham, bool parallel_stringset) {
    auto p = profiler->push("MixedWavefunction::diagonalOperator");
    DivideTasks(1000000000, 1, 1, m_head_communicator);
    auto res = Wavefunction{m_prototype};
    auto wfn = Wavefunction{m_prototype};
    res.allocate_buffer();
    wfn.allocate_buffer();
    for (const auto &bra : m_vibBasis) {
        auto iBra = m_vibBasis.index(bra);
        // Purely electronic operators
        if (NextTask(m_head_communicator)) {
            res.zero();
            for (const auto &hel : ham.elHam) {
                res.diagonalOperator(*hel.second);
            }
            ga_accumulate(m_ga_handle, iBra, res, m_elDim);
        }
        // Pure vibrational operator
        if (NextTask(m_head_communicator)) {
            res.zero();
            for (const auto &vibEl : ham.Hvib.tensor) {
                auto val = vibEl.second.oper;
                auto &vibExc = vibEl.second.exc;
                auto ket = bra.excite(vibExc);
                if (ket != bra) continue;
                res += val;
            }
            ga_accumulate(m_ga_handle, iBra, res, m_elDim);
        }
        // all mixed vibrational - electronic operators
        for (const auto &mixedTerm : ham.mixedHam) {
            for (const auto &vibEl : mixedTerm.second.tensor) {
                auto op = vibEl.second.oper;
                auto &vibExc = vibEl.second.exc;
                auto ket = bra.excite(vibExc);
                if (ket != bra) continue;
                if (!NextTask(m_head_communicator)) continue;
                res.zero();
                res.diagonalOperator(op);
                ga_accumulate(m_ga_handle, iBra, res, m_elDim);
            }
        }
    }
}

std::vector<double> MixedWavefunction::vec() const {
    if (!empty()) return {};
    std::vector<double> vec(m_dimension);
    value_type *buffer = vec.data();
    int lo=0, hi=m_dimension-1, ld = m_dimension;
    NGA_Get(m_ga_handle, &lo, &hi, buffer, &ld);
    return vec;
}

bool MixedWavefunction::compatible(const MixedWavefunction &other) const {
    bool sameSize = (m_vibBasis.vibDim() == other.m_vibBasis.vibDim());
    if (!sameSize) return sameSize;
    bool sameVibBasis = (m_vibSpace == other.m_vibSpace);
    bool sameElectronicWfn = m_prototype.compatible(other.m_prototype);
    return sameSize && sameElectronicWfn && sameVibBasis;
}

void MixedWavefunction::axpy(double a, const MixedWavefunction &x) {
    if (!compatible(x)) throw std::domain_error("Attempting to add incompatible MixedWavefunction objects");
    if (empty() || !x.empty()) throw std::runtime_error("GA not allocated");
    double b = 1.0;
    GA_Add(&a, m_ga_handle, &b, x.m_ga_handle, m_ga_handle);
}

void MixedWavefunction::scal(double a) {
    if (empty()) throw std::runtime_error("GA not allocated");
    GA_Scale(m_ga_handle, &a);
}

void MixedWavefunction::add(const MixedWavefunction &other) {
    axpy(1.0, other);
}

void MixedWavefunction::add(double a) {
    GA_Add_constant(m_ga_handle, &a);
}

void MixedWavefunction::sub(const MixedWavefunction &other) {
    axpy(-1.0, other);
}

void MixedWavefunction::sub(double a) {
    a *= -1.0;
    add(a);
}

void MixedWavefunction::recip() {
    if (empty()) throw std::runtime_error("GA not allocated");
    GA_Recip(m_ga_handle);
}

double MixedWavefunction::dot(const MixedWavefunction &other) const {
    if (!compatible(other))
        throw std::domain_error("Attempt to form scalar product between incompatible MixedWavefunction objects");
    if (empty()) throw std::runtime_error("GA not allocated");
    return GA_Ddot(m_ga_handle, other.m_ga_handle);
}

void MixedWavefunction::zero() {
    if (empty()) allocate_buffer();
    GA_Zero(m_ga_handle);
}

void MixedWavefunction::set(double val) {
    GA_Fill(m_ga_handle, &val);
}

void MixedWavefunction::set(size_t ind, double val) {
    if (ind >= m_dimension) throw std::logic_error("Out of bounds");
    if (empty()) throw std::logic_error("GA not allocated");
    int i = ind;
    int * subs = &i;
    NGA_Scatter(m_ga_handle, &val, &subs, 1);
}

MixedWavefunction &MixedWavefunction::operator*=(const double &value) {
    scal(value);
    return *this;
}

MixedWavefunction &MixedWavefunction::operator+=(const MixedWavefunction &other) {
    add(other);
    return *this;
}

MixedWavefunction &MixedWavefunction::operator-=(const MixedWavefunction &other) {
    sub(other);
    return *this;
}

MixedWavefunction &MixedWavefunction::operator+=(double value) {
    add(value);
    return *this;
}

MixedWavefunction &MixedWavefunction::operator-=(double value) {
    sub(value);
    return *this;
}

MixedWavefunction &MixedWavefunction::operator-() {
    scal(-1.0);
    return *this;
}

MixedWavefunction &MixedWavefunction::operator/=(const MixedWavefunction &other) {
    if (!compatible(other)) throw std::domain_error("Attempting to divide incompatible MixedWavefunction objects");
    if (empty() || !other.empty()) throw std::runtime_error("GA not allocated");
    GA_Elem_divide(m_ga_handle, other.m_ga_handle, m_ga_handle);
    return *this;
}

void MixedWavefunction::times(const MixedWavefunction *a, const MixedWavefunction *b) {
    if (a == nullptr || b == nullptr) throw std::logic_error("Vectors cannot be null");
    if (!compatible(*a) || !compatible(*b))
        throw std::domain_error("Attempting to divide incompatible MixedWavefunction objects");
    if (empty() || !a->empty() || b->empty()) throw std::runtime_error("GA not allocated");
    GA_Elem_multiply(a->m_ga_handle, b->m_ga_handle, m_ga_handle);
}

// this[i] = a[i]/(b[i]+shift)
void MixedWavefunction::divide(const MixedWavefunction *a, const MixedWavefunction *b,
                               double shift, bool append, bool negative) {
    if (a == nullptr || b == nullptr) throw std::logic_error("Vectors cannot be null");
    if (!compatible(*a) || !compatible(*b))
        throw std::domain_error("Attempting to divide incompatible MixedWavefunction objects");
    if (empty() || !a->empty() || b->empty()) throw std::runtime_error("GA not allocated");
    shift = negative ? -shift : shift;
    auto b_shifted = MixedWavefunction(*b);
    b_shifted.add(shift);
    if (!append)
        GA_Elem_divide(a->m_ga_handle, b_shifted.m_ga_handle, m_ga_handle);
    else {
        GA_Elem_divide(a->m_ga_handle, b_shifted.m_ga_handle, b_shifted.m_ga_handle);
        add(b_shifted);
    }
}


double operator*(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    return w1.dot(w2);
}

MixedWavefunction operator+(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    MixedWavefunction result = w1;
    return result += w2;
}

MixedWavefunction operator-(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    MixedWavefunction result = w1;
    return result -= w2;
}

MixedWavefunction operator/(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    MixedWavefunction result = w1;
    return result /= w2;
}

MixedWavefunction operator*(const MixedWavefunction &w1, const double &value) {
    MixedWavefunction result = w1;
    return result *= value;
}

MixedWavefunction operator*(const double &value, const MixedWavefunction &w1) {
    MixedWavefunction result = w1;
    return result *= value;
}

}  // namespace gci
