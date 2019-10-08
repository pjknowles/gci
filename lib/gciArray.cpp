#include "gciArray.h"

#include "gci.h"

#include <ga.h>
#include <cfloat>
#include <numeric>

namespace gci {

int get_communicator_size(MPI_Comm comm) {
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

int get_communicator_rank(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

Array::Array(MPI_Comm comm) : m_communicator(comm), m_comm_rank(get_communicator_rank(comm)),
                              m_comm_size(get_communicator_size(comm)), m_dimension(0), m_ga_handle(0), m_ga_pgroup(0),
                              m_ga_chunk(1), m_ga_allocated(false) { }

Array::Array(size_t dimension, MPI_Comm comm)
        : m_communicator(comm), m_comm_rank(get_communicator_rank(comm)), m_comm_size(get_communicator_size(comm)),
          m_dimension(dimension), m_ga_handle(0), m_ga_chunk(1), m_ga_pgroup(0), m_ga_allocated(false) { }

Array::Array(const Array &source)
        : m_communicator(source.m_communicator), m_comm_rank(source.m_comm_rank),
          m_comm_size(source.m_comm_size), m_dimension(source.m_dimension), m_ga_handle(0),
          m_ga_chunk(source.m_ga_chunk), m_ga_pgroup(0), m_ga_allocated(false) {
    *this = source;
}

Array &Array::operator=(const Array &source) noexcept {
    m_dimension = source.m_dimension;
    m_communicator = source.m_communicator;
    copy_buffer(source);
    return *this;
}

Array::~Array() {
    if (!empty()) GA_Destroy(m_ga_handle);
}


bool Array::empty() const {
    return !m_ga_allocated;
}

void Array::allocate_buffer() {
    if (!empty()) return;
    if (_ga_pgroups.find(m_communicator) == _ga_pgroups.end()) {
//  global processor ranks
        int loc_size, glob_size, glob_rank;
        MPI_Comm_size(m_communicator, &loc_size);
        MPI_Comm_size(mpi_comm_compute, &glob_size);
        MPI_Comm_rank(mpi_comm_compute, &glob_rank);
        auto glob_ranks = std::vector<int>(loc_size);
        MPI_Allgather(&glob_rank, 1, MPI_INT, glob_ranks.data(), 1, MPI_INT, m_communicator);
// create new GA processor group
        m_ga_pgroup = GA_Pgroup_create(glob_ranks.data(), loc_size);
        _ga_pgroups[m_communicator] = m_ga_pgroup;
    } else m_ga_pgroup = _ga_pgroups[m_communicator];
    m_ga_handle = NGA_Create_handle();
    NGA_Set_pgroup(m_ga_handle, m_ga_pgroup);
    auto dims = (int) m_dimension;
    NGA_Set_data(m_ga_handle, 1, &dims, C_DBL);
    NGA_Set_array_name(m_ga_handle, (char *) "Array");
    NGA_Set_chunk(m_ga_handle, &m_ga_chunk);
    auto succ = GA_Allocate(m_ga_handle);
    if (!succ) GA_Error((char *) "Failed to allocate", 1);
    m_ga_allocated = true;
}

void Array::sync() {
    GA_Pgroup_sync(m_ga_pgroup);
}

void Array::copy_buffer(const Array &source) {
    if (source.empty()) return;
    if (m_communicator != source.m_communicator)
        GA_Error((char *) "trying to copy arrays with different communicators", 1);
    if (size() != source.size()) GA_Error((char *) "trying to copy arrays of different size", 1);
    if (empty()) {
        m_ga_handle = GA_Duplicate(source.m_ga_handle, (char *) "Array");
        m_ga_chunk = source.m_ga_chunk;
        m_ga_allocated = true;
        m_ga_pgroup = source.m_ga_pgroup;
    }
    GA_Copy(source.m_ga_handle, m_ga_handle);
}

/*
 * Each process searches for lowest n numbers in it's local GA buffer, if it has one.
 * Then, sends them to root who does the final sort
 */
std::vector<size_t> Array::minlocN(size_t n) const {
    if (empty()) return {};
    auto iproc = GA_Nodeid();
    auto nproc = GA_Nnodes();
    // get distribution of GA local to this process, and access the buffer
    int lo, hi, ld;
    NGA_Distribution(m_ga_handle, iproc, &lo, &hi);
    double *buffer = nullptr;
    NGA_Access(m_ga_handle, &lo, &hi, &buffer, &ld);
    if (buffer == nullptr) GA_Error((char *) "Failed to access local GA buffer", 1);
    // search for n lowest numbers
    auto length = hi - lo + 1;
    auto nmin = length > n ? n : length;
    std::vector<size_t> minInd(n, 0);
    std::vector<double> minVal(n, DBL_MAX);
    auto ptr_to_min = std::min_element(buffer, buffer + length);
    size_t min_ind = (ptr_to_min - buffer);
    double prev_min = *ptr_to_min;
    minInd[0] = lo + min_ind;
    minVal[0] = prev_min;
    for (int i = 1; i < nmin; ++i) {
        ptr_to_min = std::min_element(buffer, buffer + length,
                                      [prev_min](double el, double smallest) {return el < smallest && el > prev_min || smallest <= prev_min;});
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
    MPI_Barrier(m_communicator);
    if (iproc == 0) {
        MPI_Gather(MPI_IN_PLACE, n, MPI_SIZE_T, minInd.data(), n, MPI_SIZE_T, 0, m_communicator);
        MPI_Gather(MPI_IN_PLACE, n, MPI_DOUBLE, minVal.data(), n, MPI_DOUBLE, 0, m_communicator);
        std::vector<double> sort_permutation(minVal.size());
        std::iota(sort_permutation.begin(), sort_permutation.end(), 0);
        std::sort(sort_permutation.begin(), sort_permutation.end(),
                  [minVal](const auto &i1, const auto &i2) {return minVal[i1] < minVal[i2];});
        std::transform(sort_permutation.begin(), sort_permutation.begin() + n, minInd.begin(),
                       [minInd](auto i) {return minInd[i];});
        minInd.resize(n);
    } else {
        MPI_Gather(minInd.data(), n, MPI_SIZE_T, minInd.data(), n * nproc, MPI_SIZE_T, 0, m_communicator);
        MPI_Gather(minVal.data(), n, MPI_DOUBLE, minVal.data(), n * nproc, MPI_DOUBLE, 0, m_communicator);
    }
    MPI_Bcast(minInd.data(), n, MPI_SIZE_T, 0, m_communicator);
    return minInd;
}


double Array::at(size_t ind) const {
    if (ind >= m_dimension) GA_Error((char *) "Out of bounds", 1);
    if (empty()) GA_Error((char *) "GA was not allocated", 1);
    double buffer;
    int lo = ind, high = ind, ld = 1;
    NGA_Get(m_ga_handle, &lo, &high, &buffer, &ld);
    return buffer;
}

std::vector<double> Array::vec() const {
    if (empty()) return {};
    std::vector<double> vec(m_dimension);
    double *buffer = vec.data();
    int lo = 0, hi = m_dimension - 1, ld = m_dimension;
    NGA_Get(m_ga_handle, &lo, &hi, buffer, &ld);
    return vec;
}

std::vector<double> Array::get(int lo, int hi) {
    if (empty()) return {};
    int ld, n = hi - lo + 1;
    auto data = std::vector<double>(n);
    NGA_Get(m_ga_handle, &lo, &hi, data.data(), &ld);
    return data;
}

void Array::put(int lo, int hi, double *data) {
    if (empty()) GA_Error((char *) "Attempting to put data into an empty GA", 1);
    int ld;
    NGA_Put(m_ga_handle, &lo, &hi, data, &ld);
}

std::vector<double> Array::gather(std::vector<int> &indices) const {
    int n = indices.size();
    std::vector<double> data(n, 0.);
    int **subsarray = new int *[n];
    for (int i = 0; i < n; ++i) {
        subsarray[i] = &(indices.at(i));
    }
    NGA_Gather(m_ga_handle, data.data(), subsarray, n);
    delete[] subsarray;
    return data;
}

void Array::scatter(std::vector<int> &indices, std::vector<double> &vals) {
    int n = indices.size();
    int **subsarray = new int *[n];
    for (int i = 0; i < n; ++i) {
        subsarray[i] = &(indices.at(i));
    }
    NGA_Scatter(m_ga_handle, vals.data(), subsarray, n);
    delete[] subsarray;
}

void Array::scatter_acc(std::vector<int> &indices, std::vector<double> &vals, double alpha) {
    int n = indices.size();
    int **subsarray = new int *[n];
    for (int i = 0; i < n; ++i) {
        subsarray[i] = &(indices.at(i));
    }
    NGA_Scatter_acc(m_ga_handle, vals.data(), subsarray, n, &alpha);
    delete[] subsarray;
}

bool Array::compatible(const Array &other) const {
    return m_dimension == other.m_dimension;
}

void Array::axpy(double a, const Array &x) {
    if (!compatible(x)) GA_Error((char *) "Attempting to add incompatible Array objects", 1);
    if (empty() || x.empty()) GA_Error((char *) "GA not allocated", 1);
    double one = 1.0;
    GA_Add(&one, m_ga_handle, &a, x.m_ga_handle, m_ga_handle);
}

void Array::axpy(double a, const std::map<size_t, double> &x) {
    // TODO do I really need to cast away the const to get values from a map???
    int n = x.size();
    std::vector<int> indices;
    std::vector<double> vals;
    indices.reserve(n);
    vals.reserve(n);
    for (const auto &xx : x) {
        indices.push_back(xx.first);
        vals.push_back(xx.second);
    }
    scatter_acc(indices, vals, a);
}

void Array::scal(double a) {
    if (empty()) GA_Error((char *) "GA not allocated", 1);
    GA_Scale(m_ga_handle, &a);
}

void Array::add(const Array &other) {
    axpy(1.0, other);
}

void Array::add(double a) {
    GA_Add_constant(m_ga_handle, &a);
}

void Array::sub(const Array &other) {
    axpy(-1.0, other);
}

void Array::sub(double a) {
    add(-a);
}

void Array::recip() {
    if (empty()) GA_Error((char *) "GA not allocated", 1);
    GA_Recip(m_ga_handle);
}

double Array::dot(const Array &other) const {
    if (!compatible(other)) GA_Error((char *) "Attempt to form scalar product between incompatible Array objects", 1);
    if (empty()) GA_Error((char *) "GA not allocated", 1);
    return GA_Ddot(m_ga_handle, other.m_ga_handle);
}

double Array::dot(const std::map<size_t, double> &other) const {
    double result = 0;
    if (m_comm_rank == 0) {
        auto n = other.size();
        std::vector<int> indices;
        std::vector<double> other_vals;
        indices.reserve(n);
        other_vals.reserve(n);
        for (const auto &item : other) {
            indices.push_back(item.first);
            other_vals.push_back(item.second);
        }
        auto this_vals = gather(indices);
        result = std::inner_product(this_vals.begin(), this_vals.end(), other_vals.begin(), 0.);
    }
    MPI_Bcast(&result, 1, MPI_DOUBLE, 0, m_communicator);
    return result;
}

void Array::zero() {
    if (empty()) allocate_buffer();
    GA_Zero(m_ga_handle);
}

void Array::set(double val) {
    if (empty()) GA_Error((char *) "GA not allocated", 1);
    GA_Fill(m_ga_handle, &val);
}

void Array::set(size_t ind, double val) {
    auto data = std::vector<double>{val};
    put(ind, ind, data.data());
}

Array &Array::operator*=(double value) {
    scal(value);
    return *this;
}

Array &Array::operator+=(const Array &other) {
    add(other);
    return *this;
}

Array &Array::operator-=(const Array &other) {
    sub(other);
    return *this;
}

Array &Array::operator+=(double value) {
    add(value);
    return *this;
}

Array &Array::operator-=(double value) {
    sub(value);
    return *this;
}

Array &Array::operator-() {
    scal(-1.0);
    return *this;
}

Array &Array::operator/=(const Array &other) {
    if (!compatible(other)) GA_Error((char *) "Attempting to divide incompatible Array objects", 1);
    if (empty() || other.empty()) GA_Error((char *) "GA not allocated", 1);
    GA_Elem_divide(m_ga_handle, other.m_ga_handle, m_ga_handle);
    return *this;
}

void Array::times(const Array *a, const Array *b) {
    if (a == nullptr || b == nullptr) GA_Error((char *) "Vectors cannot be null", 1);
    if (!compatible(*a) || !compatible(*b)) GA_Error((char *) "Attempting to divide incompatible Array objects", 1);
    if (empty() || a->empty() || b->empty()) GA_Error((char *) "GA not allocated", 1);
    GA_Elem_multiply(a->m_ga_handle, b->m_ga_handle, m_ga_handle);
}

// this[i] = a[i]/(b[i]+shift)
void Array::divide(const Array *a, const Array *b,
                   double shift, bool append, bool negative) {
    if (a == nullptr || b == nullptr) GA_Error((char *) "Vectors cannot be null", 1);
    if (!compatible(*a) || !compatible(*b)) GA_Error((char *) "Attempting to divide incompatible Array objects", 1);
    if (empty() || a->empty() || b->empty()) GA_Error((char *) "GA not allocated", 1);
    shift = negative ? -shift : shift;
    auto b_shifted = Array(*b);
    b_shifted.add(shift);
    if (append) {
        GA_Elem_divide(a->m_ga_handle, b_shifted.m_ga_handle, b_shifted.m_ga_handle);
        add(b_shifted);
    } else
        GA_Elem_divide(a->m_ga_handle, b_shifted.m_ga_handle, m_ga_handle);
}


double operator*(const Array &w1, const Array &w2) {
    return w1.dot(w2);
}

Array operator+(const Array &w1, const Array &w2) {
    Array result = w1;
    return result += w2;
}

Array operator-(const Array &w1, const Array &w2) {
    Array result = w1;
    return result -= w2;
}

Array operator/(const Array &w1, const Array &w2) {
    Array result = w1;
    return result /= w2;
}

Array operator*(const Array &w1, const double &value) {
    Array result = w1;
    return result *= value;
}

Array operator*(const double &value, const Array &w1) {
    Array result = w1;
    return result *= value;
}

} // namespace gci
