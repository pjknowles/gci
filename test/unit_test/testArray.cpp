#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <algorithm>

#include <gci.h>
#include <gciArray.h>
#include <numeric>
#include <ga.h>

using ::testing::Pointwise;
using ::testing::DoubleEq;
using ::testing::ContainerEq;
using ::testing::Each;

namespace gci {
int n_lock_calls = 0;

class Lock {
public:
    explicit Lock(int mutex = 0) : mutex(mutex) {
        GA_Lock(mutex);
        ++n_lock_calls;
    }

    ~Lock() {GA_Unlock(mutex);}

    int mutex;
};

class ArrayWithMutexF : public ::testing::Test {
public:
    ArrayWithMutexF() {sync();}

};

TEST_F(ArrayWithMutexF, constructor_empty) {
    auto l = Lock();
    auto a = Array();
}

TEST_F(ArrayWithMutexF, constructor_with_commun) {
    auto l = Lock();
    auto a = Array(mpi_comm_compute);
}

TEST_F(ArrayWithMutexF, constructor_with_dimension) {
    auto l = Lock();
    int dim = 100;
    auto a = Array(dim, mpi_comm_compute);
}

TEST_F(ArrayWithMutexF, constructor_copy) {
    auto l = Lock();
    int dim = 100;
    auto a = Array(dim, mpi_comm_compute);
    auto b = Array(a);
}

class ArrayInitializationF : public ::testing::Test, public Array {
public:
    ArrayInitializationF() : Array(dim, mpi_comm_compute) {sync();};
    static const int dim = 100;
};

TEST_F(ArrayInitializationF, size) {
    auto l = Lock();
    ASSERT_EQ(m_dimension, size());
}

TEST_F(ArrayInitializationF, empty) {
    auto l = Lock();
    EXPECT_FALSE(m_ga_allocated);
    ASSERT_TRUE(empty());
}

TEST_F(ArrayInitializationF, allocate_buffer) {
    allocate_buffer();
    auto l = Lock();
    EXPECT_TRUE(m_ga_allocated);
    ASSERT_FALSE(empty());
}

TEST_F(ArrayInitializationF, zero) {
    zero(false, true);
    auto l = Lock();
    EXPECT_TRUE(m_ga_allocated);
    ASSERT_FALSE(empty());
}

TEST_F(ArrayInitializationF, vec) {
    zero(true, true);
    auto data = vec();
    auto empty_vec = std::vector<double>(data.size(), 0.);
    auto l = Lock();
    ASSERT_THAT(empty_vec, Pointwise(DoubleEq(), empty_vec));
}

TEST_F(ArrayInitializationF, get) {
    zero(true, true);
    auto data = get(0, 10);
    auto ref_vec = std::vector<double>(data.size(), 0.);
    {
        auto l = Lock();
        ASSERT_THAT(ref_vec, Pointwise(DoubleEq(), data));
    }
    data = get(0, dim - 1);
    auto same_as_vec = vec();
    auto l = Lock();
    ASSERT_THAT(data, Pointwise(DoubleEq(), same_as_vec));
}

TEST_F(ArrayInitializationF, put) {
    allocate_buffer();
    auto range = std::vector<double>(dim);
    std::iota(range.begin(), range.end(), 0);
    if (m_comm_rank == 0)
        put(0, dim - 1, range.data());
    sync();
    auto from_ga_buffer = get(0, dim - 1);
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), range));
}

TEST_F(ArrayInitializationF, set) {
    zero(true, true);
    auto ref_values = std::vector<double>(dim, 0.);
    auto from_ga_buffer = vec();
    {
        auto l = Lock();
        ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
    }
    set(42.0, true, true);
    std::fill(ref_values.begin(), ref_values.end(), 42.0);
    from_ga_buffer = vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
}


class ArrayRangeF : public ::testing::Test, public Array {
public:
    //! Stores a range in the buffer {0, 1, 2, 3, 4, 5, .., dim-1}
    ArrayRangeF() : Array((size_t) dim, mpi_comm_compute), p_rank(GA_Nodeid()), p_size(GA_Nnodes()) {
        allocate_buffer();
        values.resize(dim);
        std::iota(values.begin(), values.end(), 1.);
        if (p_rank == 0)
            put(0, (int) values.size() - 1, values.data());
        sub_indices.reserve(sub_dim);
        sub_values.reserve(sub_dim);
        for (auto el : {0, 1, 11, 31, 40, 99}) {
            sub_indices.push_back(el);
            sub_values.push_back(values[el]);
        }
        sync();
    }


    static const int dim = 100;
    static const int sub_dim = 6;
    int p_rank, p_size;
    std::vector<double> values;
    std::vector<int> sub_indices;
    std::vector<double> sub_values;
};

TEST_F(ArrayRangeF, gather) {
    auto from_ga_buffer = gather(sub_indices);
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), sub_values));
}

TEST_F(ArrayRangeF, scatter) {
    zero(true, true);
    scatter(sub_indices, sub_values);
    auto from_ga_buffer = gather(sub_indices);
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), sub_values));
}

TEST_F(ArrayRangeF, scatter_acc) {
    auto l = Lock();
    scatter_acc(sub_indices, sub_values, 1.0);
    auto from_ga_buffer = gather(sub_indices);
    auto ref_values = sub_values;
    for (auto &el : ref_values) {
        el *= 2;
    }
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), ref_values));
    scatter_acc(sub_indices, sub_values, -1.0);
}

TEST_F(ArrayRangeF, at) {
    auto from_ga_buffer = std::vector<double>();
    for (int i = 0; i < dim; ++i) {
        from_ga_buffer.push_back(this->at(i));
    }
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
}

TEST_F(ArrayRangeF, minlocN) {
    int n = 10;
    auto ref_minloc_ind = std::vector<size_t>(n);
    std::iota(ref_minloc_ind.begin(), ref_minloc_ind.end(), 0);
    auto minloc_ind = minlocN(n);
    auto l = Lock();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
}

TEST_F(ArrayRangeF, minlocN_reverse) {
    int n = 10;
    std::reverse(values.begin(), values.end());
    if (p_rank == 0)
        put(0, dim - 1, values.data());
    sync();
    auto ref_minloc_ind = std::vector<size_t>(n);
    std::iota(ref_minloc_ind.rbegin(), ref_minloc_ind.rend(), dim - n);
    auto minloc_ind = minlocN(n);
    auto l = Lock();
    ASSERT_THAT(minloc_ind, ContainerEq(ref_minloc_ind));
}

TEST_F(ArrayRangeF, scal_double) {
    double alpha = 1.5;
    scal(alpha, false, true);
    for (auto &el : values) {
        el *= alpha;
    }
    auto from_ga_buffer = vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
}

TEST_F(ArrayRangeF, add_double) {
    double alpha = 100.;
    add(alpha, false, true);
    for (auto &el : values) {
        el += alpha;
    }
    auto from_ga_buffer = vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
}

TEST_F(ArrayRangeF, sub_double) {
    double alpha = 100.;
    sub(alpha, false, true);
    for (auto &el : values) {
        el -= alpha;
    }
    auto from_ga_buffer = vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
}

TEST_F(ArrayRangeF, recip) {
    recip(false, true);
    for (auto &el : values) {
        el = 1. / el;
    }
    auto from_ga_buffer = vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer, Pointwise(DoubleEq(), values));
}

class ArrayCollectiveOpF : public ::testing::Test {
public:
    ArrayCollectiveOpF() : alpha(1.0), beta(2.0), p_rank(GA_Nodeid()), p_size(GA_Nnodes()), a(dim, mpi_comm_compute),
                           b(dim, mpi_comm_compute) {
        a.allocate_buffer();
        b.allocate_buffer();
        range_alpha.resize(dim);
        range_beta.resize(dim);
        std::iota(range_alpha.begin(), range_alpha.end(), 0);
        for (const auto &el : range_alpha) {
            range_beta.push_back(el * beta);
        }
        auto n = dim / every_nth_element;
        auto indices = std::vector<size_t>(n);
        for (int i = 0; i < n; ++i) {
            indices[i] = range_beta[i * every_nth_element];
        }
        for (const auto &el : indices) {
            sparse_array.emplace(el, range_beta[el]);
        }
        sync();
    }

    static const int dim = 100;
    const double alpha;
    const double beta;
    std::vector<double> range_alpha;
    std::vector<double> range_beta;
    int p_rank, p_size;
    Array a;
    Array b;
    static const int every_nth_element = 5; // period for selection of sparse array elements
    std::map<size_t, double> sparse_array; // sparse selection of range_beta elements

};

TEST_F(ArrayCollectiveOpF, add) {
    a.set(alpha);
    b.set(beta);
    a.add(b);
    sync();
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha + beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
}

TEST_F(ArrayCollectiveOpF, sub) {
    a.set(alpha);
    b.set(beta);
    a.sub(b);
    sync();
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha - beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
}

TEST_F(ArrayCollectiveOpF, axpy_Reference) {
    double scale = -3.0;
    a.set(alpha);
    b.set(beta);
    a.axpy(scale, b);
    sync();
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha + scale * beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
}

TEST_F(ArrayCollectiveOpF, axpy_Pointer) {
    double scale = -3.0;
    a.set(alpha);
    b.set(beta);
    a.axpy(scale, &b, false, true);
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Each(DoubleEq(alpha + scale * beta)));
    ASSERT_THAT(from_ga_buffer_b, Each(DoubleEq(beta)));
}

TEST_F(ArrayCollectiveOpF, axpy_map) {
    a.set(alpha);
    double scale = 5.0;
    a.axpy(scale, sparse_array, false, true);
    auto from_ga_buffer_a = a.vec();
    auto ref_vals = std::vector<double>(dim, alpha);
    for (const auto &item : sparse_array) {
        ref_vals[item.first] += scale * item.second;
    }
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_vals));
}

TEST_F(ArrayCollectiveOpF, dot_Array) {
    if (p_rank == 0) {
        a.put(0, dim - 1, range_alpha.data());
        b.put(0, dim - 1, range_beta.data());
    }
    auto ga_dot = a.dot(b, true);
    double ref_dot = std::inner_product(range_alpha.begin(), range_alpha.end(), range_beta.begin(), 0.);
    auto l = Lock();
    ASSERT_THAT(ga_dot, DoubleEq(ref_dot));
}

TEST_F(ArrayCollectiveOpF, dot_map) {
    if (p_rank == 0)
        a.put(0, dim - 1, range_alpha.data());
    auto ga_dot = a.dot(sparse_array);
    double ref_dot = 0.;
    for (const auto &item : sparse_array) {
        ref_dot += range_alpha[item.first] * item.second;
    }
    auto l = Lock();
    ASSERT_THAT(ga_dot, DoubleEq(ref_dot));
}

TEST_F(ArrayCollectiveOpF, times) {
    a.set(alpha);
    b.set(beta);
    auto c = Array{a};
    c.times(&a, &b, false, true);
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto from_ga_buffer_c = c.vec();
    auto ref_a = std::vector<double>(dim, alpha);
    auto ref_b = std::vector<double>(dim, beta);
    auto ref_c = std::vector<double>(dim, alpha * beta);
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
}

// c[i] += a[i]/(b[i]-shift)
TEST_F(ArrayCollectiveOpF, divide_append_negative) {
    double shift = 0.5;
    a.set(alpha);
    b.set(beta);
    auto c = Array{a};
    c.divide(&a, &b, shift, true, true, false, true);
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto from_ga_buffer_c = c.vec();
    auto ref_a = std::vector<double>(dim, alpha);
    auto ref_b = std::vector<double>(dim, beta);
    auto ref_c = std::vector<double>(dim, alpha + alpha / (beta - shift));
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
}

// c[i] += a[i]/(b[i]+shift)
TEST_F(ArrayCollectiveOpF, divide_append_positive) {
    double shift = 0.5;
    a.set(alpha);
    b.set(beta);
    auto c = Array{a};
    c.divide(&a, &b, shift, true, false, false, true);
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto from_ga_buffer_c = c.vec();
    auto ref_a = std::vector<double>(dim, alpha);
    auto ref_b = std::vector<double>(dim, beta);
    auto ref_c = std::vector<double>(dim, alpha + alpha / (beta + shift));
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
}

// c[i] = a[i]/(b[i]+shift)
TEST_F(ArrayCollectiveOpF, divide_overwrite_positive) {
    double shift = 0.5;
    a.set(alpha);
    b.set(beta);
    auto c = Array{a};
    c.divide(&a, &b, shift, false, false, false, true);
    auto from_ga_buffer_a = a.vec();
    auto from_ga_buffer_b = b.vec();
    auto from_ga_buffer_c = c.vec();
    auto ref_a = std::vector<double>(dim, alpha);
    auto ref_b = std::vector<double>(dim, beta);
    auto ref_c = std::vector<double>(dim, alpha / (beta + shift));
    auto l = Lock();
    ASSERT_THAT(from_ga_buffer_a, Pointwise(DoubleEq(), ref_a));
    ASSERT_THAT(from_ga_buffer_b, Pointwise(DoubleEq(), ref_b));
    ASSERT_THAT(from_ga_buffer_c, Pointwise(DoubleEq(), ref_c));
}

} // namespace gci
