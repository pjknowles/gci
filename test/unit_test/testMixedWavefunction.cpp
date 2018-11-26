#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>

#include <gciMixedOperator.h>
#include <gciMixedWavefunction.h>
#include <numeric>

using namespace gci;

TEST(MixedWavefunction, constructor) {
    std::string fcidump_name = "data/fcidump_holstein_2site";
    FCIdump fcidump{fcidump_name};
    Options opt{fcidump.data()};
    opt.addParameter("FCIDUMP", fcidump_name);
    MixedWavefunction wfn{opt};
}

class MixedWavefunction_Holstein_2site : public ::testing::Test {
public:
    MixedWavefunction_Holstein_2site() : fcidump_name("data/fcidump_holstein_2site_nmodal3"), fcidump(fcidump_name),
                                         opt(fcidump.data()), wfn(opt), ham(fcidump) { }
    std::string fcidump_name;
    FCIdump fcidump;
    Options opt;
    MixedWavefunction wfn;
    MixedOperator ham;
};

TEST_F(MixedWavefunction_Holstein_2site, minlocN) {
    wfn.allocate_buffer();
    wfn.diagonalOperator(ham);
    std::vector<std::pair<int, double>> minVals;
    for (size_t i = 0; i < wfn.dimension(); ++i) {
        minVals.emplace_back(i, wfn.at(i));
    }
    std::stable_sort(minVals.begin(), minVals.end(),
                     [](const auto &el1, const auto &el2) {return el1.second < el2.second;});
    std::vector<int> lowest_ref;
    for (size_t n = 1; n < wfn.dimension(); ++n) {
        lowest_ref.resize(n);
        lowest_ref[n-1] = minVals[n-1].first;
        auto lowest = wfn.minlocN(n);
        ASSERT_THAT(lowest, ::testing::Pointwise(::testing::DoubleEq(), lowest_ref));
    }

}

TEST(MixedWavefunction_and_MixedOperator, holstein_2site_hamiltonian_matrix_elements) {
    std::string fcidump_name = "data/fcidump_holstein_2site";
    FCIdump fcidump{fcidump_name};
    Options opt{fcidump.data()};
    opt.addParameter("FCIDUMP", fcidump_name);
    MixedWavefunction wfn{opt};
    wfn.allocate_buffer();
    MixedOperator ham{fcidump};
    ASSERT_EQ(ham.nMode, 2);
    ASSERT_THAT(ham.freq, ::testing::Pointwise(::testing::DoubleEq(), std::vector<double>{0.5, 0.5}));
    ASSERT_THAT(ham.zpe, ::testing::DoubleEq(0.5));
    ASSERT_TRUE(ham.inc_d1());
    ASSERT_FALSE(ham.inc_d2());
    ASSERT_FALSE(ham.inc_T1());
    ASSERT_FALSE(ham.inc_T2());
    ASSERT_EQ(wfn.vibSpace().nMode, 2);
    ASSERT_EQ(wfn.vibSpace().nModal, 2);
    ASSERT_EQ(wfn.vibSpace().excLvl, 2);
    ASSERT_EQ(wfn.elDim(), 2);
    ASSERT_EQ(wfn.vibDim(), 4);
    ASSERT_EQ(wfn.dimension(), 8);
    std::vector<std::vector<double>> H_ref{{0.5, 1., 0.5, 0., 0.5, 0., 0., 0.},
                                           {1., 0.5, 0., 0.5, 0., 0.5, 0., 0.},
                                           {0.5, 0., 1., 1., 0., 0., 0.5, 0.},
                                           {0., 0.5, 1., 1., 0., 0., 0., 0.5},
                                           {0.5, 0., 0., 0., 1., 1., 0.5, 0.},
                                           {0., 0.5, 0., 0., 1., 1., 0., 0.5},
                                           {0., 0., 0.5, 0., 0.5, 0., 1.5, 1.},
                                           {0., 0., 0., 0.5, 0., 0.5, 1., 1.5}};
    std::vector<std::vector<double>> H{};
    for (int i = 0; i < wfn.dimension(); ++i) H.emplace_back(wfn.dimension(), 0);
    auto action = MixedWavefunction(wfn);
    action.allocate_buffer();
    for (int i = 0; i < wfn.dimension(); ++i) {
        wfn.zero();
        wfn.set(i, 1.0);
        action.zero();
        action.operatorOnWavefunction(ham, wfn);
        for (int j = 0; j < wfn.dimension(); ++j) {
            H[i][j] = action.at(j);
        }
    }
    for (int k = 0; k < wfn.dimension(); ++k) {
        EXPECT_THAT(H[k], ::testing::Pointwise(::testing::DoubleEq(), H_ref[k])) << "Column =" << k;
    }
}
