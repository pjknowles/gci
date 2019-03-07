#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>

#include <gciMixedOperator.h>
#include <numeric>

using namespace gci;

//
class TestMixedOperator : public ::testing::Test, public MixedOperator {
public:
    TestMixedOperator() {
        freq.push_back(2.0);
        nMode = 1;
    }
};

TEST(VibOp, constructor) {
    auto vibOp = VibOp(VibOpType::HO);
    EXPECT_TRUE(vibOp.mode.empty());
    vibOp = VibOp(VibOpType::Q, {0});
    EXPECT_EQ(vibOp.mode.size(), 1);
    vibOp = VibOp(VibOpType::dQ, {0});
    EXPECT_EQ(vibOp.mode.size(), 1);
    vibOp = VibOp(VibOpType::Qsq, {0, 1});
    EXPECT_EQ(vibOp.mode.size(), 2);
}

TEST(MixedOperator, fcidump_constructor) {
    std::string fcidump = "fcidump_holstein_2site";
    FCIdump dump{fcidump};
    SymmetryMatrix::Operator Hel = constructOperator(dump);
    MixedOperator ham{FCIdump(fcidump)};
    EXPECT_EQ(ham.nMode, 2);
}

TEST(MixedOperator, expectVal) {
//    read Hamiltonian from fcidump files
}

TEST(MixedOperator, iterators) {
}

TEST_F(TestMixedOperator, O_Qsq) {
    int maxM = 8;
    auto op = VibOp(VibOpType::Qsq, {0, 0});
    for (int n = 0; n < maxM; ++n) {
        double ref = 1. / freq[0] * (n + 0.5);
        auto bra = HProduct::t_Product({{0, n}});
        auto el = O_Qsq(bra, bra, op);
        ASSERT_DOUBLE_EQ(el, ref) << "< n | Qsq | n > == 1 / w * (n + 0.5)";
    }
    for (int n = 2; n < maxM; ++n) {
        int m = n - 2;
        double ref = 0.5 / freq[0] * std::sqrt(n * (n - 1));
        auto bra = HProduct::t_Product({{0, n}});
        auto ket = HProduct::t_Product({{0, m}});
        auto el = O_Qsq(bra, ket, op);
        ASSERT_DOUBLE_EQ(el, ref) << "< n | Qsq | m > == 1 / 2 / w * sqrt(n * (n - 1))";
    }
    for (int n = 0; n < maxM - 2; ++n) {
        int m = n + 2;
        if ((n - m) % 2 != 0) continue;
        double ref = 0.5 / freq[0] * std::sqrt(m * (m - 1));
        auto bra = HProduct::t_Product({{0, n}});
        auto ket = HProduct::t_Product({{0, m}});
        auto el = O_Qsq(bra, ket, op);
        ASSERT_DOUBLE_EQ(el, ref) << "< n | Qsq | m > == 1 / 2 / w * sqrt(m * (m - 1))";
    }
}
