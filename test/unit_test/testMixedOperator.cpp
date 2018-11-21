#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>

#include <gciMixedOperator.h>
#include <numeric>

using namespace gci;

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

TEST(MixedOperator, fcidump_constructor){
    std::string fcidump = "fcidump_holstein_2site";
    FCIdump dump{fcidump};
    SymmetryMatrix::Operator Hel = constructOperator(dump);
    MixedOperator ham{FCIdump(fcidump)};
    EXPECT_EQ(ham.nMode, 2);
}

TEST(MixedOperator, expectVal){
//    read Hamiltonian from fcidump files
}

TEST(MixedOperator, iterators){
}
