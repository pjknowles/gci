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
//    read Hamiltonian from fcidump files
}

TEST(MixedOperator, expectVal){
//    read Hamiltonian from fcidump files
}

TEST(MixedOperator, iterators){
//    read Hamiltonian from fcidump files
}
