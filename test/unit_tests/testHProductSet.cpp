#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>

#include <gciHProductSet.h>
#include <numeric>

using namespace gci;

class HProductSetF : public ::testing::Test {
public:
    HProductSetF() : nMode(3), nModal(4), modeCouplingLevel(2),
                     vibSpace(nMode, nModal, modeCouplingLevel) { }

    int nMode, nModal, modeCouplingLevel;
    VibSpace vibSpace;

};

TEST_F(HProductSetF, constructor_via_VibSpace) {
    auto prodSet = HProductSet(vibSpace);
    ASSERT_EQ(prodSet.vibSpace(), vibSpace);
    int dim_CI0 = 1;
    int dim_CI1 = nMode * (nModal -1);
    int dim_CI2 = 0;
    for (int iMode = 0; iMode < nMode; ++iMode) {
        for (int iModal = 1; iModal < nModal; ++iModal) {
            for (int jMode = 0; jMode < iMode; ++jMode) {
                for (int jModal = 1; jModal < nModal; ++jModal, ++dim_CI2);
            }
        }
    }
    auto vibExcLvlDim = prodSet.vibExcLvlDim();
    ASSERT_EQ(vibExcLvlDim[0], dim_CI0);
    ASSERT_EQ(vibExcLvlDim[1], dim_CI0 + dim_CI1);
    ASSERT_EQ(vibExcLvlDim[2], dim_CI0 + dim_CI1 + dim_CI2);
    ASSERT_EQ(vibExcLvlDim[2], prodSet.vibDim());
}
