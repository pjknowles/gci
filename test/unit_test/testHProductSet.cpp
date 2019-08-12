#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <algorithm>

#include <gciHProductSet.h>
#include <numeric>

using namespace gci;

class HProductSetF : public ::testing::Test {
public:
    HProductSetF() : nMode(4), nModal(5), modeCouplingLevel(2),
                     vibSpace(nMode, nModal, modeCouplingLevel),
                     prodSet(vibSpace) {
        if (nMode < 2) throw std::logic_error("Fixture requires at least 2 modes");
        if (nModal < 5) throw std::logic_error("Fixture requires at least 4 modals");
    }

    int nMode, nModal, modeCouplingLevel;
    VibSpace vibSpace;
    HProductSet prodSet;
};

TEST_F(HProductSetF, constructor_and_VibSpace) {
    ASSERT_EQ(prodSet.vibSpace(), vibSpace);
    int dim_CI0 = 1;
    int dim_CI1 = nMode * (nModal - 1);
    int dim_CI2 = 0;
    for (int iMode = 0; iMode < nMode; ++iMode) {
        for (int iModal = 1; iModal < nModal; ++iModal) {
            for (int jMode = 0; jMode < iMode; ++jMode) {
                for (int jModal = 1; jModal < nModal; ++jModal, ++dim_CI2);
            }
        }
    }
    auto excLvlDim = prodSet.excLvlDim();
    ASSERT_EQ(excLvlDim[0], dim_CI0) << "Number of ground state products";
    ASSERT_EQ(excLvlDim[1], dim_CI1) << "Number of products with single mode excited";
    ASSERT_EQ(excLvlDim[2], dim_CI2) << "Number of products with two modes excited";
    ASSERT_EQ(std::accumulate(excLvlDim.begin(), excLvlDim.end(), 0ul), prodSet.vibDim());
}

TEST_F(HProductSetF, generateFullSpace) {
    prodSet.generateFullSpace();
    //! check for duplicates
    auto uniqueCopy = std::vector<HProduct>(prodSet.vibDim(), HProduct{});
    auto newEnd = std::copy_if(prodSet.begin(), prodSet.end(), uniqueCopy.begin(),
                               [&](const HProduct &pIn) {
                                   int count = 0;
                                   for (const auto &el : prodSet) count += pIn == el ? 1 : 0;
                                   return count == 1;
                               });
    ASSERT_EQ(std::distance(uniqueCopy.begin(), newEnd), prodSet.vibDim()) << "Some products are duplicates";
    ASSERT_EQ(uniqueCopy.end(), newEnd) << "Some products are duplicates";
    //! check excitation level is consistent
    auto maxExcLvl = std::max_element(uniqueCopy.begin(), uniqueCopy.end(),
                                      [](const auto &el1, const auto &el2) {return el1.excLvl() < el2.excLvl();});
    ASSERT_EQ(maxExcLvl->excLvl(), vibSpace.excLvl);
    auto excLevel = std::vector<size_t>(vibSpace.excLvl + 1, 0);
    std::for_each(prodSet.begin(), prodSet.end(), [&](const auto &el) {excLevel[el.excLvl()] += 1;});
    ASSERT_EQ(excLevel, prodSet.excLvlDim());
}

TEST_F(HProductSetF, index) {
    prodSet.generateFullSpace();
    // Correct index in the full vibrational space is a sequentially increasing range
    auto correctIndex = std::vector<size_t>(prodSet.vibDim());
    std::iota(correctIndex.begin(), correctIndex.end(), 0ul);
    //
    auto calculatedIndex = std::vector<size_t>(prodSet.vibDim());
    std::transform(prodSet.begin(), prodSet.end(), calculatedIndex.begin(),
                   [&](const auto &prod) {return prodSet.index(prod);});
    ASSERT_THAT(calculatedIndex, ::testing::ContainerEq(correctIndex))
                                << "Wrong index for Hartree products within the full vibrational space. ";
}


TEST_F(HProductSetF, connectedSet_HO_operator) {
    // Connected to the ground state
    prodSet.generateFullSpace();
    auto groundState = prodSet[0];
    auto connectedSet = HProductSet(groundState, vibSpace, VibOp{VibOpType::HO});
    EXPECT_TRUE(connectedSet.connectedSet());
    ASSERT_TRUE(connectedSet.vibDim() == 1) << "Harmonic Oscillator is diagonal";
    ASSERT_EQ(connectedSet[0], groundState) << "Harmonic Oscillator is diagonal";
    // Some random set
    auto randomIndex = prodSet.vibDim() / 2;
    auto refState = prodSet[randomIndex];
    connectedSet = HProductSet(refState, vibSpace, VibOp{VibOpType::HO});
    ASSERT_TRUE(connectedSet.vibDim() == 1) << "Harmonic Oscillator is diagonal";
    ASSERT_EQ(connectedSet[0], refState) << "Harmonic Oscillator is diagonal";
}

class HProductSetPF : public HProductSetF, public ::testing::WithParamInterface<VibOpType> {
};

TEST_P(HProductSetPF, connectedSet_Q_and_dQ_operators) {
    auto op = GetParam();
    prodSet.generateFullSpace();
    // Reference 1: ground state
    auto groundState = prodSet[0];
    const int mode1 = 0, mode2 = 1;
    auto referenceSet = std::vector<HProduct>(1);
    referenceSet[0] = HProduct{{{mode1, 1}}};
    auto connectedSet = HProductSet(groundState, vibSpace, VibOp{op, {mode1}});
    ASSERT_THAT(connectedSet, ::testing::Pointwise(::testing::Eq(), referenceSet));
    // Reference 2: Highest excited state
    auto refProd = HProduct({{mode1, nModal}, {mode2, nModal - 1}});
    referenceSet[0] = HProduct({{mode1, nModal - 1}, {mode2, nModal - 1}});
    connectedSet = HProductSet(refProd, vibSpace, VibOp{op, {mode1}});
    ASSERT_THAT(connectedSet, ::testing::Pointwise(::testing::Eq(), referenceSet));
}

//std::unordered_map<std::type_index, std::string> typeNames(
//        {{typeid(VibOpType::Q), "Q"}, {typeid(VibOpType::dQ), "dQ"}});

INSTANTIATE_TEST_CASE_P(HProductSetF, HProductSetPF, ::testing::Values(VibOpType::Q, VibOpType::dQ),
//                        [typeNames](const auto &op) {return typeNames[typeid(op.param)];});
                        [](const auto &op) {
                            switch (op.param) {
                                case VibOpType::Q: return "Q";
                                case VibOpType::dQ: return "dQ";
                                default: return "?";
                            }
                        });

TEST_F(HProductSetF, connectedSet_Qsq_operator) {
    prodSet.generateFullSpace();
    // Qsq = Q_A * Q_B
    const int mode1 = 0, mode2 = 1;
    auto refProd = HProduct({{mode1, 2}, {mode2, 2}});
    auto referenceSet = std::vector<HProduct>(4);
    referenceSet[0] = HProduct({{mode1, 1}, {mode2, 1}});
    referenceSet[1] = HProduct({{mode1, 1}, {mode2, 3}});
    referenceSet[2] = HProduct({{mode1, 3}, {mode2, 1}});
    referenceSet[3] = HProduct({{mode1, 3}, {mode2, 3}});
    auto connectedSet = HProductSet(refProd, vibSpace, VibOp(VibOpType::Qsq, {mode1, mode2}));
    ASSERT_THAT(connectedSet, ::testing::Pointwise(::testing::Eq(), referenceSet)) << "Qsq = Q_A * Q_B";
    // Qsq = Q_A * Q_A
    referenceSet.resize(3);
    referenceSet[0] = HProduct({{mode1, 0}, {mode2, 2}});
    referenceSet[1] = HProduct({{mode1, 2}, {mode2, 2}});
    referenceSet[2] = HProduct({{mode1, 4}, {mode2, 2}});
    connectedSet = HProductSet(refProd, vibSpace, VibOp(VibOpType::Qsq, {mode1, mode1}));
    ASSERT_THAT(connectedSet, ::testing::Pointwise(::testing::Eq(), referenceSet)) << "Qsq = Q_A * Q_A";
}