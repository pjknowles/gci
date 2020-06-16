#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <stdexcept>

#include <molpro/gci/gciHProduct.h>

using namespace gci;

class HProductF : public ::testing::Test {
public:
    HProductF() :
            nMode(4), nModal(5), excLvl(3),
            mode1(0), mode2(2), mode3(3),
            modal1(1), modal2(3), modal3(4),
            vibSpace(nMode, nModal, excLvl),
            product({{mode1, modal1}, {mode2, modal2}, {mode3, modal3}}) { };
    int nMode, nModal, excLvl;
    int mode1, mode2, mode3;
    int modal1, modal2, modal3;
    VibSpace vibSpace;
    HProduct product;
};

TEST(HProductT, order_and_check_via_constructor) {
    auto disorderedProduct = HProduct({{3, 1}, {2, 3}, {0, 0}});
    auto orderedProduct = HProduct({{2, 3}, {3, 1}});
    EXPECT_EQ(orderedProduct, disorderedProduct);
    EXPECT_THROW(HProduct({{2, 3}, {2, 4}}), std::logic_error) << "One occupied modal per mode";
};

TEST_F(HProductF, empty_and_excitedModes) {
    auto emptyProduct = HProduct{{{1, 0}, {2, 0}, {3, 0}}};
    EXPECT_TRUE(emptyProduct.empty()) << "All modes in ground state";
    //
    // This is the correct way to construct an empty HProduct (ground state)!
    //
    auto emptyProduct1 = HProduct{};
    EXPECT_TRUE(emptyProduct1.empty()) << "Default empty constructor, assumes ground state";
    //
    // The rest are just sanity checks, and their use is discouraged
    //
    auto emptyProduct2 = HProduct{{}};
    EXPECT_TRUE(emptyProduct2.empty()) << "No modals specified, assumes ground state";
    auto emptyProduct3 = HProduct{{{}}};
    EXPECT_TRUE(emptyProduct3.empty()) << "Empty modal specified, assumes ground state";
    auto emptyProduct4 = HProduct{HProduct::product_t{{}}};
    EXPECT_TRUE(emptyProduct4.empty()) << "Empty product specified, assumes ground state";
    auto p = HProduct::product_t{{{}}};
    auto emptyProduct5 = HProduct{p};
    EXPECT_TRUE(emptyProduct5.empty()) << "Empty product specified with an erroneous extra {} i.e. product_t{{{}}}";
    EXPECT_EQ((std::vector<int>{}), emptyProduct.excitedModes());
    EXPECT_EQ(product.excitedModes(), (std::vector<int>{mode1, mode2, mode3}));
}

TEST_F(HProductF, excLvl) {
    EXPECT_EQ(product.excLvl(), 3) << "Number of excited modes in this product";
    EXPECT_EQ(product.excLvl(mode1), modal1);
    EXPECT_EQ(product.excLvl(mode2), modal2);
    EXPECT_EQ(product.excLvl(mode3), modal3);
    HProduct emptyProduct({{0, 0}});
    EXPECT_EQ(emptyProduct.excLvl(), 0);
    EXPECT_EQ(emptyProduct.excLvl(42), 0) << "Any mode not in the product is assumed to be in the ground state";
}

TEST_F(HProductF, operators_equality_begin_end) {
    auto firstModal = HProduct::modal_t{mode1, modal1};
    ASSERT_EQ(product.begin()[0], firstModal);
    HProduct::product_t copyByIterators;
    ASSERT_NO_THROW((copyByIterators = HProduct::product_t(product.begin(), product.end())));
    auto copyProduct = HProduct(copyByIterators);
    EXPECT_TRUE(copyProduct == product);
    HProduct differentProduct{{{1, 1}, {2, 2}}};
    EXPECT_TRUE(differentProduct != product);
}

TEST_F(HProductF, operator_brackets) {
    // retrieve excited modes
    std::vector<int> excitedModes;
    for (int iMode = 0; iMode < product.excLvl(); ++iMode) {
        excitedModes.push_back(product[iMode][0]);
    }
    ASSERT_THAT(excitedModes, ::testing::ContainerEq(excitedModes));
}

TEST_F(HProductF, raise) {
    auto raisedProduct = product;
    raisedProduct.raise(mode1);
    EXPECT_EQ(raisedProduct.excLvl(mode1), modal1 + 1);
    // Check working with the ground state
    auto empty = HProduct{};
    empty.raise(mode1);
    EXPECT_EQ(empty.excLvl(mode1), 1) << "Raising from the ground state";
    EXPECT_EQ(empty.excLvl(), 1);
    auto singlyExcited = HProduct{HProduct::product_t{{mode1, 1}}};
    EXPECT_EQ(empty, singlyExcited) << "Ensure raised ground state is the same as constructed singly excited product";
}

TEST_F(HProductF, lower) {
    auto loweredProduct = product;
    loweredProduct.lower(mode1);
    EXPECT_EQ(loweredProduct.excLvl(mode1), modal1 - 1);
    EXPECT_EQ(loweredProduct.excLvl(), excLvl - 1) << "Lowered down to ground state";
    loweredProduct.lower(mode1);
    EXPECT_FALSE(loweredProduct.withinSpace(vibSpace)) << "Lowered below the ground state";
    // Check working with the ground state
    auto empty = HProduct{};
    empty.lower(mode1);
    EXPECT_FALSE(loweredProduct.withinSpace(vibSpace)) << "Lowered below the ground state";
}

TEST_F(HProductF, changeModal) {
    auto raisedProduct = product;
    int diff = 2;
    raisedProduct.changeModal(mode1, diff);
    EXPECT_EQ(raisedProduct.excLvl(mode1), modal1 + diff);
    auto loweredProduct = product;
    loweredProduct.changeModal(mode2, -diff);
    EXPECT_EQ(loweredProduct.excLvl(mode2), modal2 - diff);
    loweredProduct.changeModal(mode1, -diff);
    EXPECT_FALSE(loweredProduct.withinSpace(vibSpace)) << "Trying to lower below the ground state";
}

TEST_F(HProductF, withinSpace) {
    EXPECT_TRUE(product.withinSpace(vibSpace)) << "Product was hardcoded to be within the vibrational space";
    product.raise(mode3);
    auto raisedMode3 = product;
    EXPECT_FALSE(raisedMode3.withinSpace(vibSpace)) << "Mode3 is already at highest excited level";

}