#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <gciRun.h>
#include <gciPersistentOperator.h>
#include <gciUtils.h>

using ::testing::Pointwise;
using ::testing::ContainerEq;

namespace gci {

class PersistentOperatorDataF : public ::testing::Test {
public:
    PersistentOperatorDataF()
            : fname_op_rank1("test/unit_test/data/he2_1el.fcidump"),
              fname_op_rank2("test/unit_test/data/he2.fcidump"),
              fname_hdf5("test_PersistentOperatorDataF.h5") {
        if (!utils::file_exists(fname_op_rank2))
            throw std::runtime_error("PersistentOperatorDataF::PersistentOperatorDataF() test file not found");
//        else std::cout << "test file found" << std::endl;
        FCIdump dump(fname_op_rank2);
        op_rank2 = std::make_shared<SymmetryMatrix::Operator>(constructOperator(dump));
        op_rank2->m_description = "Rank 2 operator";
    }

    ~PersistentOperatorDataF() {if (utils::file_exists(fname_hdf5)) std::remove(fname_hdf5.c_str());}

    std::string fname_op_rank1;
    std::string fname_op_rank2;
    std::string fname_hdf5;
//    std::shared_ptr<SymmetryMatrix::Operator> op_rank1;
    std::shared_ptr<SymmetryMatrix::Operator> op_rank2;
};

//TEST_F(PersistentOperatorDataF, construct_from_op_rank1) {
//    if (gci::parallel_rank == 0) {
//        auto p_op = PersistentOperator(op_rank1, fname_hdf5);
//        EXPECT_EQ(p_op.description(), op_rank1->m_description);
//        EXPECT_EQ(p_op.get(), op_rank1) << "should point to the same object";
//    }
//}

TEST_F(PersistentOperatorDataF, construct_from_op_rank2) {
    if (gci::parallel_rank == 0) {
        auto p_op = PersistentOperator(op_rank2, fname_hdf5);
        EXPECT_EQ(p_op.description(), op_rank2->m_description);
        EXPECT_TRUE(utils::file_exists(fname_hdf5)) << "operator should be stored on hdf5";
        auto op = p_op.get();
        EXPECT_NE(op, op_rank2) << "get() creates a new object";
        auto bs_ref = op_rank2->bytestream();
        auto bs = op->bytestream();
        EXPECT_THAT(bs.data(), ContainerEq(bs_ref.data()))
                            << "operators should have the same bytestream representation";
    }
    GA_Sync();
}

TEST_F(PersistentOperatorDataF, construct_from_hdf5) {
    if (gci::parallel_rank == 0) {
        {
            auto p_op = PersistentOperator(op_rank2, fname_hdf5);
        }
        EXPECT_TRUE(utils::file_exists(fname_hdf5)) << "operator should be stored on hdf5";
        auto p_op = PersistentOperator(fname_hdf5, op_rank2->m_description);
        auto op = p_op.get();
        EXPECT_NE(op, op_rank2) << "get() creates a new object";
        auto bs_ref = op_rank2->bytestream();
        auto bs = op->bytestream();
        EXPECT_THAT(bs.data(), ContainerEq(bs_ref.data()))
                            << "operators should have the same bytestream representation";
    }
    GA_Sync();
}

} // namespace gci
