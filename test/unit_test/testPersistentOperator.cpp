#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <gciRun.h>
#include <gciPersistentOperator.h>
#include <gciUtils.h>

#include "parallel_utils.h"

using ::testing::Pointwise;
using ::testing::ContainerEq;

namespace gci {

class PersistentOperatorDataF : public ::testing::Test {
public:
    PersistentOperatorDataF()
            : fname_op_rank1("test/unit_test/data/he2_1el.fcidump"),
              fname_op_rank2("test/unit_test/data/he2.fcidump"),
              fname_hdf5("test_PersistentOperatorDataF.h5"),
              file_id(utils::open_hdf5_file(fname_hdf5, gci::mpi_comm_compute, true)) {
        if (!utils::file_exists(fname_op_rank2))
            throw std::runtime_error("PersistentOperatorDataF::PersistentOperatorDataF() test file not found");
        FCIdump dump1(fname_op_rank1);
        op_rank1 = std::make_shared<SymmetryMatrix::Operator>(constructOperator(dump1));
        op_rank1->m_description = "Rank 1 operator";
        FCIdump dump2(fname_op_rank2);
        op_rank2 = std::make_shared<SymmetryMatrix::Operator>(constructOperator(dump2));
        op_rank2->m_description = "Rank 2 operator";
    }

    ~PersistentOperatorDataF() {
        H5Fclose(file_id);
        if (gci::parallel_rank == 0) std::remove(fname_hdf5.c_str());
    }

    std::string fname_op_rank1;
    std::string fname_op_rank2;
    std::string fname_hdf5;
    hid_t file_id;
    std::shared_ptr<SymmetryMatrix::Operator> op_rank1;
    std::shared_ptr<SymmetryMatrix::Operator> op_rank2;
};

TEST_F(PersistentOperatorDataF, construct_from_op_rank1) {
    {
        auto l = Lock();
        auto p_op = PersistentOperator(op_rank1, op_rank1->m_description, 0, file_id, true);
        EXPECT_EQ(p_op.description(), op_rank1->m_description);
        EXPECT_EQ(p_op.get(), op_rank1) << "should point to the same object";
    }
    GA_Sync();
}

TEST_F(PersistentOperatorDataF, construct_from_op_rank2) {
    std::unique_ptr<PersistentOperator> p_op;
    if (gci::parallel_rank == 0)
        p_op = std::make_unique<PersistentOperator>(op_rank2, op_rank2->m_description, 0, file_id);
    else {
        p_op = std::make_unique<PersistentOperator>(nullptr, op_rank2->m_description, 0, file_id);
    }
    {
        auto l = Lock();
        EXPECT_EQ(p_op->description(), op_rank2->m_description);
        EXPECT_TRUE(utils::file_exists(fname_hdf5)) << "operator should be stored on hdf5";
        auto op = p_op->get();
        EXPECT_NE(op, op_rank2) << "get() creates a new object";
        //TODO check for equality of buffers
    }
    GA_Sync();
}

TEST_F(PersistentOperatorDataF, construct_from_hdf5) {
    if (gci::parallel_rank == 0)
        auto p_op = PersistentOperator(op_rank2, op_rank2->m_description, 0, file_id);
    else
        auto p_op = PersistentOperator(nullptr, op_rank2->m_description, 0, file_id);
    {
        auto l = Lock();
        EXPECT_TRUE(utils::file_exists(fname_hdf5)) << "operator should be stored on hdf5";
        auto p_op = PersistentOperator(file_id, op_rank2->m_description);
        auto op = p_op.get();
        EXPECT_NE(op, op_rank2) << "get() creates a new object";
        //TODO check for equality of buffers
    }
    GA_Sync();
}

} // namespace gci
