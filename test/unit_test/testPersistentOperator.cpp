#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/gci/gciPersistentOperator.h>
#include <molpro/gci/gciRun.h>
#include <molpro/gci/gciUtils.h>

#include "parallel_utils.h"

using ::testing::ContainerEq;
using ::testing::Pointwise;

using molpro::gci::PersistentOperator;

namespace {
auto data_path = std::string(DATA_PATH);
}
class PersistentOperatorDataF : public ::testing::Test {
public:
  PersistentOperatorDataF()
      : fname_op_rank1(data_path + "/he2_1el.fcidump"), fname_op_rank2(data_path + "/he2.fcidump"),
        fname_hdf5("test_PersistentOperatorDataF.h5"),
        file_id(molpro::gci::utils::open_hdf5_file(fname_hdf5, molpro::gci::mpi_comm_compute, true)) {
    if (!molpro::gci::utils::file_exists(fname_op_rank2))
      throw std::runtime_error("PersistentOperatorDataF::PersistentOperatorDataF() test file not found");
    molpro::FCIdump dump1(fname_op_rank1);
    op_rank1 = std::make_shared<molpro::Operator>(molpro::gci::constructOperator(dump1));
    op_rank1->m_description = "Rank 1 operator";
    molpro::FCIdump dump2(fname_op_rank2);
    op_rank2 = std::make_shared<molpro::Operator>(molpro::gci::constructOperator(dump2));
    op_rank2->m_description = "Rank 2 operator";
  }

  ~PersistentOperatorDataF() {
    H5Fclose(file_id);
    if (molpro::gci::parallel_rank == 0)
      std::remove(fname_hdf5.c_str());
  }

  std::string fname_op_rank1;
  std::string fname_op_rank2;
  std::string fname_hdf5;
  hid_t file_id;
  std::shared_ptr<molpro::Operator> op_rank1;
  std::shared_ptr<molpro::Operator> op_rank2;
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
  PersistentOperator p_op;
  if (molpro::gci::parallel_rank == 0)
    p_op = PersistentOperator(op_rank2, op_rank2->m_description, 0, file_id);
  else
    p_op = PersistentOperator(nullptr, op_rank2->m_description, 0, file_id);
  {
    auto l = Lock();
    EXPECT_EQ(p_op.description(), op_rank2->m_description);
    EXPECT_TRUE(molpro::gci::utils::file_exists(fname_hdf5)) << "operator should be stored on hdf5";
    auto op = p_op.get();
    EXPECT_NE(op, op_rank2) << "get() creates a new object";
    // TODO check for equality of buffers
  }
  GA_Sync();
}

TEST_F(PersistentOperatorDataF, construct_from_hdf5) {
  if (molpro::gci::parallel_rank == 0)
    auto p_op = PersistentOperator(op_rank2, op_rank2->m_description, 0, file_id);
  else
    auto p_op = PersistentOperator(nullptr, op_rank2->m_description, 0, file_id);
  {
    auto l = Lock();
    EXPECT_TRUE(molpro::gci::utils::file_exists(fname_hdf5)) << "operator should be stored on hdf5";
    auto p_op = PersistentOperator(file_id, op_rank2->m_description);
    auto op = p_op.get();
    EXPECT_NE(op, op_rank2) << "get() creates a new object";
    // TODO check for equality of buffers
  }
  GA_Sync();
}
