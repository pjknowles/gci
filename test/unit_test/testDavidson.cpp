#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <molpro/gci/gciDavidson.h>
#include <molpro/gci/gciMixedWavefunction.h>

#include "parallel_utils.h"

using ::testing::Pointwise;
using ::testing::DoubleEq;
using ::testing::ContainerEq;
using ::testing::Each;

namespace gci {

//!@todo there should be a temporary directory for test files
class DavidsonReadWriteF : public ::testing::Test {
public:
    DavidsonReadWriteF() : dim(100), nS(3), p_rank(GA_Nodeid()), fname("test_read_write_array.h5") {
        ref_values.resize(nS);
        for (auto i = 0ul; i < nS; ++i) {
            for (auto j = 0ul; j < dim; ++j)
                ref_values[i].push_back(i + j);
            aa_write.emplace_back(dim, mpi_comm_compute);
            aa_write.back().zero();
            if (p_rank == 0)
                aa_write.back().put(0, dim - 1, &ref_values[i][0]);
            aa_read.emplace_back(dim, mpi_comm_compute);
            aa_read.back().zero();
        }
        aa_write[0].sync();
    }

    const unsigned int dim;
    const unsigned int nS;
    int p_rank;
    std::string fname;
    std::vector<Array> aa_write;
    std::vector<Array> aa_read;
    std::vector<std::vector<double>> ref_values;
};

TEST_F(DavidsonReadWriteF, read_write_wfn) {
    run::davidson_read_write_wfn<Array>(aa_write, fname, true);
    run::davidson_read_write_wfn<Array>(aa_read, fname, false);
    {
        auto l = Lock();
        for (auto i = 0ul; i < nS; ++i) {
            auto read_values = aa_read[i].vec();
            EXPECT_THAT(read_values, Pointwise(DoubleEq(), ref_values[i]));
        }
    }
}

} // namespace gci
