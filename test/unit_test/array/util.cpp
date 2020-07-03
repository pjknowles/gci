#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "molpro/gci/array/util.h"

using molpro::gci::array::util::Distribution;
using ::testing::ContainerEq;

TEST(Distribution, constructor) {
  auto d = Distribution(5, 17);
  auto ref_distribution = std::vector<std::pair<unsigned long, size_t>>{{0, 4}, {4, 4}, {8, 3}, {11, 3}, {14, 3}};
  ASSERT_THAT(d.proc_buffer, ContainerEq(ref_distribution));
}

TEST(Distribution, process_map) {
  auto d = Distribution(3, 20);
  auto ref_lo_hi = std::vector<std::pair<unsigned long, unsigned long>>{{0, 0},   {1, 5},   {3, 7},   {5, 9},
                                                                        {13, 14}, {15, 19}, {18, 12}, {17, 6}};
  auto ref_process_pairs =
      std::vector<std::pair<int, int>>{{0, 0}, {0, 0}, {0, 1}, {0, 1}, {1, 2}, {2, 2}, {2, 1}, {2, 0}};
  auto process_pairs = std::vector<std::pair<int, int>>(ref_lo_hi.size());
  std::transform(ref_lo_hi.begin(), ref_lo_hi.end(), process_pairs.begin(),
                 [&d](auto p) { return d.process_map(p.first, p.second); });
  ASSERT_THAT(process_pairs, ContainerEq(ref_process_pairs));
}
