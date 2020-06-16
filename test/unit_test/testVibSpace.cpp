#include <gtest/gtest.h>
#include <stdexcept>

#include <molpro/gci/gciVibSpace.h>

using namespace gci;

TEST(TestVibSpace, constructor_and_equality_operator) {
  int nMode = 1, nModal = 1, modeCoupling = 1;
  auto space1 = VibSpace{nMode, nModal, modeCoupling};
  auto space2 = VibSpace{nMode, nModal, modeCoupling};
  ASSERT_TRUE(space1 == space2);
  int wrongModeCouplingLevel = 2;
  ASSERT_THROW((VibSpace{nMode, nModal, wrongModeCouplingLevel}), std::logic_error);
}
