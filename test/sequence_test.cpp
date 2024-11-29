#include <gtest/gtest.h>
#include "sequence.h"

TEST(SequenceTests, reverse) {
  Sequence s(new string("AAAATTTTCCCCGGGG"));
  Sequence rc = ~s;
  EXPECT_EQ(*s.mStr, "AAAATTTTCCCCGGGG");
  EXPECT_EQ(*rc.mStr, "CCCCGGGGAAAATTTT");
}
