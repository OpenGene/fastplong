#include <gtest/gtest.h>
#include "sequence.h"
#include "util.h"

TEST(SequenceTests, Negative) {
  Sequence s(new string("AAAATTTTCCCCGGGG"));
  Sequence rc = ~s;
  // EXPECT_EQ(s.mStr, new string("AAAATTTTCCCCGGGG"));
  // EXPECT_EQ(rc.mStr, new string("CCCCGGGGAAAATTTT"));
}
