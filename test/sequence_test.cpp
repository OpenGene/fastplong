#include <gtest/gtest.h>
#include "../src/sequence.h"
#include <algorithm>

TEST(SequenceTests, reversecomplmeent) {
  Sequence s(new string("AAAATTTTCCCCGGGG"));
  Sequence rc = ~s;
  EXPECT_EQ(*s.mStr, "AAAATTTTCCCCGGGG");
  EXPECT_EQ(*rc.mStr, "CCCCGGGGAAAATTTT");
}

TEST(SequenceTests, reversecomplement_nonpadded) {
  Sequence s(new string("AAAATTTTCCCCGGGGCG"));
  Sequence rc = ~s;
  EXPECT_EQ(*s.mStr, "AAAATTTTCCCCGGGGCG");
  EXPECT_EQ(*rc.mStr, "CGCCCCGGGGAAAATTTT");
}

TEST(SequenceTests, reversecomplement_large) {
  auto largeString = new string(10000000, 'A');

  Sequence s(largeString);
  Sequence rc = ~s;

  auto largeResultString = new string(10000000, 'T');
  EXPECT_EQ(*s.mStr, *largeString);
  EXPECT_EQ(*rc.mStr, *largeResultString);
}

