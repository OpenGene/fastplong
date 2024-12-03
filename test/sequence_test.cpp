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
  auto largeSequence = new char[10000000];
  auto largeString = new string(largeSequence, 10000000);

  std::fill(largeString->begin(), largeString->end(), 'A');
  Sequence s(largeString);
  Sequence rc = ~s;

  auto largeResult = new char[10000000];
  auto largeResultString = new string(largeResult, 10000000);
  std::fill(largeResultString->begin(), largeResultString->end(), 'T');
  EXPECT_EQ(*s.mStr, *largeString);
  EXPECT_EQ(*rc.mStr, *largeResultString);
}

