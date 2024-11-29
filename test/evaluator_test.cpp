#include <gtest/gtest.h>
#include "evaluator.h"

TEST(EvaluatorTests, int2Seq) {
    Evaluator eval(NULL);
    string s = "ATCGATCGAT";
    EXPECT_EQ(eval.int2seq(eval.seq2int(s, 0, 10, -1), 10), s);
}
