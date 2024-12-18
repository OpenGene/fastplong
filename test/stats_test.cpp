#include <gtest/gtest.h>
#include "../src/stats.h"
#include "../src/evaluator.h"
#include <algorithm>

TEST(StatsTests, summarise) {
	Read* left = new Read(
        new string("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA"),
		new string("TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAG"),
		new string("+"),
		new string("AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"));

    Options* options = new Options();
    Stats stats(options);
    stats.statRead(left);
    stats.summarize();
    EXPECT_EQ(stats.getCycles(), 102);
    EXPECT_EQ(stats.getBases(), 102);
    EXPECT_EQ(stats.getReads(), 1);
    EXPECT_EQ(stats.getQ20(), 101);
    EXPECT_EQ(stats.getQ30(), 100);
    EXPECT_EQ(stats.getGCNumber(), 35);
}