#include <gtest/gtest.h>
#include "filter.h"

TEST(FilerTest, trimAndCut) {
    Read r("@name",
            "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT",
            "+",
            "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    Options opt;
    opt.qualityCut.enabledFront = true;
    opt.qualityCut.enabledTail = true;
    opt.qualityCut.windowSizeFront = 4;
    opt.qualityCut.qualityFront = 20;
    opt.qualityCut.windowSizeTail = 4;
    opt.qualityCut.qualityTail = 20;
    Filter filter(&opt);
    int frontTrimmed = 0;
    Read* ret = filter.trimAndCut(&r, 0, 1, frontTrimmed);

    EXPECT_EQ(*ret->mSeq, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    EXPECT_EQ(*ret->mQuality, "CCCCCCCCCCC////CCCCCCCCCCCCC");
}
