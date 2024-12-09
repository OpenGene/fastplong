#include <gtest/gtest.h>
#include "../src/polyx.h"

TEST(PolyX, trimPolyX) {
    Read r("@name",
        "ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");

    FilterResult fr(NULL, false);
    PolyX::trimPolyX(&r, &fr, 10);
    // r.print();

    EXPECT_EQ(*r.mSeq, "ATTTT");
    EXPECT_EQ(fr.getTotalPolyXTrimmedReads(), 1);
    EXPECT_EQ(fr.getTotalPolyXTrimmedBases(), 51);
}
