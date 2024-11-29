#include <gtest/gtest.h>
#include "adaptertrimmer.h"

TEST(AdapterTrimmer, trimBySequenceStart) {
    Read r("@name",
        "AGGTGCTGCGCATACTTTTCCACGGGGATACTACTGGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEE///EEEEEEEEEEEEEEEEEEEEE");
    string adapter = "GCGCATACTTTTCCACGGGGATACTACTG";
    int trimmed = AdapterTrimmer::trimBySequenceStart(&r, NULL, adapter, 0.3, 0);
    EXPECT_EQ(*r.mSeq, "GGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT");

    Read r2("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAAGCGCATACTTTTCCACGGGGA",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEET");
    trimmed = AdapterTrimmer::trimBySequenceEnd(&r2, NULL, adapter, 0.3, 0);
    EXPECT_EQ(*r2.mSeq, "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA");

    /*Read read("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    vector<string> adapterList;
    adapterList.push_back("GCTAGCTAGCTAGCTA");
    adapterList.push_back("AAATTTCCCGGGAAATTTCCCGGG");
    adapterList.push_back("ATCGATCGATCGATCG");
    adapterList.push_back("AATTCCGGAATTCCGG");
    trimmed = AdapterTrimmer::trimByMultiSequences(&read, NULL, adapterList);
    if (*read.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG") {
        cerr << read.mSeq << endl;
        return false;
    }*/
}
