#include "unittest.h"
#include "fastqreader.h"
#include "read.h"
#include "filter.h"
#include "adaptertrimmer.h"
#include "polyx.h"
#include "nucleotidetree.h"
#include "evaluator.h"
#include <time.h>

UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= report(Read::test(), "Read::test");
    passed &= report(FastqReader::test(), "FastqReader::test");
    passed &= report(Filter::test(), "Filter::test");
    passed &= report(AdapterTrimmer::test(), "AdapterTrimmer::test");
    passed &= report(PolyX::test(), "PolyX::test");
    passed &= report(NucleotideTree::test(), "NucleotideTree::test");
    passed &= report(Evaluator::test(), "Evaluator::test");
    printf("\n==========================\n");
    printf("%s\n\n", passed?"ALL PASSED":"FAILED");
}

bool UnitTest::report(bool result, string message) {
    printf("%s:%s\n\n", message.c_str(), result?" PASSED":" FAILED");
    return result;
}
