#include "sequence.h"
#include <hwy/highway.h>
#include <hwy/contrib/algo/transform-inl.h>
#include <hwy/aligned_allocator.h>
#include "simdutil.h"

namespace hn = hwy::HWY_NAMESPACE;

Sequence::Sequence(){
}

Sequence::Sequence(string* seq){
    mStr = seq;
}

Sequence::~Sequence(){
    if(mStr)
        delete mStr;
}

void Sequence::print(){
    std::cerr << *mStr;
}

int Sequence::length(){
    return mStr->length();
}

string Sequence::reverseComplement(string *HWY_RESTRICT origin) {
    return hn::reverseComplement(origin);
}

Sequence Sequence::reverseComplement() {
    string* reversed = new string(Sequence::reverseComplement(mStr));
    return Sequence(reversed);
}

Sequence Sequence::operator~(){
    return reverseComplement();
}
