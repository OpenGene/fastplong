#include "sequence.h"
#include <hwy/highway.h>
#include <hwy/contrib/algo/transform-inl.h>
#include <hwy/aligned_allocator.h>
#include "util.h"

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

string Sequence::reverseComplement(string* origin) {

    auto length = origin->length();
    const hn::ScalableTag<uint8_t> d;
    const auto sequence = reinterpret_cast<const uint8_t*>(origin->c_str());
    const auto transform = [](const auto d, auto output, const auto sequence) HWY_ATTR
    {
        const auto a = hn::Set(d, 'A');
        const auto t = hn::Set(d, 'T');
        const auto c = hn::Set(d, 'C');
        const auto g = hn::Set(d, 'G');
        const auto n = hn::Set(d, 'N');
        output = hn::IfThenElse(hn::Eq(sequence, a), t, n);
        output = hn::IfThenElse(hn::Eq(sequence, t), a, output);
        output = hn::IfThenElse(hn::Eq(sequence, g), c, output);
        output = hn::IfThenElse(hn::Eq(sequence, c), g, output);
        return output;
    };
    if (length <= 1000000) {
        uint8_t output[length];
        Transform1Reversed(d, output, length, sequence, transform);
        auto retVal = reinterpret_cast<char *>(output);
        std::string reversed(retVal, length);
        return reversed;
    } else {
        const auto allocated = hwy::AllocateAligned<uint8_t>(length);
        Transform1Reversed(d, allocated.get(), length, sequence, transform);
        auto retVal = reinterpret_cast<char *>(allocated.get());
        std::string reversed(retVal, length);
        return reversed;
    }
}

Sequence Sequence::reverseComplement() {
    string* reversed = new string(Sequence::reverseComplement(mStr));
    return Sequence(reversed);
}

Sequence Sequence::operator~(){
    return reverseComplement();
}
