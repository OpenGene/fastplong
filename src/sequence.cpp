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
    auto length = origin->length();
    const hn::ScalableTag<uint8_t> d;
    const auto sequence = reinterpret_cast<const uint8_t*>(origin->c_str());
    const auto transform = [](const auto d, auto output, const auto sequence) HWY_ATTR
    {
        const auto A = hn::Set(d, 'A');
        const auto T = hn::Set(d, 'T');
        const auto C = hn::Set(d, 'C');
        const auto G = hn::Set(d, 'G');
        const auto a = hn::Set(d, 'a');
        const auto t = hn::Set(d, 't');
        const auto c = hn::Set(d, 'c');
        const auto g = hn::Set(d, 'g');
        const auto N = hn::Set(d, 'N');

        // output[i] = sequence[i] == 'A' || sequence[i] == 'a' ? 'T' : 'N'
        output = hn::IfThenElse(hn::Or(hn::Eq(sequence, A), hn::Eq(sequence, a)), T, N);
        output = hn::IfThenElse(hn::Or(hn::Eq(sequence, T), hn::Eq(sequence, t)), A, output);
        output = hn::IfThenElse(hn::Or(hn::Eq(sequence, C), hn::Eq(sequence, c)), G, output);
        output = hn::IfThenElse(hn::Or(hn::Eq(sequence, G), hn::Eq(sequence, g)), C, output);
        return output;
    };
    if (length <= 1000000) {
        uint8_t output[length];
        hn::Transform1Reversed(d, output, length, sequence, transform);
        auto retVal = reinterpret_cast<char *>(output);
        std::string reversed(retVal, length);
        return reversed;
    } else {
        const auto allocated = hwy::AllocateAligned<uint8_t>(length);
        hn::Transform1Reversed(d, allocated.get(), length, sequence, transform);
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
