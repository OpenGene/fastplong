#include "sequence.h"
#include <hwy/highway.h>
#include <hwy/contrib/algo/transform-inl.h>
#include <hwy/aligned_allocator.h>
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
    const auto sequence = reinterpret_cast<const uint8_t*>(&origin[0]);
    auto output = new uint8_t[length];
    const auto transform = [](const auto d, auto output, const auto sequence) HWY_ATTR
    {
        const auto a = hn::Set(d, 'A');
        const auto t = hn::Set(d, 'T');
        const auto c = hn::Set(d, 'C');
        const auto g = hn::Set(d, 'G');
        output = hn::IfThenElse(hn::Eq(sequence, a), t, output);
        output = hn::IfThenElse(hn::Eq(sequence, t), a, output);
        output = hn::IfThenElse(hn::Eq(sequence, c), g, output);
        output = hn::IfThenElse(hn::Eq(sequence, g), c, output);
        return output;
    };

    const hn::ScalableTag<uint8_t> d;
    Transform1(d, output, length, sequence, transform);

    auto retVal = reinterpret_cast<char *>(output);
    std::string reversed(retVal, length);
    return reversed;
}

Sequence Sequence::reverseComplement() {
    string* reversed = new string(Sequence::reverseComplement(mStr));
    return Sequence(reversed);
}

Sequence Sequence::operator~(){
    return reverseComplement();
}
