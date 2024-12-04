#pragma once
#include <string>
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwy {
namespace HWY_NAMESPACE {

template <class D, class Func, typename T = hwy::HWY_NAMESPACE::TFromD<D>>
void Transform1Reversed(D d, T* HWY_RESTRICT inout, size_t count,
                const T* HWY_RESTRICT in1, const Func& func) {
  const size_t N = hwy::HWY_NAMESPACE::Lanes(d);

  size_t idx = 0;
  if (count >= N) {
    for (; idx <= count - N; idx += N) {
      const Vec<D> v = LoadU(d, inout + idx);
      const Vec<D> v1 = LoadU(d, in1 + idx);
      StoreU(Reverse(d, func(d, v, v1)), d, inout + count - N - idx);
    }
  }

  // `count` was a multiple of the vector length `N`: already done.
  if (HWY_UNLIKELY(idx == count)) return;

  const size_t remaining = count - idx;
  HWY_DASSERT(0 != remaining && remaining < N);
  const Vec<D> v = MaskedLoad(FirstN(d, remaining), d, inout + idx);
  const Vec<D> v1 = MaskedLoad(FirstN(d, remaining), d, in1 + idx);
  StoreN(
      SlideDownLanes(
          d, Reverse(d, func(d, v, v1)), N - remaining),
      d, inout, remaining);
}

std::string reverseComplement(std::string *HWY_RESTRICT origin) {

    auto length = origin->length();
    const ScalableTag<uint8_t> d;
    const auto sequence = reinterpret_cast<const uint8_t*>(origin->c_str());
    const auto transform = [](const auto d, auto output, const auto sequence) HWY_ATTR
    {
        const auto a = Set(d, 'A');
        const auto t = Set(d, 'T');
        const auto c = Set(d, 'C');
        const auto g = Set(d, 'G');
        const auto n = Set(d, 'N');
        output = IfThenElse(Eq(sequence, a), t, n);
        output = IfThenElse(Eq(sequence, t), a, output);
        output = IfThenElse(Eq(sequence, g), c, output);
        output = IfThenElse(Eq(sequence, c), g, output);
        return output;
    };
    if (length <= 1000000) {
        uint8_t output[length];
        Transform1Reversed(d, output, length, sequence, transform);
        auto retVal = reinterpret_cast<char *>(output);
        std::string reversed(retVal, length);
        return reversed;
    } else {
        const auto allocated = AllocateAligned<uint8_t>(length);
        Transform1Reversed(d, allocated.get(), length, sequence, transform);
        auto retVal = reinterpret_cast<char *>(allocated.get());
        std::string reversed(retVal, length);
        return reversed;
    }
}

}
}
HWY_AFTER_NAMESPACE();
