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
  const Vec<D> v = LoadN(d, inout + idx, remaining);
  const Vec<D> v1 = LoadN(d, in1 + idx, remaining);
  StoreN(
      SlideDownLanes(
          d, Reverse(d, func(d, v, v1)), N - remaining),
      d, inout, remaining);
}

}
}
HWY_AFTER_NAMESPACE();
