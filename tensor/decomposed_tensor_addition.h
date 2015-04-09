#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORADDITION_H
#define TCC_TENSOR_DECOMPOSEDTENSORADDITION_H

#include "decomposed_tensor.h"

namespace tcc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
add(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r) {
    if(l.ndecomp() >= 2){
        if(r.ndecomp() >= 2){
            auto const &l_extent = l.tensor(0).range().size();
            auto const &r_extent = r.tensor(1).range().size();
            const auto l_rank = l.rank();
            const auto r_rank = r.rank();
            const auto out_rank = l_rank + r_rank;
            TA::Range l_new_range{l_extent[0], out_rank};
            TA::Range r_new_range{out_rank, r_extent[1], r_extent[2]};

            TA::Tensor<T> L(l_new_range);
            TA::Tensor<T> R(r_new_range);

            // copy into left
            const auto l_vol = L.range().volume();
            auto data = L.data();
            auto ll_data = l.tensor(0).data();
            auto rl_data = r.tensor(0).data();
            // Access data in stride
            for(auto i = 0ul; i < l_vol; i += out_rank){
                std::copy(ll_data, ll_data + l_rank, data + i);
                ll_data += l_rank;
                std::copy(rl_data, rl_data + r_rank, data + i + l_rank);
                rl_data += r_rank;
            }

            // copy into R
            data = R.data();
            auto lr_data = l.tensor(1).data();
            auto rr_data = r.tensor(1).data();
            auto old_l_vol = l.tensor(1).range().volume();
            auto old_r_vol = r.tensor(1).range().volume();
            // R is copied by rows so we don't have to loop.
            std::copy(lr_data, lr_data + old_l_vol, data);
            std::copy(rr_data, rr_data + old_r_vol, data + old_l_vol);
            return DecomposedTensor<T>(l.cut(), std::move(L), std::move(R));
        }
    }
}

template <typename T>
DecomposedTensor<T>
add(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r,
    TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
add(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r,
    const T factor) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
add(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r,
    const T factor, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
add(DecomposedTensor<T> const &l, const T factor) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
add(DecomposedTensor<T> const &l, const T factor, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>&
add_to(DecomposedTensor<T> &l, DecomposedTensor<T> const &r) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>&
add_to(DecomposedTensor<T> &l, DecomposedTensor<T> const &r, const T factor) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &
add_to(DecomposedTensor<T> &l, const T factor) {
    assert(false);
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORADDITION_H
