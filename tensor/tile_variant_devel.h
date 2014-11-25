#pragma once
#ifndef TCC_TENSOR_TILEVARIANT_DEVEL_H
#define TCC_TENSOR_TILEVARIANT_DEVEL_H

#include "../include/btas.h"

namespace tcc {
namespace tensor {

template <typename T>
class TileVariantDevel {
  public:
    TileVariantDevel(btas::Tensor<T> t) : tensor_{std::move(t)} {}

    btas::Tensor<T> const &tensor() const { return tensor_; }
    btas::Tensor<T> &tensor() { return tensor_; }

  private:
    btas::Tensor<T> tensor_;
};

} // namespace tensor
} // namespace tcc


#endif /* end of include guard: TCC_TENSOR_TILEVARIANT_DEVEL_H */
