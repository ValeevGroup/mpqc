#pragma once
#ifndef TCC_TENSOR_DENSITYTENSOR_H
#define TCC_TENSOR_DENSITYTENSOR_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

namespace tcc {
namespace tensor {

template <typename T>
class DensityTensor : public TA::Tensor<T> {
  public:
    operator TA::Tensor<T>() const {return static_cast<TA::Tensor<T>>(*this);}

    DensityTensor() = default;
    DensityTensor(DensityTensor const &) = default;
    DensityTensor(DensityTensor &&) = default;
    DensityTensor &operator=(DensityTensor const &) = default;
    DensityTensor &operator=(DensityTensor &&) = default;

    DensityTensor(TA::Range const &r) : TA::Tensor<T>(r) {}
    DensityTensor(TA::Range const &r, T value) : TA::Tensor<T>(r, value) {}
    DensityTensor(TA::Tensor<T> const &t) : TA::Tensor<double>(t) {}
};

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DENSITYTENSOR_H
