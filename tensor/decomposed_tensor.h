#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSOR_H
#define TCC_TENSOR_DECOMPOSEDTENSOR_H

#include "../common/namespaces.h"
#include "../include/tiledarray.h"
#include <vector>

namespace tcc {
namespace tensor {

template <typename T>
class DecomposedTensor {
    double cut_ = 1e-7;
    std::vector<TA::Tensor<T>> tensors_;

    // For shorting serialization function name
    template <typename U>
    using m_enable_if_t = typename madness::enable_if<U>::type;

  public:
    using numeric_type = T;
    using value_type = TA::Tensor<T>;

    DecomposedTensor() = default;
    ~DecomposedTensor() = default;
    DecomposedTensor(DecomposedTensor const &) = default;
    DecomposedTensor(DecomposedTensor &&) = default;
    DecomposedTensor &operator=(DecomposedTensor const &) = default;
    DecomposedTensor &operator=(DecomposedTensor &&) = default;

    DecomposedTensor(double c) : cut_(c) {}
    DecomposedTensor(double c, std::vector<TA::Tensor<T>> ts)
            : cut_(c), tensors_{std::move(ts)} {}

    template <typename... Tensors>
    DecomposedTensor(double c, Tensors &&... ts)
            : cut_(c), tensors_{std::forward<Tensors>(ts)...} {}

    double cut() const { return cut_; }
    std::size_t ndecomp() const { return tensors_.size(); }

    std::vector<TA::Tensor<T>> &tensors() { return tensors_; }
    std::vector<TA::Tensor<T>> const &tensors() const { return tensors_; }
    TA::Tensor<T> &tensor(std::size_t i) { return tensors_[i]; }
    TA::Tensor<T> const &tensor(std::size_t i) const { return tensors_[i]; }


    // I wish this didn't have to be here, but for ease of writing it does
    template <typename Archive>
    m_enable_if_t<madness::archive::is_input_archive<Archive>>
    serialize(Archive &ar) {
        assert(false);
    }

    template <typename Archive>
    m_enable_if_t<madness::archive::is_output_archive<Archive>>
    serialize(Archive &ar) {
        assert(false);
    }
};


} // namespace tensor
} // namespace tcc

#endif //  TCC_TENSOR_DECOMPOSEDTENSOR_H
