
#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_H_

#include <iostream>
#include <tiledarray.h>
#include <vector>

namespace mpqc {
namespace tensor {

/// Class to hold decomposition of a tensor.  The class can hold an arbitrary
/// number of decomposition, but the math functions in the helper headers
/// will initially be written to only handle a two way decomposition.
template <typename T>
class DecomposedTensor {
  double cut_ = 1e-6;
  std::vector<TA::Tensor<T>> tensors_;

 public:
  using numeric_type = T;
  using value_type = TA::Tensor<T>;

  DecomposedTensor() : tensors_(){};
  // ~DecomposedTensor() = default;
  // DecomposedTensor(DecomposedTensor const &) = default;
  // DecomposedTensor(DecomposedTensor &&) = default;
  // DecomposedTensor &operator=(DecomposedTensor const &) = default;
  // DecomposedTensor &operator=(DecomposedTensor &&) = default;

  DecomposedTensor(double c) : cut_(c) {}
  DecomposedTensor(double c, std::vector<TA::Tensor<T>> ts)
      : cut_(c), tensors_{std::move(ts)} {}

  template <typename... Tensors>
  DecomposedTensor(double c, Tensors &&... ts)
      : cut_(c), tensors_({std::forward<Tensors>(ts)...}) {}

  // Will default init cut, this is here to assist with set_all_local
  // functionality
  DecomposedTensor(TA::Range const &r, T value)
      : tensors_{TA::Tensor<T>(r, value)} {}

  double cut() const { return cut_; }
  std::size_t ndecomp() const { return tensors_.size(); }
  bool empty() const { return tensors_.empty(); }

  // only works for two way atm Get the right dimension of the first tensor.
  std::size_t rank() const {
    assert(!empty());
    return tensors_[0].range().extent()[1];
  }

  std::vector<std::size_t> orders() const {
    std::vector<std::size_t> o;
    o.reserve(ndecomp());
    for (auto i = 0ul; i < ndecomp(); ++i) {
      o.push_back(tensors_[i].range().rank());
    }
    return o;
  }

  std::vector<TA::Tensor<T>> &tensors() { return tensors_; }
  std::vector<TA::Tensor<T>> const &tensors() const { return tensors_; }

  TA::Tensor<T> &tensor(std::size_t i) { return tensors_[i]; }
  TA::Tensor<T> const &tensor(std::size_t i) const { return tensors_[i]; }

  TA::Range const &range(std::size_t i) const { return tensors_[i].range(); }

  template <typename Archive>
  typename std::enable_if<
      madness::archive::is_output_archive<Archive>::value>::type
  serialize(Archive &ar) {
    double thresh = cut();
    ar &thresh;
    std::size_t ntensors = tensors_.size();
    ar &ntensors;
    for (auto const &t : tensors_) {
      ar &t;
    }
  }

  template <typename Archive>
  typename std::enable_if<
      madness::archive::is_input_archive<Archive>::value>::type
  serialize(Archive &ar) {
    ar &cut_;
    std::size_t ntensors = 0;
    ar &ntensors;
    for (auto i = 0ul; i < ntensors; ++i) {
      TA::Tensor<T> temp;
      ar &temp;
      tensors_.push_back(temp);
    }
  }
};

extern template class DecomposedTensor<double>;
extern template DecomposedTensor<double>::DecomposedTensor();

template <typename T>
std::ostream &operator<<(std::ostream &os, DecomposedTensor<T> const &t) {
  std::cout << "DecomposedTensor: Cut = " << t.cut()
            << ". Times Decomposed = " << t.ndecomp() << ". Tensors:\n";
  for (auto i = 0ul; i < t.ndecomp(); ++i) {
    std::cout << "\t" << t.tensor(i) << std::endl;
  }
  return os;
}

}  // namespace tensor
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_H_
