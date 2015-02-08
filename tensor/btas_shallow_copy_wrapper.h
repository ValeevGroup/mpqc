#pragma once
#ifndef TCC_TENSOR_BTASSHALLOWCOPYWRAPPER_H
#define TCC_TENSOR_BTASSHALLOWCOPYWRAPPER_H

#include "../include/tiledarray.h"
#include "../include/btas.h"

#include <memory>
#include <numeric>

namespace tcc {
namespace tensor {

// I don't really want these types to be public at the moment, this whole
// interface things isn't quite worked out at the moment and will need to
// be refined.
namespace {
template <std::size_t N>
using RangeNd = btas::RangeNd<CblasRowMajor, std::array<long, N>>;

template <std::size_t N>
using TensorType = btas::Tensor<double, RangeNd<N>>;
}

/// ShallTensor wraps a btas::Tensor<double, SpecialRange> with a shared_ptr.
/// At the moment it only holds the type, it is not currently usable in
/// expressions
template <std::size_t N, typename Range = TiledArray::Range>
class ShallowTensor {
  public:
    using eval_type = ShallowTensor;
    using value_type = TensorType<N>;
    using range_type = Range;
    using numeric_type = double;
    using size_type = std::size_t;

    ShallowTensor() = default;
    ~ShallowTensor() = default;
    ShallowTensor(ShallowTensor const &) = default;
    ShallowTensor(ShallowTensor &&) = default;
    ShallowTensor &operator=(ShallowTensor const &) = default;
    ShallowTensor &operator=(ShallowTensor &&) = default;

    ShallowTensor(Range range) : tensor_{}, range_(std::move(range)) {}
    ShallowTensor(Range range, value_type tensor)
        : tensor_{std::make_shared<value_type>(std::move(tensor))},
          range_(std::move(range)) {}

    ShallowTensor clone() const { return ShallowTensor(range_, *tensor_); }

    value_type const &tensor() const {
        assert(!empty());
        return *tensor_;
    }
    value_type &tensor() {
        assert(!empty());
        return *tensor_;
    }

    numeric_type norm() const {
        const auto size = tensor_->size();
        const auto begin = tensor_->data();
        return std::sqrt(
            std::accumulate(begin, begin + size, 0.0,
                            [](double out, double x) { return out + x * x; }));
    }

    bool empty() const { return !tensor_; }
    Range const &range() const { return range_; }

    template <typename Archive>
    void serialize(Archive &ar) {
        assert(false);
    }

  private:
    std::shared_ptr<value_type> tensor_;
    Range range_;
};

} // namespace tensor
} // namespace tcc


#endif /* end of include guard: TCC_TENSOR_BTASSHALLOWCOPYWRAPPER_H */
