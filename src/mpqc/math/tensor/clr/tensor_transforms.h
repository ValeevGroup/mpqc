
#ifndef MPQC_TENSOR_TENSORTRANSFORMS_H
#define MPQC_TENSOR_TENSORTRANSFORMS_H

#include <tiledarray.h>
#include "mpqc/math/tensor/clr/decomposed_tensor.cpp"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.cpp"
#include "mpqc/math/tensor/clr/tile.h"

namespace mpqc {
namespace tensor {

class TaToDecompTensor {
    double clr_thresh_;
    bool compress_;

    Tile<DecomposedTensor<double>> rank_three_convert(TA::Tensor<double> &&t) {
        auto range = t.range();

        auto const extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];
        const auto k = extent[2];
        auto local_range = TA::Range{i, j, k};

        auto tensor = TA::Tensor<double>(local_range, t.data());

        DecomposedTensor<double> dc_tile(clr_thresh_, std::move(tensor));

        if (compress_ && clr_thresh_ != 0.0) {
            auto lr_tile = tensor::algebra::two_way_decomposition(dc_tile);

            if (!lr_tile.empty()) {
                dc_tile = std::move(lr_tile);
            }
        }

        return Tile<DecomposedTensor<double>>(range, std::move(dc_tile));
    }

    Tile<DecomposedTensor<double>> rank_two_convert(TA::Tensor<double> &&t) {
        auto range = t.range();

        const auto extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];

        auto local_range = TA::Range{i, j};
        auto tensor = TA::Tensor<double>(local_range, t.data());

        auto new_tile
              = DecomposedTensor<double>(clr_thresh_, std::move(tensor));

        if (compress_ && clr_thresh_ != 0.0) {
            auto test = tensor::algebra::two_way_decomposition(new_tile);
            if (!test.empty()) {
                new_tile = std::move(test);
            }
        }

        return Tile<DecomposedTensor<double>>(range, std::move(new_tile));
    }

  public:
    using result_type = Tile<DecomposedTensor<double>>;

    TaToDecompTensor(double thresh, bool compress = true)
            : clr_thresh_(thresh), compress_(compress) {}

    Tile<DecomposedTensor<double>> operator()(TA::Tensor<double> &&t) {
        if (t.range().rank() == 3) {
            return rank_three_convert(std::move(t));
        } else if (t.range().rank() == 2) {
            return rank_two_convert(std::move(t));
        } else {
            throw std::invalid_argument("Input tensor to TaToDecompTensor did "
                                        "did not have the correct number of "
                                        "dimensions.");
        }
        return Tile<DecomposedTensor<double>>();
    }

    Tile<DecomposedTensor<double>>
    operator()(madness::Future<TA::Tensor<double>> t_fut) {
        // Just force the future
        TA::Tensor<double> t = t_fut.get();

        if (t.range().rank() == 3) {
            return rank_three_convert(std::move(t));
        } else if (t.range().rank() == 2) {
            return rank_two_convert(std::move(t));
        } else {
            throw std::invalid_argument("Input tensor to TaToDecompTensor did "
                                        "did not have the correct number of "
                                        "dimensions.");
        }
        return Tile<DecomposedTensor<double>>();
    }

    Tile<DecomposedTensor<double>> operator()(TA::Tensor<double> const &t_ref) {
        // Just make a copy of the reference
        TA::Tensor<double> t = t_ref.clone();

        if (t.range().rank() == 3) {
            return rank_three_convert(std::move(t));
        } else if (t.range().rank() == 2) {
            return rank_two_convert(std::move(t));
        } else {
            throw std::invalid_argument("Input tensor to TaToDecompTensor did "
                                        "did not have the correct number of "
                                        "dimensions.");
        }
        return Tile<DecomposedTensor<double>>();
    }

};

class DecompToTaTensor {
  public:
    TA::Tensor<double>
    operator()(Tile<DecomposedTensor<double>> const &t) const {
        auto tensor = tensor::algebra::combine(t.tile());
        auto range = t.range();
        return TA::Tensor<double>(range, tensor.data());
    }
};


} // namespace tensor
} // namespace mpqc

#endif // MPQC_TENSOR_TENSORTRANSFORMS_H
