#pragma once
#include "../tensor/decomp_algebra.h"
#include "../tensor/decomp_tensor.h"

namespace tcc {
namespace integrals {
namespace compute_functors {

class BtasToDecompTensor {
  public:
    BtasToDecompTensor(double c) : cut_{c} {}

    template <std::size_t N>
    tcc::tensor::DecompTensor operator()(tensor::ShallowTensor<N> const bt) {
        TATensor tensor{bt.range(), bt.tensor().data()};
        auto lr_data = tensor::algebra::ColPivotedQr(tensor, cut_);
        const auto Nelems = tensor.range().volume();
        if (Nelems > lr_data.L_.range().volume()
                     + lr_data.R_.range().volume()) {
            return tcc::tensor::DecompTensor{lr_data.product_range_, lr_data.L_,
                                             lr_data.R_, cut_};
        } else {
            return tcc::tensor::DecompTensor{tensor.range(), tensor, cut_};
        }
    }

  private:
    double cut_;
};

class TaToDecompTensor {
  public:
    TaToDecompTensor(double c) : cut_{c} {}

    tcc::tensor::DecompTensor operator()(TATensor const &t) {
        auto lr_data = tensor::algebra::ColPivotedQr(t, cut_);
        const auto Nelems = t.range().volume();
        if (Nelems > lr_data.L_.range().volume()
                     + lr_data.R_.range().volume()) {
            return tcc::tensor::DecompTensor{lr_data.product_range_, lr_data.L_,
                                             lr_data.R_, cut_};
        } else {
            return tcc::tensor::DecompTensor{t.range(), t, cut_};
        }
    }

  private:
    double cut_;
};

} // namespace compute_functors
} // namespace integrals
} // namespace tcc
