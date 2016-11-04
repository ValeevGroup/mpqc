

#include "mpqc/math/external/eigen/eigen.h"
#include <tiledarray.h>

using MatrixI =
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct ShapeTracker {
  MatrixI Cshape_;
  MatrixI Fshape_;
  bool init = false;

  void update(TA::Tensor<float> const &Cshape,
              TA::Tensor<float> const &Fshape) {
    if (!init) {
      auto c_extent = Cshape.range().extent();
      auto f_extent = Fshape.range().extent();

      Cshape_.resize(c_extent[0], c_extent[1] * c_extent[2]);
      Cshape_.setZero();
      Fshape_.resize(f_extent[0], f_extent[1] * f_extent[2]);
      Fshape_.setZero();
      init = true;
    }

    auto cd = Cshape.data();
    auto fd = Fshape.data();

    auto cd_ = Cshape_.data();
    auto fd_ = Fshape_.data();

    for (auto i = 0ul; i < Cshape.range().volume(); ++i) {
      if (*(cd + i) > TiledArray::SparseShape<float>::threshold()) {
        *(cd_ + i) = 1;
      }
      if (*(fd + i) > TiledArray::SparseShape<float>::threshold()) {
        *(fd_ + i) = 1;
      }
    }
  }
};
