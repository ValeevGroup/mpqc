#ifndef MPQC_UTILITY_ARRAYINFO_TENSOR_STORE_H
#define MPQC_UTILITY_ARRAYINFO_TENSOR_STORE_H
#include "../../../../include/tiledarray.h"
#include <string>
#include <array>
#include <iostream>
#include <ostream>

namespace mpqc {
namespace util {

template <typename T>
void write_tensor3D_to_file(TiledArray::Tensor<T> const &t,
                            std::ofstream &file) {
  auto lo = t.range().lobound_data();
  auto hi = t.range().upbound_data();

  for (auto i = lo[0]; i < hi[0]; ++i) {
    for (auto j = lo[1]; j < hi[1]; ++j) {
      for (auto k = lo[2]; k < hi[2]; ++k) {
          const auto val = t(i,j,k);
          if(val >= 1e-16){
            file << i << " " << j << " " << k << " " << val << "\n";
          } else {
            file << i << " " << j << " " << k << " " << 0.0 << "\n";
          }
      }
    }
  }
}

inline void write_array_tuple3D(
    TiledArray::DistArray<TiledArray::Tensor<double>,
                          TiledArray::SparsePolicy> const &A,
    std::string const &output_file_name) {
  const auto volume = A.trange().tiles().volume();

  std::ofstream outfile(output_file_name);

  // Inefficent on multiple nodes
  TiledArray::Tensor<double> dummy;
  for (auto i = 0ul; i < volume; ++i) {
    if (!A.is_zero(i)) {
      dummy = A.find(i);
    } else {
      dummy = TiledArray::Tensor<double>(A.trange().make_tile_range(i), 0.0);
    }
    write_tensor3D_to_file(dummy, outfile);
  }
  outfile.close();
}

inline void write_shape_tuple3D(
    TiledArray::DistArray<TiledArray::Tensor<double>,
                          TiledArray::SparsePolicy> const &A,
    std::string const &output_file_name) {

  std::ofstream outfile(output_file_name);

  write_tensor3D_to_file(A.get_shape().data(), outfile);

  outfile.close();
}

}  //  namespace util
}  // namespace mpqc
#endif  // MPQC_UTILITY_ARRAYINFO_TENSOR_STORE_H
