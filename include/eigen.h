#pragma once
#ifndef TCC_INCLUDE_EIGEN_H
#define TCC_INCLUDE_EIGEN_H

#pragma GCC diagnostic push
#pragma GCC system_header
#include <Eigen/Dense>
#pragma GCC diagnostic pop

template <typename T>
using RowMatrix
    = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using RowMatrixXd = RowMatrix<double>;

namespace madness {
namespace archive {

template <class>
class archive_array;
template <class T>
inline archive_array<T> wrap(const T*, unsigned int);
template <class Archive, typename Data>
struct ArchiveStoreImpl;
template <class Archive, typename Data>
struct ArchiveLoadImpl;

template <class Archive, typename _T>
struct ArchiveStoreImpl<Archive, RowMatrix<_T>> {
  static inline void store(const Archive& ar, const RowMatrix<_T>& t) {
    ar & t.rows() & t.cols();
    if (t.size())
      ar & madness::archive::wrap(t.data(), t.size());
  }
};

template <class Archive, typename _T>
struct ArchiveLoadImpl<Archive, RowMatrix<_T>> {
  static inline void load(const Archive& ar, RowMatrix<_T>& t) {
    typename RowMatrix<_T>::Index nrows, ncols;
    ar & nrows & ncols;
    t.resize(nrows, ncols);
    if (t.size())
      ar & madness::archive::wrap(t.data(), t.size());
  }
};

}  // namespace archive
}  // namespace madness

#endif // TCC_INCLUDE_EIGEN_H
