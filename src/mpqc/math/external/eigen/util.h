#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace detail {

/*!
 * \brief This sorts eigenvalues and eigenvectors
 * in ascending order of the real parts of eigenvalues
 *
 * @param[in,out] eigVal the vector of complex eigenvalues
 * @param[in,out] eigVec the complex matrix consisting of complex eigenvectors
 */
void sort_eigen(VectorZ &eigVal, MatrixZ &eigVec);

/**
 *
 * @tparam _Scalar
 * @tparam _Rows
 * @tparam _Cols
 * @tparam _Options
 * @tparam _MaxRows
 * @tparam _MaxCols
 * @param[in,out] matrix an Eigen::Matrix object
 * @param[in] comparison_tolerance the comparison tolerance; values \f$ a,b \f$ are equal if \f$ ||a|-|b||/(|a|+|b|) \leq \epsilon \f$ where \f$ \epsilon \f$ is @c comparison_tolerance
 * @note A matrix has canonical column phase if the first element with largest absolute magnitude in each column is positive
 */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto canonical_column_phase(Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                            _Scalar comparison_tolerance) {
  assert(comparison_tolerance >= 0);
  // returns true if |a| < |b| (equality is determined by comparison_tolerance)
  auto abs_less = [comparison_tolerance](auto a, auto b) {
    auto abs_a = std::abs(a);
    auto abs_b = std::abs(b);
    if (abs_a + abs_b != 0.0) {
      auto ref_diff = std::abs(abs_a - abs_b) / (abs_a + abs_b);
      return ref_diff < comparison_tolerance ? false : (abs_a < abs_b);
    } else
      return false;  // 0 < 0 is false
  };

  const int ncols = matrix.cols();
  const int nrows = matrix.rows();
  if (nrows > 0) {
    for (int c = 0; c != ncols; ++c) {
      auto col = matrix.col(c);
      auto max_elem = col(0);
      for (int r = 1; r < nrows; ++r) {
        auto elem = col(r);
        if (abs_less(max_elem, elem)) {
          max_elem = elem;
        }
      }
      if (max_elem < 0) {
        col *= -1;
      }
    }
  }

}


/** Creates an Eigen::IOFormat object for printing Eigen matrices in format understood by comma initialization operator.
 *
 * @param[in] matrix an Eigen::Matrix object
 * @param[in] precision number of digits to print; default (-1) is to use full precision
 * @param[in] align_cols a flag that controls whether columns are aligned
 * @return an Eigen::IOFormat object; to use a format object @c fmt to print matrix @c M do @code std::cout << M.format(fmt) << std::endl; @endcode
 * @sa mpqc::detail::write_commainit
 * @note Eigen::IOFormat is described at https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
 * @note comma initialization is described at https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
 * */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto format_commainit(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                      int precision = -1, bool align_cols = true) {
  return matrix.format(Eigen::IOFormat(precision > 0 ? precision : Eigen::FullPrecision,
                                       align_cols ? 0 : Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";"));
}

/** Writes an Eigen matrix to a stream in format understood by comma initialization operator.
 *
 * @param[in,out] os a std::ostream object reference
 * @param[in] matrix an Eigen::Matrix object
 * @param[in] name a string identifying this matrix
 * @param[in] precision number of digits to print; default (-1) is to use full precision
 * @param[in] align_cols a flag that controls whether columns are aligned
 * @note comma initialization is described at https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
 * */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void write_commainit(std::ostream& os,
                     const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                     const std::string name = "matrix",
                     int precision = -1, bool align_cols = true) {
  os << name << "(" << matrix.rows() << "," << matrix.cols() << "); "
     << name << format_commainit(matrix, precision, align_cols);
}

/** Creates an Eigen::IOFormat object for printing Eigen matrices in format useful for static initialization of C++ multidimentional arrays.
 *
 * @param[in] matrix an Eigen::Matrix object
 * @return an Eigen::IOFormat object; to use a format object @c fmt to print matrix @c M do @code std::cout << M.format(fmt) << std::endl; @endcode
 * @note Eigen::IOFormat is described at https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
 * */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto format_cpp(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix) {
  return matrix.format(Eigen::IOFormat(Eigen::FullPrecision, 0, ", ", ",\n", "{", "}", "{", "}"));
}

}  // namespace detail
}  // namespace mpqc

#endif // MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
