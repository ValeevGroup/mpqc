#ifndef PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_
#define PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace pbc {
namespace detail {

/*!
 * \brief This reduces the lattice range of a rank-2 array whose
 * row index is orbital basis in the reference unit cell and
 * column index is orbital basis in a range of unit cells.
 *
 * \param arg_array argument array
 * \param arg_range lattice range of argument array
 * \param result_range lattice range of result array
 * \return a rank-2 array with reduced lattice range in columns
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> reduced_size_array(
    TA::DistArray<Tile, Policy> const &arg_array, Vector3i const &arg_range,
    Vector3i const &result_range) {
  // check that argument lattice range is greater than result lattice range
  assert(arg_range(0) >= result_range(0) && arg_range(1) >= result_range(1) &&
         arg_range(2) >= result_range(2));
  auto &world = arg_array.world();

  using ::mpqc::lcao::detail::direct_ord_idx;
  using ::mpqc::lcao::detail::direct_3D_idx;
  using ::mpqc::lcao::detail::extend_trange1;

  // # of lattices corresponding to lattice ranges
  const auto arg_range_size =
      1 + direct_ord_idx(arg_range(0), arg_range(1), arg_range(2), arg_range);
  const auto result_range_size =
      1 + direct_ord_idx(result_range(0), result_range(1), result_range(2),
                         result_range);

  auto arg_tiles_range = arg_array.trange().tiles_range();
  // check that argument array is a rank-2 tensor
  assert(arg_tiles_range.rank() == 2u);
  // check that argument array and argument lattice range match
  assert(arg_tiles_range.extent(1) ==
         arg_tiles_range.extent(0) * arg_range_size);

  // make tiled ranges of the result array
  const auto tr0 = arg_array.trange().data()[0];
  const auto tr1 = extend_trange1(tr0, result_range_size);

  const auto ext = tr0.extent();
  const auto arg_eig = array_ops::array_to_eigen(arg_array);
  RowMatrixXd result_eig(ext, ext * result_range_size);
  for (auto result_ord = 0; result_ord != result_range_size; ++result_ord) {
    const auto result_3D = direct_3D_idx(result_ord, result_range);
    const auto arg_ord =
        direct_ord_idx(result_3D(0), result_3D(1), result_3D(2), arg_range);

    result_eig.block(0, result_ord * ext, ext, ext) =
        arg_eig.block(0, arg_ord * ext, ext, ext);
  }

  return array_ops::eigen_to_array<Tile, Policy>(world, result_eig, tr0, tr1);
}

/*!
 * \brief This enlarges the lattice range of a rank-2 array whose
 * row index is orbital basis in the reference unit cell and
 * column index is orbital basis in a range of unit cells.
 *
 * \param arg_array argument array
 * \param arg_range lattice range of argument array
 * \param result_range lattice range of result array
 * \return a rank-2 array with enlarged lattice range in columns
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> enlarged_size_array(
    TA::DistArray<Tile, Policy> const &arg_array, Vector3i const &arg_range,
    Vector3i const &result_range) {
  // check that argument lattice range is smaller than result lattice range
  assert(arg_range(0) <= result_range(0) && arg_range(1) <= result_range(1) &&
         arg_range(2) <= result_range(2));
  auto &world = arg_array.world();

  using ::mpqc::lcao::detail::direct_ord_idx;
  using ::mpqc::lcao::detail::direct_3D_idx;
  using ::mpqc::lcao::detail::extend_trange1;

  // # of lattices corresponding to lattice ranges
  const auto arg_range_size = 1 + direct_ord_idx(arg_range, arg_range);
  const auto result_range_size = 1 + direct_ord_idx(result_range, result_range);

  auto arg_tiles_range = arg_array.trange().tiles_range();
  // check that argument array is a rank-2 tensor
  assert(arg_tiles_range.rank() == 2u);
  // check that argument array and argument lattice range match
  assert(arg_tiles_range.extent(1) ==
         arg_tiles_range.extent(0) * arg_range_size);

  // make tiled ranges of the result array
  const auto tr0 = arg_array.trange().data()[0];
  const auto tr1 = extend_trange1(tr0, result_range_size);

  const auto ext = tr0.extent();
  const auto arg_eig = array_ops::array_to_eigen(arg_array);
  RowMatrixXd result_eig(ext, ext * result_range_size);
  result_eig.setZero();
  for (auto arg_ord = 0; arg_ord != arg_range_size; ++arg_ord) {
    const auto arg_3D = direct_3D_idx(arg_ord, arg_range);
    const auto result_ord = direct_ord_idx(arg_3D, result_range);
    result_eig.block(0, result_ord * ext, ext, ext) =
        arg_eig.block(0, arg_ord * ext, ext, ext);
  }

  return array_ops::eigen_to_array<Tile, Policy>(world, result_eig, tr0, tr1);
}

/*!
 * \brief This performs dot product of two matrices.
 * The left and right matrices may have different # of columns due to different
 * lattice ranges.
 *
 * \param L left array
 * \param R right array
 * \param L_range lattice range of the left array
 * \param R_range lattice range of the right array
 * \return
 */
template <typename Tile, typename Policy>
double dot_product(TA::DistArray<Tile, Policy> const &L,
                   TA::DistArray<Tile, Policy> const &R,
                   Vector3i const &L_range, Vector3i const &R_range) {
  double result;

  if (L_range == R_range) {
    result = L("m, n") * R("m, n");
  } else if (L_range(0) <= R_range(0) && L_range(1) <= R_range(1) &&
             L_range(2) <= R_range(2)) {
    auto reduced_R = reduced_size_array(R, R_range, L_range);
    result = L("m, n") * reduced_R("m, n");
  } else if (L_range(0) >= R_range(0) && L_range(1) >= R_range(1) &&
             L_range(2) >= R_range(2)) {
    auto reduced_L = reduced_size_array(L, L_range, R_range);
    result = reduced_L("m, n") * R("m, n");
  } else
    throw "invalid lattice ranges";

  return result;
}

/*!
 * \brief This performs addition of two matrices.
 * The left and right matrices may have different # of columns due to different
 * lattice ranges.
 *
 * \param L left array
 * \param R right array
 * \param L_range lattice range of the left array
 * \param R_range lattice range of the right array
 * \return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> add(TA::DistArray<Tile, Policy> const &L,
                                TA::DistArray<Tile, Policy> const &R,
                                Vector3i const &L_range,
                                Vector3i const &R_range,
                                const double preL = 1.0,
                                const double preR = 1.0) {
  TA::DistArray<Tile, Policy> result;
  if (L_range == R_range) {
    result("m, n") = preL * L("m, n") + preR * R("m, n");
  } else if (L_range(0) <= R_range(0) && L_range(1) <= R_range(1) &&
             L_range(2) <= R_range(2)) {
    auto enlarged_L = enlarged_size_array(L, L_range, R_range);
    result("m, n") = preL * enlarged_L("m, n") + preR * R("m, n");
  } else if (L_range(0) >= R_range(0) && L_range(1) >= R_range(1) &&
             L_range(2) >= R_range(2)) {
    auto enlarged_R = enlarged_size_array(R, R_range, L_range);
    result("m, n") = preL * L("m, n") + preR * enlarged_R("m, n");
  } else {
    throw "invalid lattice ranges";
  }

  return result;
}

}  // namespace detail
}  // namespace pbc
}  // namespace mpqc

#endif  // PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_
