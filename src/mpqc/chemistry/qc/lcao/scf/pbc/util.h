#ifndef PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_
#define PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace pbc {
namespace detail {

/*!
 * \brief This reduces the lattice range of a rank-2 array whose
 * row index is orbital basis in the reference lattice and
 * column index is orbital basis in a range of lattices.
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

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
using Basis = ::mpqc::lcao::gaussian::Basis;
/*!
 * \brief This computes non-negligible shell pair list; ; shells \c i and \c j
 * form a non-negligible pair if they share a center or the Frobenius norm of
 * their overlap is greater than threshold
 * \param basis1 a basis
 * \param basis2 a basis
 * \param threshold
 *
 * \return a list of pairs with
 * key: shell index
 * mapped value: a vector of shell indices
 */
shellpair_list_t parallel_compute_shellpair_list(
    madness::World &world, const Basis &basis1, const Basis &basis2,
    double threshold = 1e-12,
    double target_precision = std::numeric_limits<double>::epsilon());

}  // namespace detail
}  // namespace pbc
}  // namespace mpqc

#endif  // PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_
