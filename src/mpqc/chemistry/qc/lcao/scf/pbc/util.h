#ifndef PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_
#define PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/math/external/tiledarray/util.h"

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

  using ::mpqc::detail::direct_ord_idx;
  using ::mpqc::detail::direct_3D_idx;
  using ::mpqc::detail::extend_trange1;

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

  using ::mpqc::detail::direct_ord_idx;
  using ::mpqc::detail::direct_3D_idx;
  using ::mpqc::detail::extend_trange1;

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

  using ::mpqc::detail::direct_ord_idx;
  auto L_size = 1 + direct_ord_idx(L_range, L_range);
  auto R_size = 1 + direct_ord_idx(R_range, R_range);
  auto L_tiles_range = L.trange().tiles_range();
  auto R_tiles_range = R.trange().tiles_range();
  assert(L_tiles_range.extent(1) == L_tiles_range.extent(0) * L_size);
  assert(R_tiles_range.extent(1) == R_tiles_range.extent(0) * R_size);

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
  } else {
    auto x = std::min(L_range(0), R_range(0));
    auto y = std::min(L_range(1), R_range(1));
    auto z = std::min(L_range(2), R_range(2));
    Vector3i min_range(x, y, z);
    auto reduced_L = reduced_size_array(L, L_range, min_range);
    auto reduced_R = reduced_size_array(R, R_range, min_range);
    result = reduced_L("m, n") * reduced_R("m, n");
  }
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
  using ::mpqc::detail::direct_ord_idx;
  auto L_size = 1 + direct_ord_idx(L_range, L_range);
  auto R_size = 1 + direct_ord_idx(R_range, R_range);
  auto L_tiles_range = L.trange().tiles_range();
  auto R_tiles_range = R.trange().tiles_range();
  assert(L_tiles_range.extent(1) == L_tiles_range.extent(0) * L_size);
  assert(R_tiles_range.extent(1) == R_tiles_range.extent(0) * R_size);

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
    auto x = std::max(L_range(0), R_range(0));
    auto y = std::max(L_range(1), R_range(1));
    auto z = std::max(L_range(2), R_range(2));
    Vector3i max_range(x, y, z);
    auto enlarged_L = enlarged_size_array(L, L_range, max_range);
    auto enlarged_R = enlarged_size_array(R, R_range, max_range);
    result("m, n") = preL * enlarged_L("m, n") + preR * enlarged_R("m, n");
  }

  return result;
}

/*!
 * \brief This prints by-unit-cell Frobenius norms and Infinity norms of a
 * matrix M(μ_0, ν_R).
 *
 * \param M_array a Tiled Array object
 * \param max_lattice_range max range of lattice vector R, indicated by the
 * index of the farthest unit cell (x, y, z) where x, y, and z are integers
 * \param name name of the matrix
 */
template <typename Tile, typename Policy>
void print_norms_by_unit_cell(TA::DistArray<Tile, Policy> const &M_array,
                              Vector3i const &max_lattice_range,
                              std::string const &name) {
  using ::mpqc::detail::direct_3D_idx;
  using ::mpqc::detail::direct_ord_idx;

  const auto elements_range = M_array.trange().elements_range();
  const auto ext0 = elements_range.extent(0);
  const auto ext1 = elements_range.extent(1);
  const auto max_lattice_size = 1 + direct_ord_idx(max_lattice_range, max_lattice_range);
  assert(ext1 / ext0 == max_lattice_size);

  const auto M_eigen = array_ops::array_to_eigen(M_array);
  ExEnv::out0() << "\nNorms of matrix " << name << ":\n";
  for (auto uc_ord = 0ul; uc_ord != max_lattice_size; ++uc_ord) {
    const auto uc_3D = direct_3D_idx(uc_ord, max_lattice_range);
    const auto block = M_eigen.block(0, uc_ord * ext0, ext0, ext0);
    const auto norm_frobenius = block.template lpNorm<2>();
    const auto norm_infty = block.template lpNorm<Eigen::Infinity>();
    ExEnv::out0() << "unit cell (" << uc_3D.transpose()
                  << "), Frobenius norm = " << norm_frobenius
                  << ", Infinity norm = " << norm_infty
                  << std::endl;
  }
}

/*!
 * \brief This extracts a sliced matrix at k = \c k_ord from a large matrix
 * \c M (μ, ν; k). Note that \c M is stored as a 2d array with \cnrows =
 * \c size(obs) and \c ncols = \c size(obs) * \c size(k_points).
 *
 * \param M a Tiled Array object M(μ, ν; k=all)
 * \param k_ord the ordinal number of a specific k point
 * \param nk the total number of k points in each direction
 * \return a Tiled Array object M(μ, ν; k=k_ord)
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> slice_array_at_k(
    const TA::DistArray<Tile, Policy> &M, const size_t k_ord,
    const Vector3i &nk) {
  using ::mpqc::detail::k_ord_idx;

  const auto k_size = 1 + k_ord_idx(nk(0) - 1, nk(1) - 1, nk(2) - 1, nk);
  assert(k_ord >= 0 && k_ord < k_size && "ordinal # of k is out of range");

  auto tr0 = M.trange().data()[0];
  auto tr1 = M.trange().data()[1];
  assert(tr1.extent() == tr0.extent() * k_size);

  const auto tr0_upper = tr0.tiles_range().second;

  typedef std::vector<size_t> block;
  block low{0, tr0_upper * k_ord};
  block up{tr0_upper, tr0_upper * (k_ord + 1)};

  TA::DistArray<Tile, Policy> result;
  result("i, j") = M("i, j").block(low, up);
  return result;
}

/*!
 * \brief This truncates lattice range of an array D(μ, ν_R). If norms of all
 * tiles within a unit cell are below \c threshold, this unit cell will be
 * removed and the array size is shrinked.
 * \param D Tiled array object
 * \param RD_max the original lattice range of the input array
 * \param threshold
 * \return the updated lattice range
 */
template <typename Tile, typename Policy>
Vector3i truncate_lattice_range(const TA::DistArray<Tile, Policy> &D,
                                const Vector3i &RD_max,
                                const double threshold) {
  ExEnv::out0() << "\tUser specified lattice range = " << RD_max.transpose()
                << std::endl;

  using ::mpqc::detail::direct_3D_idx;
  using ::mpqc::detail::direct_ord_idx;

  const auto &Dnorms = D.shape().data();
  std::vector<Vector3i> sig_density_uc_list;
  const auto RD_size = 1 + direct_ord_idx(RD_max, RD_max);

  const auto ntiles = D.trange().tiles_range().extent(0);

  // determine significant unit cells
  for (auto RD_ord = 0ul; RD_ord != RD_size; ++RD_ord) {
    const auto RD_3D = direct_3D_idx(RD_ord, RD_max);
    auto is_significant = false;
    for (auto mu = 0ul; mu != ntiles; ++mu) {
      for (auto nu = 0ul; nu != ntiles; ++nu) {
        const auto nu_in_D = nu + RD_ord * ntiles;
        std::array<size_t, 2> idx{{mu, nu_in_D}};
        if (Dnorms(idx) >= threshold) {
          is_significant = true;
          sig_density_uc_list.emplace_back(RD_3D);
          break;
        }
      }
      if (is_significant) {
        break;
      }
    }
  }

  // renew lattice range based on the list of significant unit cells
  auto x = 0;
  auto y = 0;
  auto z = 0;
  for (const auto &RD_3D : sig_density_uc_list) {
    x = std::max(x, RD_3D(0));
    y = std::max(y, RD_3D(1));
    z = std::max(z, RD_3D(2));
  }
  Vector3i new_RD_max = {x, y, z};
  ExEnv::out0() << "\tUpdated lattice range = " << new_RD_max.transpose()
                << std::endl;

  return new_RD_max;
}

}  // namespace detail
}  // namespace pbc
}  // namespace mpqc

#endif  // PQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_UTIL_H_
