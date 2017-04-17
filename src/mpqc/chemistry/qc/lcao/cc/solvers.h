#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_

#include "mpqc/chemistry/qc/cc/solvers.h"
#include "mpqc/chemistry/qc/lcao/factory/factory.h"

namespace mpqc {
namespace lcao {
namespace cc {

namespace detail {

template <typename Tile, typename Policy,
          typename EigenVectorX =
              Eigen::Matrix<typename Tile::element_type, Eigen::Dynamic, 1>>
TA::DistArray<Tile, Policy> jacobi_update_t2_abij(
    const TA::DistArray<Tile, Policy>& r2_abij, const EigenVectorX& ens_occ,
    const EigenVectorX& ens_uocc) {
  auto denom = [ens_occ, ens_uocc](Tile& result_tile, const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];
    const auto i0 = result_tile.range().lobound()[2];
    const auto in = result_tile.range().upbound()[2];
    const auto j0 = result_tile.range().lobound()[3];
    const auto jn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::scalar_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens_uocc[a];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens_uocc[b];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens_occ[i];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens_occ[j];
            const auto e_iajb = e_i + e_j - e_a - e_b;
            const auto old = arg_tile[tile_idx];
            const auto result_abij = old / (e_iajb);
            const auto abs_result_abij = std::abs(result_abij);
            norm += abs_result_abij * abs_result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto delta_t2_abij = TA::foreach (r2_abij, denom);
  delta_t2_abij.world().gop.fence();
  return delta_t2_abij;
}

template <typename Tile, typename Policy,
          typename EigenVectorX =
              Eigen::Matrix<typename Tile::element_type, Eigen::Dynamic, 1>>
TA::DistArray<Tile, Policy> jacobi_update_t1_ai(
    const TA::DistArray<Tile, Policy>& r1_ai, const EigenVectorX& ens_occ,
    const EigenVectorX& ens_uocc) {
  auto denom = [ens_occ, ens_uocc](Tile& result_tile, const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto i0 = result_tile.range().lobound()[1];
    const auto in = result_tile.range().upbound()[1];

    auto tile_idx = 0;
    typename Tile::scalar_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens_uocc[a];
      for (auto i = i0; i < in; ++i, ++tile_idx) {
        const auto e_i = ens_occ[i];
        const auto e_ia = e_i - e_a;
        const auto old = arg_tile[tile_idx];
        const auto result_ai = old / (e_ia);
        const auto abs_result_ai = std::abs(result_ai);
        norm += abs_result_ai * abs_result_ai;
        result_tile[tile_idx] = result_ai;
      }
    }
    return std::sqrt(norm);
  };

  auto delta_t1_ai = TA::foreach (r1_ai, denom);
  delta_t1_ai.world().gop.fence();
  return delta_t1_ai;
}

}  // namespace detail

/// JacobiDIISSolver updates the CC T amplitudes using standard Jacobi+DIIS
template <typename T>
class JacobiDIISSolver : public ::mpqc::cc::DIISSolver<T, T> {
 public:
  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object; it will be queried for all keywords of ::mpqc::cc::DIISSolver<T,T> .
   */
  // clang-format on

  JacobiDIISSolver(const KeyVal& kv,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> f_ii,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> f_aa)
      : ::mpqc::cc::DIISSolver<T, T>(kv) {
    std::swap(f_ii_, f_ii);
    std::swap(f_aa_, f_aa);
  }
  virtual ~JacobiDIISSolver() = default;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 1> f_ii_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> f_aa_;

  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    t1("a,i") += detail::jacobi_update_t1_ai(r1, f_ii_, f_aa_)("a,i");
    t2("a,b,i,j") += detail::jacobi_update_t2_abij(r2, f_ii_, f_aa_)("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }
};

/// PNOSolver updates the CC T amplitudes using standard Jacobi+DIIS in PNO
/// space
template <typename T, typename DT>
class PNOSolver : public ::mpqc::cc::DIISSolver<T, T> {
 public:
  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object; it will be queried for all keywords of ::mpqc::cc::DIISSolver , as well
   * as the following additional keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | pno_method | string | standard | The PNO construction method. Valid values are: \c standard . |
   * | tpno | double | 1e-8 | The PNO construction threshold. This non-negative integer specifies the screening threshold for the eigenvalues of the pair density. Setting this to zero will cause the full (untruncated) set of PNOs to be used. |
   * | tosv | double | 1e-9 | The OSV construction threshold. This non-negative integer specifies the screening threshold for the eigenvalues of the pair density of the diagonal pairs. Setting this to zero will cause the full (untruncated) set of OSVs to be used. |
   */
  // clang-format on
  PNOSolver(const KeyVal& kv, Factory<T,DT>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv), factory_(factory),
      pno_method_(kv.value<std::string>("pno_method", "standard")),
      tpno_(kv.value<double>("tpno", 1.e-8)) {
    // compute and store PNOs truncated with threshold tpno_
    assert(false && "not yet implemented");
    // for singles: store OSVs (PNOs for diagonal pairs) screened with threshold tosv_
    assert(false && "not yet implemented");
    // for every pair must compute diagonals of the Fock matrix in the PNO basis for the Jacobi update
    assert(false && "not yet implemented");
    // for singles: same, but in the OSV basis
    assert(false && "not yet implemented");
  }
  virtual ~PNOSolver() = default;

  /// @return PNO truncation threshold
  double tpno() const { return tpno_; }
  /// @return OSV truncation threshold
  double tosv() const { return tosv_; }

 private:
  // note: not update_only since DIIS is done in the PNO subspace, must be followed by
  //       backtransform updated amplitudes to full space
  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    // transform r1 and r2 to OSV and PNO basis, respectively
    assert(false && "not yet implemented");
    // apply jacobi update
    assert(false && "not yet implemented");
    // DIIS extrapolate in OSV/PNO basis
    assert(false && "not yet implemented");
    // back-transform extrapolated t1 and t2 from OSV and PNO to full virtual space
    assert(false && "not yet implemented");
  }
  Factory<T, DT>& factory_;
  std::string pno_method_;  //!< the PNO construction method
  double tpno_;  //!< the PNO truncation threshold
  double tosv_;  //!< the OSV (diagonal PNO) truncation threshold
};

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_ */
