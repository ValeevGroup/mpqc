#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_

#include "mpqc/chemistry/qc/cc/solvers.h"
#include "mpqc/chemistry/qc/lcao/factory/factory.h"

using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;

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

// // JacobiDIISSolver updates the CC T amplitudes using standard Jacobi+DIIS
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
/// @warning This class assumes that the 1- and 2-body amplitudes/residuals
///          given to Solver::update() are laid out as "a,i" and "a,b,i,j",
///          respectively
template <typename T>
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
  PNOSolver(const KeyVal& kv, Factory<T>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv),
        factory_(factory),
        pno_method_(kv.value<std::string>("pno_method", "standard")),
        tpno_(kv.value<double>("tpno", 1.e-8)),
        tosv_(kv.value<double>("tosv", 1.e-9)) {
    // compute and store PNOs truncated with threshold tpno_
    // store PNOs for diagonal pair as OSVs truncated with threshold tosv_

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

    auto& fac = factory_;
    auto& world = fac.world();
    auto& ofac = fac.orbital_registry();

    auto nocc = ofac.retrieve("m").rank();
    // ExEnv::out0() << "nocc = " << nocc;
    auto nocc_act = ofac.retrieve("i").rank();
    auto nvir = ofac.retrieve("a").rank();
    auto nfzc = nocc - nocc_act;

    // Form Fock array
    auto F = fac.compute(L"(p|F|q)");

    // Select just diagonal elements of Fock aray and transform
    // to Eigen vector; use for computing PNOs
    Eigen::VectorXd eps_p = TA::array_to_eigen(F).diagonal();
    auto eps_o = eps_p.segment(nfzc, nocc_act);
    auto eps_v = eps_p.tail(nvir);

    // Transform entire Fock array to Eigen Matrix
    Eigen::MatrixXd F_all = TA::array_to_eigen(F);

    // Select just the occupied portion of the Fock matrix
    F_occ_act_ = F_all.block(nfzc, nfzc, nocc_act, nocc_act);

    // Select just the unoccupied portion of the Fock matrix
    Eigen::MatrixXd F_uocc = F_all.block(nocc, nocc, nvir, nvir);

    // Compute all K_aibj
    auto K = fac.compute(L"(a i|G|b j)");
    const auto ktrange = K.trange();

    // zero out amplitudes
    if (!T_.is_initialized()) {
      T_ = Array(world, K.trange(), K.shape());
      T_.fill(0.0);
    }

    // For storing PNOs and and the Fock matrix in the PNO basis
    pnos_.resize(nocc_act * nocc_act);
    F_pno_diag_.resize(nocc_act * nocc_act);

    // For storing OSVs (PNOs when i = j) and the Fock matrix in
    // the OSV basis
    osvs_.resize(nocc_act);
    F_osv_diag_.resize(nocc_act);

    // Loop over each pair of occupieds to form PNOs
    for (int i = 0; i < nocc_act; ++i) {
      double eps_i = eps_o[i];

      for (int j = 0; j < nocc_act; ++j) {
        double eps_j = eps_o[j];
        int delta_ij = (i == j) ? 1 : 0;
        std::array<int, 4> tile_ij = {{0, i, 0, j}};
        std::array<int, 4> tile_ji = {{0, j, 0, i}};
        const auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
        const auto ord_ji = ktrange.tiles_range().ordinal(tile_ji);
        TA::TensorD K_ij = K.find(ord_ij);
        TA::TensorD K_ji = K.find(ord_ji);
        auto ext_ij = K_ij.range().extent_data();
        auto ext_ji = K_ji.range().extent_data();
        Eigen::MatrixXd K_ij_mat =
            TA::eigen_map(K_ij, ext_ij[0] * ext_ij[1], ext_ij[2] * ext_ij[3]);
        Eigen::MatrixXd K_ji_mat =
            TA::eigen_map(K_ji, ext_ji[0] * ext_ji[1], ext_ji[2] * ext_ji[3]);

        Eigen::MatrixXd T_ij(nvir, nvir);
        Eigen::MatrixXd T_ji(nvir, nvir);
        Eigen::MatrixXd T_tilde_ij(nvir, nvir);

        for (int a = 0; a < nvir; ++a) {
          double eps_a = eps_v[a];
          for (int b = 0; b < nvir; ++b) {
            double eps_b = eps_v[b];

            T_ij(a, b) = -K_ij_mat(a, b) / (eps_a + eps_b - eps_i - eps_j);
            T_ji(a, b) = -K_ji_mat(a, b) / (eps_a + eps_b - eps_i - eps_j);
          }
        }
        T_tilde_ij = 2 * T_ij - T_ji;
        Eigen::MatrixXd D_ij =
            (T_tilde_ij.adjoint() * T_ij + T_tilde_ij * T_ij.adjoint()) /
            (1 + delta_ij);

        // Diagonalize D_ij to get PNOs and corresponding occupation numbers.
        es.compute(D_ij);
        Eigen::MatrixXd pno_ij = es.eigenvectors();
        auto occ_ij = es.eigenvalues();

        // truncate PNOs
        size_t pnodrop = 0;
        for (size_t i = 0; i != occ_ij.cols(); ++i) {
          if (!(occ_ij(i) >= tpno_))
            ++pnodrop;
          else
            break;
        }
        const auto npno = nvir - pnodrop;

        // Store truncated PNOs
        // pnos[i*nocc_act + j] = pno_ij.block(0,pnodrop,nvir,npno);
        Eigen::MatrixXd pno_trunc = pno_ij.block(0, pnodrop, nvir, npno);
        pnos_[i * nocc_act + j] = pno_trunc;

        // Transform the Fock matrix to the PNO basis and store

        // Originally was transforming F_uocc using untruncated PNO matrix:
        // F_pno[i*nocc_act + j] = pno_ij.transpose() * F_uocc * pno_ij;

        // Now transforming using truncated PNO matrix; store just diag
        // elements:
        F_pno_diag_[i * nocc_act + j] =
            (pno_trunc.transpose() + F_uocc * pno_trunc).diagonal();

        auto nosv = 0;

        // truncate OSVs
        if (i == j) {
          size_t osvdrop = 0;
          for (size_t i = 0; i != occ_ij.cols(); ++i) {
            if (!(occ_ij(i) >= tosv_))
              ++osvdrop;
            else
              break;
          }
          nosv = nvir - osvdrop;

          // Store truncated OSVs
          // osvs[i] = pno_ij.block(0, osvdrop, nvir, nosv);
          Eigen::MatrixXd osv_trunc = pno_ij.block(0, osvdrop, nvir, nosv);
          osvs_[i] = osv_trunc;

          // Store Fock matrix transformed to diagonal OSV basis
          F_osv_diag_[i] =
              (osv_trunc.transpose() * F_uocc * osv_trunc).diagonal();
        }
      }
    }
  }
  virtual ~PNOSolver() = default;

  /// @return PNO truncation threshold
  double tpno() const { return tpno_; }
  /// @return OSV truncation threshold
  double tosv() const { return tosv_; }

 private:

  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space
  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    t1("a,i") += jacobi_update_t1(r1, F_occ_act_, F_osv_diag_, osvs_)("a,i");
    t2("a,b,i,j") += jacobi_update_t2(r2, F_occ_act_, F_pno_diag_, pnos_)("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }

  /// Overrides DIISSolver::update() .
  /// @note call DIISSolver::update(), then backtransform updated PNO amplitudes to the full space
  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    ::mpqc::cc::DIISSolver<T, T>::update(t1, t2, r1, r2);

    // T r1_copy = r1;
    // T r2_copy = r2;
    // T1T2<T, T> r(r1_copy, r2_copy);
    // T1T2<T, T> t(t1, t2);
    // diis_.extrapolate(t, r);
    // t1 = t.t1();
    // t2 = t.t2();


    // back-transform extrapolated t1 and t2 from OSV and PNO to full virtual
    // space
    ////////////////////////////////// again, this must be done by applying
    /// lambdas to delta_t1_ai and delta_t2_abij one tile at a time
    //assert(false && "not yet implemented");

  }

  void back_transform(T& t1, T& t2) {
    t1 = back_transform_t1(t1, osvs_);
    t2 = back_transform_t2(t2, pnos_);
  }


  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t2(
      const TA::DistArray<Tile, Policy>& r2_abij, const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_pno_diag,
      const std::vector<Eigen::MatrixXd>& pnos) {

    auto nocc_act = F_occ_act.rows();

    auto update2 = [F_occ_act, F_pno_diag, pnos, nocc_act](Tile& result_tile,
                                                       const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i and j indices
      const auto i = result_tile.range().lobound()[1];
      const auto j = result_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      Eigen::MatrixXd pno_ij = pnos[i * nocc_act + j];

      // Extent data of tile
      const auto ext = result_tile.range().extent_data();

      // Convert data in tile to Eigen Matrix
      // const Eigen::MatrixXd r2 = TA::eigen_map(result_tile, ext[2]*ext[3],
      // ext[0]*ext[1]);

      // Convert data in tile to Eigen::Map and transform to PNO basis
      const Eigen::MatrixXd r2_pno =
          pno_ij.transpose() *
          TA::eigen_map(result_tile, ext[0] * ext[1], ext[2] * ext[3]) * pno_ij;

      // Select correct vector containing diagonal elements of Fock matrix in
      // PNO basis
      const Eigen::VectorXd ens_uocc = F_pno_diag[i * nocc_act + j];

      // Determine number of PNOs
      const auto npno = ens_uocc.rows();

      // Select e_i and e_j
      const auto e_i = F_occ_act(i, i);
      const auto e_j = F_occ_act(j, j);

      auto tile_idx = 0;
      typename Tile::scalar_type norm = 0.0;
      for (auto a = 0; a < npno; ++a) {
        const auto e_a = ens_uocc[a];
        for (auto b = 0; b < npno; ++b, ++tile_idx) {
          const auto e_b = ens_uocc[b];
          const auto e_iajb = e_i + e_j - e_a - e_b;
          const auto old = arg_tile[tile_idx];
          const auto result_abij = old / e_iajb;
          const auto abs_result_abij = std::abs(result_abij);
          norm += abs_result_abij * abs_result_abij;
          result_tile[tile_idx] = result_abij;
        }
      }

      return std::sqrt(norm);
    };

    auto delta_t2_abij = TA::foreach (r2_abij, update2);
    delta_t2_abij.world().gop.fence();
    return delta_t2_abij;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t1(
      const TA::DistArray<Tile, Policy>& r1_ai, const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_osv_diag,
      const std::vector<Eigen::MatrixXd>& osvs) {
    auto update1 = [F_occ_act, F_osv_diag, osvs](Tile& result_tile,
                                             const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i index
      const auto i = result_tile.range().lobound()[1];

      // Select appropriate matrix of OSVs
      Eigen::MatrixXd osv_i = osvs[i];

      // Extent data of tile
      const auto ext = result_tile.range().extent_data();

      // Convert data in tile to Eigen Matrix
      // const Eigen::MatrixXd r1 = TA::eigen_map(result_tile, ext[1], ext[0]);

      // Convert data in tile to Eigen::Map and transform to OSV basis
      const Eigen::MatrixXd r1_osv =
          osv_i.transpose() * TA::eigen_map(result_tile, ext[1], ext[0]) *
          osv_i;

      // Select correct vector containing diagonal elements of Fock matrix in
      // OSV basis
      const Eigen::VectorXd ens_uocc = F_osv_diag[i];

      // Determine number of OSVs
      const auto nosv = ens_uocc.rows();

      // Select e_i
      const auto e_i = F_occ_act(i, i);

      auto tile_idx = 0;
      typename Tile::scalar_type norm = 0.0;
      for (auto a = 0; a < nosv; ++a) {
        const auto e_a = ens_uocc[a];
        const auto e_ia = e_i - e_a;
        const auto old = arg_tile[tile_idx];
        const auto result_ai = old / e_ia;
        const auto abs_result_ai = std::abs(result_ai);
        norm += abs_result_ai * abs_result_ai;
        result_tile[tile_idx] = result_ai;
      }
      return std::sqrt(norm);
    };

    auto delta_t1_ai = TA::foreach (r1_ai, update1);
    delta_t1_ai.world().gop.fence();
    return delta_t1_ai;
  }


  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> back_transform_t1(
      const TA::DistArray<Tile, Policy>& t1,
      const std::vector<Eigen::MatrixXd>& osvs) {
    auto back1 = [osvs](Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      const auto ext = result_tile.range().extent_data();

      // determine i index
      const auto i = result_tile.range().lobound()[1];


      // back transform
      const Eigen::MatrixXd result_tile_tr =
          osvs[i] * TA::eigen_map(result_tile, ext[1], ext[0]) *
          osvs[i].transpose();


      const auto rows = result_tile_tr.rows();
      const auto cols = result_tile_tr.cols();

      // Replace the original tile with the back-transformed tile
      for (auto r=0; r<rows; ++r) {
        for (auto c=0; c<cols; ++c) {
          auto idx = r*cols + c;
          result_tile[idx] = result_tile_tr(r,c);
        }
      }

      return 0;

    };

    auto transformed_t1 = TA::foreach (t1, back1);
    transformed_t1.world().gop.fence();
    return transformed_t1;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> back_transform_t2(
      const TA::DistArray<Tile, Policy>& t2,
      const std::vector<Eigen::MatrixXd>& pnos) {
    auto back2 = [pnos](Tile& result_tile, const Tile& arg_tile) {

      const auto nocc_act = (pnos.size())/2;

      result_tile = Tile(arg_tile.range());

      const auto ext = result_tile.range().extent_data();

      // determine i and j indices
      const auto i = result_tile.range().lobound()[1];
      const auto j = result_tile.range().lobound()[3];


      // back transform
      const Eigen::MatrixXd result_tile_tr =
          pnos[i*nocc_act + j] * TA::eigen_map(result_tile, ext[0] * ext[1], ext[2] * ext[3]) * 
          (pnos[i*nocc_act + j]).transpose();


      const auto rows = result_tile_tr.rows();
      const auto cols = result_tile_tr.cols();

      // Replace the original tile with the back-transformed tile
      for (auto r=0; r<rows; ++r) {
        for (auto c=0; c<cols; ++c) {
          auto idx = r*cols + c;
          result_tile[idx] = result_tile_tr(r,c);
        }
      }

      return 0;

    };

    auto transformed_t2 = TA::foreach (t2, back2);
    transformed_t2.world().gop.fence();
    return transformed_t2;
  }


  Factory<T>& factory_;
  std::string pno_method_;  //!< the PNO construction method
  double tpno_;             //!< the PNO truncation threshold
  double tosv_;             //!< the OSV (diagonal PNO) truncation threshold
  Array T_;

  Eigen::MatrixXd F_occ_act_;
  // For storing PNOs and and the Fock matrix in the PNO basis
  std::vector<Eigen::MatrixXd> pnos_;
  std::vector<Eigen::VectorXd> F_pno_diag_;

  // For storing OSVs (PNOs when i = j) and the Fock matrix in
  // the OSV basis
  std::vector<Eigen::MatrixXd> osvs_;
  std::vector<Eigen::VectorXd> F_osv_diag_;


};

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_ */
