#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_

// set to 1, must have libint-2.4.0-beta.2
#define PRODUCE_PNO_MOLDEN_FILES 1
#if PRODUCE_PNO_MOLDEN_FILES
#include "libint2/lcao/molden.h"
#endif

#include "mpqc/chemistry/molecule/common.h"
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
class PNOSolver : public ::mpqc::cc::DIISSolver<T, T>,
                  public madness::WorldObject<PNOSolver<T>> {
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
        madness::WorldObject<PNOSolver<T>>(factory.world()),
        factory_(factory),
        pno_method_(kv.value<std::string>("pno_method", "standard")),
        pno_canonical_(kv.value<std::string>("pno_canonical", "false")),
        tpno_(kv.value<double>("tpno", 1.e-8)),
        tosv_(kv.value<double>("tosv", 1.e-9)) {
    // part of WorldObject initialization
    this->process_pending();

    // compute and store PNOs truncated with threshold tpno_
    // store PNOs for diagonal pair as OSVs truncated with threshold tosv_

    // Check that tiling is done appropriately
    if (kv.exists("occ_block_size")) {
      int occ_block_size_ = (kv.value<int>("occ_block_size"));
      if (occ_block_size_ != 1) {
        throw InputError("occ_block_size must be set to 1 in the input file.");
      }
    } else {
      throw InputError("occ_block_size was not specified in the input file.");
    }

    if (kv.exists("unocc_block_size")) {
      int unocc_block_size_ = (kv.value<int>("unocc_block_size"));
      if (unocc_block_size_ < 1000000000) {
        throw InputError(
            "unocc_block_size must be greater than or equal to 1000000000 in "
            "the input file.");
      }
    } else {
      throw InputError("unocc_block_size was not specified in the input file.");
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

    auto& fac = factory_;
    auto& world = fac.world();
    auto& ofac = fac.orbital_registry();

    auto nocc = ofac.retrieve("m").rank();
    auto nocc_act = ofac.retrieve("i").rank();
    nocc_act_ = nocc_act;
    auto nvir = ofac.retrieve("a").rank();
    auto nfzc = nocc - nocc_act;

    // Form Fock array
    auto F = fac.compute(L"<p|F|q>[df]");

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
    // auto K = fac.compute(L"(a b|G|i j)");
    auto K = fac.compute(L"<a b|G|i j>[df]");
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

#if PRODUCE_PNO_MOLDEN_FILES
    // prepare to Molden
    const auto libint2_atoms = to_libint_atom(fac.atoms()->atoms());
    const auto C_i_eig = TA::array_to_eigen(ofac.retrieve("i").coefs());
    const auto C_a_eig = TA::array_to_eigen(ofac.retrieve("a").coefs());
    const auto libint2_shells = fac.basis_registry()->retrieve(L"μ")->flattened_shells();

    // write out active occupied orbitals
    auto occs = Eigen::VectorXd::Constant(C_i_eig.cols(), 2.0);
    libint2::molden::Export xport(libint2_atoms, libint2_shells, C_i_eig, occs);
    xport.write("occ.molden");
#endif

    // Loop over each pair of occupieds to form PNOs
    for (int i = 0; i < nocc_act; ++i) {
      double eps_i = eps_o[i];

      for (int j = 0; j < nocc_act; ++j) {
        double eps_j = eps_o[j];
        int delta_ij = (i == j) ? 1 : 0;
        std::array<int, 4> tile_ij = {{0, 0, i, j}};
        std::array<int, 4> tile_ji = {{0, 0, j, i}};
        const auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
        const auto ord_ji = ktrange.tiles_range().ordinal(tile_ji);
        TA::TensorD K_ij = K.find(ord_ij);
        TA::TensorD K_ji = K.find(ord_ji);
        auto ext_ij = K_ij.range().extent_data();
        auto ext_ji = K_ji.range().extent_data();
        Eigen::MatrixXd K_ij_mat =
            TA::eigen_map(K_ij, ext_ij[0] * ext_ij[2], ext_ij[1] * ext_ij[3]);
        Eigen::MatrixXd K_ji_mat =
            TA::eigen_map(K_ji, ext_ji[0] * ext_ji[2], ext_ji[1] * ext_ji[3]);

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

        // Eq. 23, JCP 128 034106 (2013)
        T_tilde_ij = 4 * T_ij - 2 * T_ji;
        Eigen::MatrixXd D_ij =
            (T_tilde_ij.transpose() * T_ij + T_tilde_ij * T_ij.transpose()) /
            (1.0 + delta_ij);

        // Diagonalize D_ij to get PNOs and corresponding occupation numbers.
        es.compute(D_ij);
        Eigen::MatrixXd pno_ij = es.eigenvectors();
        auto occ_ij = es.eigenvalues();

        // truncate PNOs
        size_t pnodrop = 0;
        if (tpno_ != 0.0) {
          for (size_t k = 0; k != occ_ij.rows(); ++k) {
            if (!(occ_ij(k) >= tpno_))
              ++pnodrop;
            else
              break;
          }
        }
        const auto npno = nvir - pnodrop;

        // Store truncated PNOs
        // pnos[i*nocc_act + j] = pno_ij.block(0,pnodrop,nvir,npno);
        Eigen::MatrixXd pno_trunc = pno_ij.block(0, pnodrop, nvir, npno);
        // pnos_[ij] = pno_trunc;
        pnos_[i * nocc_act + j] = pno_trunc;

#if PRODUCE_PNO_MOLDEN_FILES
        // write PNOs to Molden
        {
          Eigen::MatrixXd molden_coefs(C_i_eig.rows(), 2 + pno_trunc.cols());
          molden_coefs.col(0) = C_i_eig.col(i);
          molden_coefs.col(1) = C_i_eig.col(j);
          molden_coefs.block(0, 2, C_i_eig.rows(), pno_trunc.cols()) = C_a_eig * pno_trunc;

          Eigen::VectorXd occs(2 + pno_trunc.cols());
          occs.setZero();
          occs[0] = 2.0;
          occs[1] = 2.0;

          Eigen::VectorXd evals(2 + pno_trunc.cols());
          evals(0) = 0.0;
          evals(1) = 0.0;
          evals.tail(pno_trunc.cols()) = occ_ij.tail(pno_trunc.cols());

          libint2::molden::Export xport(libint2_atoms, libint2_shells,
                                        molden_coefs, occs, evals);
          xport.write(std::string("pno_") + std::to_string(i) + "_" +
                      std::to_string(j) + ".molden");
        }
#endif

        // Transform F to PNO space
        Eigen::MatrixXd F_pno_ij = pno_trunc.transpose() * F_uocc * pno_trunc;

        // Store just the diagonal elements of F_pno_ij
        // F_pno_diag_[ij] = F_pno_ij.diagonal();
        F_pno_diag_[i * nocc_act + j] = F_pno_ij.diagonal();

        /////// Transform PNOs to canonical PNOs if pno_canonical_ == true

        if (pno_canonical_ == "true" && npno > 0) {
          // Compute eigenvectors of F in PNO space
          es.compute(F_pno_ij);
          Eigen::MatrixXd pno_transform_ij = es.eigenvectors();

          // Transform pno_ij to canonical PNO space; pno_ij -> can_pno_ij
          Eigen::MatrixXd can_pno_ij = pno_trunc * pno_transform_ij;

          // Replace standard with canonical PNOs
          // pnos_[ij] = can_pno_ij;
          // F_pno_diag_[ij] = es.eigenvalues();
          pnos_[i * nocc_act + j] = can_pno_ij;
          F_pno_diag_[i * nocc_act + j] = es.eigenvalues();
        }

        // truncate OSVs

        // auto nosv = 0;

        // auto osvdrop = 0;
        if (i == j) {
          size_t osvdrop = 0;
          if (tosv_ != 0.0) {
            for (size_t k = 0; k != occ_ij.rows(); ++k) {
              if (!(occ_ij(k) >= tosv_))
                ++osvdrop;
              else
                break;
            }
          }
          const auto nosv = nvir - osvdrop;
          if (nosv == 0) {  // all OSV truncated indicates total nonsense
            throw LimitExceeded<size_t>("all OSVs truncated", __FILE__,
                                        __LINE__, 1, 0);
          }

          // Store truncated OSVs
          Eigen::MatrixXd osv_trunc = pno_ij.block(0, osvdrop, nvir, nosv);
          osvs_[i] = osv_trunc;

          // Transform F to OSV space
          Eigen::MatrixXd F_osv_i = osv_trunc.transpose() * F_uocc * osv_trunc;

          // Store just the diagonal elements of F_osv_i
          F_osv_diag_[i] = F_osv_i.diagonal();

          /////// Transform OSVs to canonical OSVs if pno_canonical_ == true
          if (pno_canonical_ == "true") {
            // Compute eigenvectors of F in OSV space
            es.compute(F_osv_i);
            Eigen::MatrixXd osv_transform_i = es.eigenvectors();

            // Transform osv_i to canonical OSV space: osv_i -> can_osv_i
            Eigen::MatrixXd can_osv_i = osv_trunc * osv_transform_i;

            // Replace standard with canonical OSVs
            osvs_[i] = can_osv_i;
            F_osv_diag_[i] = es.eigenvalues();
          }
        }
      }
    }
    auto sum_osv = 0;
    for (int i = 0; i < nocc_act; ++i) {
      sum_osv += osvs_[i].cols();
    }
    auto ave_nosv = sum_osv / nocc_act;
    ExEnv::out0() << "The average number of OSVs is " << ave_nosv << std::endl;

    auto sum_pno = 0;
    for (int i = 0; i < nocc_act; ++i) {
      for (int j = 0; j < nocc_act; ++j) {
        sum_pno += pnos_[i * nocc_act + j].cols();
      }
    }
    auto ave_npno = sum_pno / (nocc_act * nocc_act);
    ExEnv::out0() << "The average number of PNOs is " << ave_npno << std::endl;
  }
  virtual ~PNOSolver() = default;

  /// @return PNO truncation threshold
  double tpno() const { return tpno_; }
  /// @return OSV truncation threshold
  double tosv() const { return tosv_; }

  const auto& pno(int i, int j) const { return pnos_[i * nocc_act_ + j]; }
  const auto& osv(int i) const { return osvs_[i]; }

 private:
  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space
  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    auto delta_t1_ai = jacobi_update_t1(r1, F_occ_act_, F_osv_diag_, osvs_);
    auto delta_t2_abij = jacobi_update_t2(r2, F_occ_act_, F_pno_diag_, pnos_);
    t1("a,i") += delta_t1_ai("a,i");
    t2("a,b,i,j") += delta_t2_abij("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }

  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    update_only(t1, t2, r1, r2);
    T r1_osv = osv_transform_ai(r1, osvs_);
    T r2_pno = pno_transform_abij(r2, pnos_);
    mpqc::cc::T1T2<T, T> r(r1_osv, r2_pno);
    mpqc::cc::T1T2<T, T> t(t1, t2);
    this->diis().extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t2(
      const TA::DistArray<Tile, Policy>& r2_abij,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_pno_diag,
      const std::vector<Eigen::MatrixXd>& pnos) {

    auto update2 = [F_occ_act, F_pno_diag, pnos, this](
                       Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd pno_ij = pnos[ij];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to PNO basis
      const Eigen::MatrixXd r2_pno =
          pno_ij.transpose() *
          TA::eigen_map(arg_tile, ext[0] * ext[2], ext[1] * ext[3]) * pno_ij;

      // Create a matrix delta_t2_pno to hold updated values of delta_t2 in PNO
      // basis this matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::MatrixXd delta_t2_pno = r2_pno;

      // Select correct vector containing diagonal elements of Fock matrix in
      // PNO basis
      const Eigen::VectorXd& ens_uocc = F_pno_diag[ij];

      // Determine number of PNOs
      const auto npno = ens_uocc.rows();

      // Determine number of uocc
      const auto nuocc = pno_ij.rows();

      // Select e_i and e_j
      const auto e_i = F_occ_act(i, i);
      const auto e_j = F_occ_act(j, j);

      for (auto a = 0; a < npno; ++a) {
        const auto e_a = ens_uocc[a];
        for (auto b = 0; b < npno; ++b) {
          const auto e_b = ens_uocc[b];
          const auto e_abij = e_i + e_j - e_a - e_b;
          const auto r_abij = r2_pno(a, b);
          delta_t2_pno(a, b) = r_abij / e_abij;
        }
      }

      // Back transform delta_t2_pno to full space
      Eigen::MatrixXd delta_t2_full =
          pno_ij * delta_t2_pno * pno_ij.transpose();

      // Convert delta_t2_full to tile and compute norm
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nuocc; ++r) {
        for (auto c = 0; c < nuocc; ++c) {
          const auto idx = r * nuocc + c;
          const auto elem = delta_t2_full(r, c);
          const auto abs_elem = std::abs(elem);
          norm += abs_elem * abs_elem;
          result_tile[idx] = elem;
        }
      }

      return std::sqrt(norm);
    };

    auto delta_t2_abij = TA::foreach(r2_abij, update2);
    delta_t2_abij.world().gop.fence();
    return delta_t2_abij;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t1(
      const TA::DistArray<Tile, Policy>& r1_ai,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_osv_diag,
      const std::vector<Eigen::MatrixXd>& osvs) {
    auto update1 = [F_occ_act, F_osv_diag, osvs](Tile& result_tile,
                                                 const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i index
      const auto i = arg_tile.range().lobound()[1];

      // Select appropriate matrix of OSVs
      Eigen::MatrixXd osv_i = osvs[i];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to OSV basis
      const Eigen::VectorXd r1_osv =
          osv_i.transpose() * TA::eigen_map(arg_tile, ext[0], ext[1]);

      // Create a matrix delta_t1_osv to hold updated values of delta t1 in OSV
      // basis this matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::VectorXd delta_t1_osv = r1_osv;

      // Select correct vector containing diagonal elements of Fock matrix in
      // OSV basis
      const Eigen::VectorXd& ens_uocc = F_osv_diag[i];

      // Determine number of OSVs
      const auto nosv = ens_uocc.rows();

      // Determine number of uocc
      const auto nuocc = osv_i.rows();

      // Select e_i
      const auto e_i = F_occ_act(i, i);

      for (auto a = 0; a < nosv; ++a) {
        const auto e_a = ens_uocc[a];
        const auto e_ai = e_i - e_a;
        const auto r_ai = r1_osv(a);
        delta_t1_osv(a) = r_ai / e_ai;
      }

      // Back transform delta_t1_osv to full space
      // Eigen::MatrixXd delta_t1_full = osv_i * delta_t1_osv *
      // osv_i.transpose();
      Eigen::VectorXd delta_t1_full = osv_i * delta_t1_osv;

      // Convert delta_t1_full to tile and compute norm
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nuocc; ++r) {
        const auto elem = delta_t1_full(r);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto delta_t1_ai = TA::foreach (r1_ai, update1);
    delta_t1_ai.world().gop.fence();
    return delta_t1_ai;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> pno_transform_abij(
      const TA::DistArray<Tile, Policy>& abij,
      const std::vector<Eigen::MatrixXd>& pnos) {

    auto tform = [pnos, this](
        Tile& result_tile, const Tile& arg_tile) {

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      const auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd pno_ij = pnos[ij];
      const auto nuocc = pno_ij.rows();
      const auto npno = pno_ij.cols();

      // Convert data in tile to Eigen::Map and transform to PNO basis
      const Eigen::MatrixXd result_eig =
          pno_ij.transpose() * TA::eigen_map(arg_tile, nuocc, nuocc) * pno_ij;

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{npno,npno,1l,1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < npno; ++r) {
        for (auto c = 0; c < npno; ++c) {
          const auto idx = r * npno + c;
          const auto elem = result_eig(r, c);
          const auto abs_elem = std::abs(elem);
          norm += abs_elem * abs_elem;
          result_tile[idx] = elem;
        }
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach(abij, tform);
    result.world().gop.fence();
    return result;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> osv_transform_ai(
      const TA::DistArray<Tile, Policy>& ai,
      const std::vector<Eigen::MatrixXd>& osvs) {

    auto tform = [osvs, this](
        Tile& result_tile, const Tile& arg_tile) {

      // determine i index
      const auto i = arg_tile.range().lobound()[1];

      // Select appropriate matrix of OSVs
      Eigen::MatrixXd osv_i = osvs[i];
      const auto nuocc = osv_i.rows();
      const auto nosv = osv_i.cols();

      // Convert data in tile to Eigen::Map and transform to OSV basis
      const Eigen::MatrixXd result_eig =
          osv_i.transpose() * TA::eigen_map(arg_tile, nuocc, 1);

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{nosv,1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nosv; ++r) {
        const auto elem = result_eig(r, 0);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach(ai, tform);
    result.world().gop.fence();
    return result;
  }

  // squared norm of 1-body residual in OSV subspace
  struct R1SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R1SquaredNormReductionOp(PNOSolver<T>* solver) : solver_(solver) {}

    // Reduction functions
    // Make an empty result object
    result_type operator()() const { return 0; }
    // Post process the result (no operation, passthrough)
    const result_type& operator()(const result_type& result) const {
      return result;
    }
    void operator()(result_type& result, const result_type& arg) const {
      result += arg;
    }
    /// Reduce an argument pair
    void operator()(result_type& result, const argument_type& arg) const {
      const auto i = arg.range().lobound()[1];
      const auto nuocc = arg.range().extent_data()[0];
      const Eigen::MatrixXd arg_osv =
          TA::eigen_map(arg, 1, nuocc) * solver_->osv(i);
      result += arg_osv.squaredNorm();
    }

    PNOSolver<T>* solver_;
  };  // R1SquaredNormReductionOp

  // squared norm of 2-body residual in PNO subspace
  struct R2SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R2SquaredNormReductionOp(PNOSolver<T>* solver) : solver_(solver) {}

    // Reduction functions
    // Make an empty result object
    result_type operator()() const { return 0; }
    // Post process the result (no operation, passthrough)
    const result_type& operator()(const result_type& result) const {
      return result;
    }
    void operator()(result_type& result, const result_type& arg) const {
      result += arg;
    }
    /// Reduce an argument pair
    void operator()(result_type& result, const argument_type& arg) const {
      const auto i = arg.range().lobound()[2];
      const auto j = arg.range().lobound()[3];
      const auto nuocc = arg.range().extent_data()[0];
      const Eigen::MatrixXd arg_pno = solver_->pno(i, j).transpose() *
                                      TA::eigen_map(arg, nuocc, nuocc) *
                                      solver_->pno(i, j);
      result += arg_pno.squaredNorm();
    }

    PNOSolver<T>* solver_;
  };  // R2SquaredNormReductionOp

 public:
  /// Overrides Solver<T,T>::error()
  virtual double error(const T& r1, const T& r2) override {
    R1SquaredNormReductionOp op1(this);
    R2SquaredNormReductionOp op2(this);
    return sqrt(r1("a,i").reduce(op1).get() + r2("a,b,i,j").reduce(op2).get()) /
           (size(r1) + size(r2));
  }

 private:
  Factory<T>& factory_;
  std::string pno_method_;     //!< the PNO construction method
  std::string pno_canonical_;  //!< whether or not to canonicalize PNO/OSV
  double tpno_;                //!< the PNO truncation threshold
  double tosv_;                //!< the OSV (diagonal PNO) truncation threshold
  int nocc_act_;               //!< the number of active occupied orbitals
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


////// PSVO solver ////////

/// PSVOSolver updates the CC T amplitudes using standard Jacobi+DIIS in PSVO
/// space
/// @warning This class assumes that the 1- and 2-body amplitudes/residuals
///          given to Solver::update() are laid out as "a,i" and "a,b,i,j",
///          respectively
template <typename T>
class PSVOSolver : public ::mpqc::cc::DIISSolver<T, T>,
                  public madness::WorldObject<PSVOSolver<T>> {
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
  PSVOSolver(const KeyVal& kv, Factory<T>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv),
        madness::WorldObject<PSVOSolver<T>>(factory.world()),
        factory_(factory),
        tpsvo_(kv.value<double>("tpsvo", 1.e-5)) {
    // part of WorldObject initialization
    this->process_pending();

    // Check that tiling is done appropriately
    if (kv.exists("occ_block_size")) {
      int occ_block_size_ = (kv.value<int>("occ_block_size"));
      if (occ_block_size_ != 1) {
        throw InputError("occ_block_size must be set to 1 in the input file.");
      }
    } 
    else {
      throw InputError("occ_block_size was not specified in the input file.");
    }

    if (kv.exists("unocc_block_size")) {
      int unocc_block_size_ = (kv.value<int>("unocc_block_size"));
      if (unocc_block_size_ < 1000000000) {
        throw InputError(
            "unocc_block_size must be greater than or equal to 1000000000 in "
            "the input file.");
      }
    } 
    else {
      throw InputError("unocc_block_size was not specified in the input file.");
    }

    // Compute integrals

    auto& fac = factory_;
    auto& world = fac.world();
    auto& ofac = fac.orbital_registry();

    auto nocc = ofac.retrieve("m").rank();
    auto nocc_act = ofac.retrieve("i").rank();
    nocc_act_ = nocc_act;
    auto nvir = ofac.retrieve("a").rank();
    auto nfzc = nocc - nocc_act;

    // Form Fock array
    auto F = fac.compute(L"<p|F|q>[df]");

    // Select just diagonal elements of Fock aray and transform
    // to Eigen vector; use for computing PSVOs
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
    // auto K = fac.compute(L"(a b|G|i j)");
    auto K = fac.compute(L"<a b|G|i j>[df]");
    const auto ktrange = K.trange();

    // zero out amplitudes
    if (!T_.is_initialized()) {
      T_ = Array(world, K.trange(), K.shape());
      T_.fill(0.0);
    }


    // For storing PSVOs and and the Fock matrix in the PSVO basis
    l_psvos_.resize(nocc_act * nocc_act);
    r_psvos_.resize(nocc_act * nocc_act);
    F_l_psvo_diag_.resize(nocc_act * nocc_act);
    F_r_psvo_diag_.resize(nocc_act * nocc_act);

    // Loop over each pair of occupieds to form amplitude matrices
    for (int i = 0; i < nocc_act; ++i) {
      double eps_i = eps_o[i];

      for (int j = 0; j < nocc_act; ++j) {
        double eps_j = eps_o[j];
        int delta_ij = (i == j) ? 1 : 0;
        std::array<int, 4> tile_ij = {{0, 0, i, j}};
        std::array<int, 4> tile_ji = {{0, 0, j, i}};
        const auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
        const auto ord_ji = ktrange.tiles_range().ordinal(tile_ji);
        TA::TensorD K_ij = K.find(ord_ij);
        TA::TensorD K_ji = K.find(ord_ji);
        auto ext_ij = K_ij.range().extent_data();
        auto ext_ji = K_ji.range().extent_data();
        Eigen::MatrixXd K_ij_mat =
            TA::eigen_map(K_ij, ext_ij[0] * ext_ij[2], ext_ij[1] * ext_ij[3]);
        Eigen::MatrixXd K_ji_mat =
            TA::eigen_map(K_ji, ext_ji[0] * ext_ji[2], ext_ji[1] * ext_ji[3]);

        Eigen::MatrixXd T_ij(nvir, nvir);

        for (int a = 0; a < nvir; ++a) {
          double eps_a = eps_v[a];
          for (int b = 0; b < nvir; ++b) {
            double eps_b = eps_v[b];

            T_ij(a, b) = -K_ij_mat(a, b) / (eps_a + eps_b - eps_i - eps_j);

          } // for each b
        } // for each a

        // Compute SVD of each T_ij matrix
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T_ij, Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto sing_vals = svd.singularValues();

        // truncate PSVOs
        size_t psvo_drop = 0;
        if (tpsvo_ != 0.0) {
          for (size_t k = 0; k != sing_vals.rows(); ++k) {
            if (!(sing_vals(k) >= tpsvo_))
              ++psvo_drop;
            else
              break;
          } // for each k
        } // if tpsvo != 0

        const auto npsvo = nvir - psvo_drop;

        // store truncated PSVOs
        Eigen::MatrixXd r_psvo_trunc = svd.matrixV().block(0, psvo_drop, nvir, npsvo);
        Eigen::MatrixXd l_psvo_trunc = svd.matrixU().block(0, psvo_drop, nvir, npsvo);
        r_psvos_[i*nocc_act + j] = r_psvo_trunc;
        l_psvos_[i*nocc_act + j] = l_psvo_trunc;

        // transform F to right PSVO space and store just diagonal elements
        Eigen::MatrixXd F_r_psvo = svd.matrixV().transpose() * F_uocc * svd.matrixV();
        F_r_psvo_diag_[i*nocc_act + j] = F_r_psvo.diagonal();

        // transform F to left PSVO space and store just diagonal elements
        Eigen::MatrixXd F_l_psvo = svd.matrixU().transpose() * F_uocc.transpose() * svd.matrixU();
        F_l_psvo_diag_[i*nocc_act + j] = F_l_psvo.diagonal();


      } // for each j
    } // for each i

    // Compute average number of PSVOs per pair and print out

    auto sum_psvo = 0;
    for (int i=0; i<nocc_act; ++i) {
      for (int j=0; j<nocc_act; ++j) {
        sum_psvo += r_psvos_[i*nocc_act + j].cols();
      }
    } 
    auto ave_npsvo = sum_psvo / (nocc_act * nocc_act);
    ExEnv::out0() << "The average number of PSVOs is " << ave_npsvo << std::endl;
  } // PSVOSolver

  virtual ~PSVOSolver() = default;

  /// @return PSVO truncation threshold
  double tpsvo() const { return tpsvo_; }

  const auto& l_psvo(int i, int j) const { return l_psvos_[i*nocc_act_ + j]; }
  const auto& r_psvo(int i, int j) const { return r_psvos_[i*nocc_act_ + j]; }

private:
  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space
  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    auto delta_t1_ai = jacobi_update_t1(r1, F_occ_act_, F_l_psvo_diag_, 
                                        F_r_psvo_diag_, l_psvos_, r_psvos_);

    auto delta_t2_abij = jacobi_update_t2(r2, F_occ_act_, F_l_psvo_diag_,
                                          F_r_psvo_diag_, l_psvos_, r_psvos_);
    t1("a,i") += delta_t1_ai("a,i");
    t2("a,b,i,j") += delta_t2_abij("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }

  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    update_only(t1, t2, r1, r2);
    T r1_psvo = psvo_transform_ai(r1, l_psvos_, r_psvos_);
    T r2_psvo = psvo_transform_abij(r2, l_psvos_, r_psvos_);
    mpqc::cc::T1T2<T, T> r(r1_psvo, r2_psvo);
    mpqc::cc::T1T2<T, T> t(t1, t2);
    this->diis().extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t2(
      const TA::DistArray<Tile, Policy>& r2_abij,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_l_psvo_diag,
      const std::vector<Eigen::VectorXd>& F_r_psvo_diag,
      const std::vector<Eigen::MatrixXd>& l_psvos,
      const std::vector<Eigen::MatrixXd>& r_psvos) {

    auto update2 = [F_occ_act, F_l_psvo_diag, F_r_psvo_diag, l_psvos, r_psvos, this](
                       Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd l_psvo_ij = l_psvos[ij];
      Eigen::MatrixXd r_psvo_ij = r_psvos[ij];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to PNO basis
      const Eigen::MatrixXd r2_psvo =
          l_psvo_ij.adjoint() *
          TA::eigen_map(arg_tile, ext[0] * ext[2], ext[1] * ext[3]) * r_psvo_ij;

      // Create a matrix delta_t2_pno to hold updated values of delta_t2 in PNO
      // basis this matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::MatrixXd delta_t2_psvo = r2_psvo;

      // Select correct vector containing diagonal elements of Fock matrix in
      // PNO basis
      const Eigen::VectorXd& l_uocc = F_l_psvo_diag[ij];
      const Eigen::VectorXd& r_uocc = F_r_psvo_diag[ij];

      // Determine number of PNOs
      const auto npno = l_uocc.rows();

      // Determine number of uocc
      const auto nuocc = l_psvo_ij.rows();

      // Select e_i and e_j
      const auto e_i = F_occ_act(i, i);
      const auto e_j = F_occ_act(j, j);

      for (auto a = 0; a < npno; ++a) {
        const auto e_a = l_uocc[a];
        for (auto b = 0; b < npno; ++b) {
          const auto e_b = r_uocc[b];
          const auto e_abij = e_i + e_j - e_a - e_b;
          const auto r_abij = r2_psvo(a, b);
          delta_t2_psvo(a, b) = r_abij / e_abij;
        }
      }

      // Back transform delta_t2_psvo to full space
      Eigen::MatrixXd delta_t2_full =
          r_psvo_ij * delta_t2_psvo * l_psvo_ij.adjoint();

      // Convert delta_t2_full to tile and compute norm
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nuocc; ++r) {
        for (auto c = 0; c < nuocc; ++c) {
          const auto idx = r * nuocc + c;
          const auto elem = delta_t2_full(r, c);
          const auto abs_elem = std::abs(elem);
          norm += abs_elem * abs_elem;
          result_tile[idx] = elem;
        }
      }

      return std::sqrt(norm);
    };

    auto delta_t2_abij = TA::foreach(r2_abij, update2);
    delta_t2_abij.world().gop.fence();
    return delta_t2_abij;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t1(
      const TA::DistArray<Tile, Policy>& r1_ai,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_r_psvo_diag,
      const std::vector<Eigen::MatrixXd>& r_psvos) {
    auto update1 = [F_occ_act, F_r_psvo_diag, r_psvos, this](
                      Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i index
      const auto i = arg_tile.range().lobound()[1];

      // Select appropriate matrix of PSVOs
      Eigen::MatrixXd r_psvo_i = r_psvos[i*nocc_act_ + i];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to PSVO basis
      const Eigen::VectorXd r1_psvo =
          r_psvo_i.transpose() * TA::eigen_map(arg_tile, ext[0], ext[1]);

      // Create a matrix delta_t1_osv to hold updated values of delta t1 in OSV
      // basis this matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::VectorXd delta_t1_psvo = r1_psvo;

      // Select correct vector containing diagonal elements of Fock matrix in
      // PSVO basis
      const Eigen::VectorXd& r_uocc = F_r_psvo_diag[i*nocc_act_ + i];

      // Determine number of OSVs
      const auto nosv = r_uocc.rows();

      // Determine number of uocc
      const auto nuocc = r_psvo_i.rows();

      // Select e_i
      const auto e_i = F_occ_act(i, i);

      for (auto a = 0; a < nosv; ++a) {
        const auto e_a = r_uocc[a];
        const auto e_ai = e_i - e_a;
        const auto r_ai = r1_psvo(a);
        delta_t1_psvo(a) = r_ai / e_ai;
      }

      // Back transform delta_t1_osv to full space
      // Eigen::MatrixXd delta_t1_full = osv_i * delta_t1_osv *
      // osv_i.transpose();
      Eigen::VectorXd delta_t1_full = r_psvo_i * delta_t1_psvo;

      // Convert delta_t1_full to tile and compute norm
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nuocc; ++r) {
        const auto elem = delta_t1_full(r);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto delta_t1_ai = TA::foreach (r1_ai, update1);
    delta_t1_ai.world().gop.fence();
    return delta_t1_ai;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> psvo_transform_abij(
      const TA::DistArray<Tile, Policy>& abij,
      const std::vector<Eigen::MatrixXd>& l_psvos,
      const std::vector<Eigen::MatrixXd>& r_psvos) {

    auto tform = [l_psvos, r_psvos, this](
        Tile& result_tile, const Tile& arg_tile) {

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      const auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd l_psvo_ij = l_psvos[ij];
      Eigen::MatrixXd r_psvo_ij = r_psvos[ij];
      const auto nuocc = l_psvo_ij.rows();
      const auto npsvo = l_psvo_ij.cols();

      // Convert data in tile to Eigen::Map and transform to PNO basis
      const Eigen::MatrixXd result_eig =
          l_psvo_ij.transpose() * TA::eigen_map(arg_tile, nuocc, nuocc) * r_psvo_ij;

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{npsvo,npsvo,1l,1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < npsvo; ++r) {
        for (auto c = 0; c < npsvo; ++c) {
          const auto idx = r * npsvo + c;
          const auto elem = result_eig(r, c);
          const auto abs_elem = std::abs(elem);
          norm += abs_elem * abs_elem;
          result_tile[idx] = elem;
        }
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach(abij, tform);
    result.world().gop.fence();
    return result;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> psvo_transform_ai(
      const TA::DistArray<Tile, Policy>& ai,
      const std::vector<Eigen::MatrixXd>& r_psvos) {

    auto tform = [r_psvos, this](
        Tile& result_tile, const Tile& arg_tile) {

      // determine i index
      const auto i = arg_tile.range().lobound()[1];

      // Select appropriate matrix of OSVs
      Eigen::MatrixXd r_psvo_i = r_psvos[i];
      const auto nuocc = r_psvo_i.rows();
      const auto npsvo = r_psvo_i.cols();

      // Convert data in tile to Eigen::Map and transform to OSV basis
      const Eigen::MatrixXd result_eig =
          r_psvo_i.transpose() * TA::eigen_map(arg_tile, nuocc, 1);

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{npsvo,1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < npsvo; ++r) {
        const auto elem = result_eig(r, 0);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach(ai, tform);
    result.world().gop.fence();
    return result;
  }

  private:
  Factory<T>& factory_;
  double tpsvo_;        //!< the truncation threshold for PSVOs
  int nocc_act_;        //!< the number of active occupied orbitals
  Array T_;

  Eigen::MatrixXd F_occ_act_;

  // For storing PSVOs and and the Fock matrix in the PSVO basis
  std::vector<Eigen::MatrixXd> l_psvos_;
  std::vector<Eigen::MatrixXd> r_psvos_;
  std::vector<Eigen::VectorXd> F_l_psvo_diag_;
  std::vector<Eigen::VectorXd> F_r_psvo_diag_;




}; // class: PSVOSolver



}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_ */
