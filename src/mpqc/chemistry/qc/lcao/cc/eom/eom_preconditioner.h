//
// Created by Chong Peng on 8/9/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_PRECONDITIONER_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_PRECONDITIONER_H_

#include "mpqc/chemistry/qc/lcao/cc/solvers.h"
#include "mpqc/math/linalg/davidson_diag.h"

namespace mpqc {
namespace lcao {
namespace cc {

/// preconditioner in DavidsonDiag for EOM-EA-CCSD, approximate the diagonal
/// H_bar matrix
template <typename Array>
class EEPred : public DavidsonDiagPred<::mpqc::cc::TPack<Array>> {
 public:
  using element_type = typename Array::element_type;
  using Tile = typename Array::value_type;

  EEPred(const EigenVector<element_type> &eps_O,
         const EigenVector<element_type> &eps_V)
      : eps_o_(eps_O), eps_v_(eps_V) {}

  // default constructor
  EEPred() : eps_o_(), eps_v_() {}

  ~EEPred() = default;

  /// override the abstract virtual function
  virtual void operator()(
      const EigenVector<element_type> &e,
      std::vector<::mpqc::cc::TPack<Array>> &guess) const override {
    std::size_t n_roots = e.size();
    TA_ASSERT(n_roots == guess.size());
    for (std::size_t i = 0; i < n_roots; i++) {
      compute(e[i], guess[i]);
    }
  }

  void compute(const element_type &e, ::mpqc::cc::TPack<Array> &guess) const {
    const auto &eps_v = this->eps_v_;
    const auto &eps_o = this->eps_o_;

    auto task1 = [&eps_v, &eps_o, e](Tile &result_tile) {
      const auto &range = result_tile.range();
      float norm = 0.0;
      for (const auto &i : range) {
        const auto result = result_tile[i] / (e + eps_o[i[1]] - eps_v[i[0]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    auto task2 = [&eps_v, &eps_o, e](Tile &result_tile) {
      const auto &range = result_tile.range();
      float norm = 0.0;
      for (const auto &i : range) {
        const auto result = result_tile[i] / (e - eps_v[i[0]] - eps_v[i[1]] +
                                              eps_o[i[2]] + eps_o[i[3]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    TA::foreach_inplace(guess[0], task1);
    TA::foreach_inplace(guess[1], task2);

    guess[0].world().gop.fence();
  }

 private:
  /// diagonal of F_ij matrix
  EigenVector<element_type> eps_o_;
  /// diagonal of F_ab matrix
  EigenVector<element_type> eps_v_;
};

// preconditioner of EOM-EA-CCSD approximate the diagonal of H matrix
template <typename Array>
class EAPred : public DavidsonDiagPred<::mpqc::cc::TPack<Array>> {
 public:
  using element_type = typename Array::element_type;
  using Tile = typename Array::value_type;
  /// diagonal of F_ij matrix
  EigenVector<element_type> eps_o;
  /// diagonal of F_ab matrix
  EigenVector<element_type> eps_v;

  EAPred(const EigenVector<element_type> &eps_O,
         const EigenVector<element_type> &eps_V)
      : eps_o(eps_O), eps_v(eps_V) {}

  // default constructor
  EAPred() = default;
  ~EAPred() = default;

  /// override the abstract virtual function
  virtual void operator()(
      const EigenVector<element_type> &e,
      std::vector<::mpqc::cc::TPack<Array>> &guess) const override {
    std::size_t n_roots = e.size();
    TA_ASSERT(n_roots == guess.size());
    for (std::size_t i = 0; i < n_roots; i++) {
      compute(e[i], guess[i]);
    }
  }

  void compute(const element_type &e, ::mpqc::cc::TPack<Array> &guess) const {
    const auto &eps_v = this->eps_v;
    const auto &eps_o = this->eps_o;

    auto task1 = [&eps_v, e](Tile &result_tile) {
      const auto &range = result_tile.range();
      float norm = 0.0;
      for (const auto &i : range) {
        const auto result = result_tile[i] / (e - eps_v[i[0]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    auto task2 = [&eps_v, &eps_o, e](Tile &result_tile) {
      const auto &range = result_tile.range();
      float norm = 0.0;
      for (const auto &i : range) {
        const auto result =
            result_tile[i] / (e - eps_v[i[0]] - eps_v[i[1]] + eps_o[i[2]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    TA::foreach_inplace(guess[0], task1);
    TA::foreach_inplace(guess[1], task2);

    guess[0].world().gop.fence();
  }
};

/// preconditioner for EOM-IP-CCSD, approximate the diagonal of H matrix
template <typename Array>
class IPPred : public DavidsonDiagPred<::mpqc::cc::TPack<Array>> {
 public:
  using element_type = typename Array::element_type;
  using Tile = typename Array::value_type;
  /// diagonal of F_ij matrix
  EigenVector<element_type> eps_o;
  /// diagonal of F_ab matrix
  EigenVector<element_type> eps_v;

  IPPred(const EigenVector<element_type> &eps_O,
         const EigenVector<element_type> &eps_V)
      : eps_o(eps_O), eps_v(eps_V) {}

  // default constructor
  IPPred() = default;
  ~IPPred() = default;

  /// override the abstract virtual function
  virtual void operator()(
      const EigenVector<element_type> &e,
      std::vector<::mpqc::cc::TPack<Array>> &guess) const override {
    std::size_t n_roots = e.size();
    TA_ASSERT(n_roots == guess.size());
    for (std::size_t i = 0; i < n_roots; i++) {
      compute(e[i], guess[i]);
    }
  }

  void compute(const element_type &e, ::mpqc::cc::TPack<Array> &guess) const {
    const auto &eps_v = this->eps_v;
    const auto &eps_o = this->eps_o;

    auto task1 = [&eps_o, e](Tile &result_tile) {
      const auto &range = result_tile.range();
      float norm = 0.0;
      for (const auto &i : range) {
        const auto result = result_tile[i] / (e + eps_o[i[0]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    auto task2 = [&eps_v, &eps_o, e](Tile &result_tile) {
      const auto &range = result_tile.range();
      float norm = 0.0;
      for (const auto &i : range) {
        const auto result =
            result_tile[i] / (e - eps_v[i[0]] + eps_o[i[1]] + eps_o[i[2]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    TA::foreach_inplace(guess[0], task1);
    TA::foreach_inplace(guess[1], task2);

    guess[0].world().gop.fence();
  }
};

/// PNO preconditioner for EOM-EE-CCSD
template <typename Array>
class PNOEEPred : public DavidsonDiagPred<::mpqc::cc::TPack<Array>> {
 public:
  using Tile = typename Array::value_type;
  using Matrix = RowMatrix<typename Tile::numeric_type>;
  using Vector = EigenVector<typename Tile::numeric_type>;

  PNOEEPred() = default;

  PNOEEPred(const Array &T2, std::size_t n_roots, const Vector &eps_o,
            const Matrix &F_uocc, double tpno, double tosv, bool pno_canonical)
      : tpno_(tpno),
        tosv_(tosv),
        pno_canonical_(pno_canonical),
        n_roots_(n_roots),
        eps_o_(eps_o) {
    init_reblock(T2);

    // use first excited state amplitude to initialize PNOs
    auto T_reblock = detail::reblock_t2(T2, reblock_i_, reblock_a_);
    auto D = detail::construct_density(T_reblock);

    detail::construct_pno(D, F_uocc,
                          tpno_, tosv_,
                          pnos_, npnos_, F_pno_diag_,
                          osvs_, nosvs_, F_osv_diag_, pno_canonical_);
  }

  ~PNOEEPred() = default;

  virtual void operator()(
      const EigenVector<typename Tile::numeric_type> &e,
      std::vector<::mpqc::cc::TPack<Array>> &guess) const override {
    TA_ASSERT(e.size() == guess.size());
    //    TA_ASSERT(e.size() == n_roots_);

    // precondition
    const auto n = guess.size();
    for (std::size_t j = 0; j < n; ++j) {
      compute(e[j], guess[j]);
    }
  }

  /// override the default norm function
  virtual typename Tile::numeric_type norm(
      const ::mpqc::cc::TPack<Array> &t1t2) const override {
    auto r1_reblock = detail::reblock_t1(t1t2[0], reblock_i_, reblock_a_);
    auto r2_reblock = detail::reblock_t2(t1t2[1], reblock_i_, reblock_a_);

    detail::R1SquaredNormReductionOp<Array> op1(osvs_);
    detail::R2SquaredNormReductionOp<Array> op2(pnos_);

    return sqrt(r1_reblock("a,i").reduce(op1).get() +
                r2_reblock("a,b,i,j").reduce(op2).get()) /
           (size(r1_reblock) + size(r2_reblock));
  }

  /// return pno truncation threshold
  double tpno() const { return tpno_; }

  /// return osv truncation threshold
  double tosv() const { return tosv_; }

  void compute(const typename Tile::numeric_type &e,
               ::mpqc::cc::TPack<Array> &guess) const {
    // reblock
    guess[0] = detail::reblock_t1(guess[0], reblock_i_, reblock_a_);
    guess[1] = detail::reblock_t2(guess[1], reblock_i_, reblock_a_);

    // pno update
    guess[0] =
        detail::pno_jacobi_update_t1(guess[0], eps_o_, F_osv_diag_, osvs_, -e);
    guess[1] =
        detail::pno_jacobi_update_t2(guess[1], eps_o_, F_pno_diag_, pnos_, -e);

    // unblock
    guess[0] = detail::unblock_t1(guess[0], reblock_i_, reblock_a_);
    guess[1] = detail::unblock_t2(guess[1], reblock_i_, reblock_a_);

    // change the sign
    guess[0]("a,i") = -guess[0]("a,i");
    guess[1]("a,b,i,j") = -guess[1]("a,b,i,j");
  }

  void init_reblock(const Array &T2) {
    // initialize reblock array
    std::size_t n_unocc = T2.trange().dim(0).extent();
    std::size_t n_occ = T2.trange().dim(2).extent();
    auto &world = T2.world();

    const TA::TiledRange1 occ_row = T2.trange().dim(3);

    std::vector<std::size_t> occ_blocks;
    for (std::size_t i = 0; i <= n_occ; ++i) {
      occ_blocks.push_back(i);
    }
    const TA::TiledRange1 occ_col =
        TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());

    reblock_i_ = mpqc::array_ops::create_diagonal_array_from_eigen<
        Tile, TA::detail::policy_t<Array>>(world, occ_row, occ_col, 1.0);

    const TA::TiledRange1 uocc_row = T2.trange().dim(0);

    std::vector<std::size_t> uocc_blocks{0, n_unocc};
    const TA::TiledRange1 uocc_col =
        TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());

    reblock_a_ = mpqc::array_ops::create_diagonal_array_from_eigen<
        Tile, TA::detail::policy_t<Array>>(world, uocc_row, uocc_col, 1.0);
  }

 protected:
  double tpno_;
  double tosv_;
  bool pno_canonical_;
  std::size_t n_roots_;
  // diagonal of F_ij matrix
  Vector eps_o_;

  Array reblock_i_;
  Array reblock_a_;

  /// pnos for excited states
  std::vector<Matrix> pnos_;
  /// # of pnos for each pair
  std::vector<int> npnos_;
  /// diagonal of F_ab in pno basis for excited states
  std::vector<Vector> F_pno_diag_;
  /// osvs for excited states
  std::vector<Matrix> osvs_;
  /// # of osvs for each orbital
  std::vector<int> nosvs_;
  /// diagonal of F_ab in osv basis for excited states
  std::vector<Vector> F_osv_diag_;
};

/// state average PNO preconditioner for EOM-CCSD, the only difference is the
/// constructor

template <typename Array>
class StateAveragePNOEEPred : public PNOEEPred<Array> {
 public:
  using typename PNOEEPred<Array>::Matrix;
  using typename PNOEEPred<Array>::Vector;
  StateAveragePNOEEPred() = default;

  StateAveragePNOEEPred(const std::vector<Array> &T2, std::size_t n_roots,
                        const Vector &eps_o, const Matrix &F_uocc, double tpno,
                        double tosv, bool pno_canonical)
      : PNOEEPred<Array>() {
    this->tpno_ = tpno;
    this->tosv_ = tosv;
    this->pno_canonical_ = pno_canonical;
    this->n_roots_ = n_roots;
    this->eps_o_ = eps_o;

    this->init_reblock(T2[0]);

    // compute average of D
    Array D;
    {
      Array Ds;
      for (std::size_t i = 0; i < n_roots; i++) {
        auto T_reblock =
            detail::reblock_t2(T2[i], this->reblock_i_, this->reblock_a_);
        Ds = detail::construct_density(T_reblock);

        if (i == 0) {
          D("a,b,i,j") = Ds("a,b,i,j");
        } else {
          D("a,b,i,j") += Ds("a,b,i,j");
        }
      }

      D("a,b,i,j") = typename Array::element_type(1.0 / n_roots) * D("a,b,i,j");
    }

    detail::construct_pno(D, F_uocc,
                          this->tpno_, this->tosv_,
                          this->pnos_, this->npnos_, this->F_pno_diag_,
                          this->osvs_, this->nosvs_, this->F_osv_diag_,
                          this->pno_canonical_);
  }
};

/// state merged PNO preconditioner for EOM-CCSD, the only difference is the
/// constructor

template <typename Array>
class StateMergedPNOEEPred : public PNOEEPred<Array> {
public:
  using typename PNOEEPred<Array>::Matrix;
  using typename PNOEEPred<Array>::Vector;
  StateMergedPNOEEPred() = default;

  StateMergedPNOEEPred(const std::vector<Array> &T2, std::size_t n_roots,
                        const Vector &eps_o, const Matrix &F_uocc, double tpno,
                        double tosv, bool pno_canonical)
      : PNOEEPred<Array>() {
    this->tpno_ = tpno;
    this->tosv_ = tosv;
    this->pno_canonical_ = pno_canonical;
    this->n_roots_ = n_roots;
    this->eps_o_ = eps_o;

    this->init_reblock(T2[0]);

    // compute average of D
    Array D;
    {
      Array Ds;
      for (std::size_t i = 0; i < n_roots; i++) {
        auto T_reblock =
            detail::reblock_t2(T2[i], this->reblock_i_, this->reblock_a_);
        Ds = detail::construct_density(T_reblock);

        if (i == 0) {
          D("a,b,i,j") = Ds("a,b,i,j");
        } else {
          D("a,b,i,j") += Ds("a,b,i,j");
        }
      }
    }

    detail::construct_pno(D, F_uocc,
                          this->tpno_, this->tosv_,
                          this->pnos_, this->npnos_, this->F_pno_diag_,
                          this->osvs_, this->nosvs_, this->F_osv_diag_,
                          this->pno_canonical_);
  }
};

/// State Specific PNO preconditioner for EOM-CCSD

template <typename Array>
class StateSpecificPNOEEPred
    : public DavidsonDiagPred<::mpqc::cc::TPack<Array>> {
 public:
  using Tile = typename Array::value_type;
  using Matrix = RowMatrix<typename Tile::numeric_type>;
  using Vector = EigenVector<typename Tile::numeric_type>;

  StateSpecificPNOEEPred(const std::vector<Array> &T2, const Vector &eps_o,
                         const Matrix &F_uocc, double tpno, double tosv,
                         bool pno_canonical)
      : eps_o_(eps_o), tpno_(tpno), tosv_(tosv), pno_canonical_(pno_canonical) {
    auto &world = T2[0].world();
    n_roots_ = T2.size();

    // resize for n roots
    pnos_.resize(n_roots_);
    npnos_.resize(n_roots_);
    F_pno_diag_.resize(n_roots_);
    osvs_.resize(n_roots_);
    nosvs_.resize(n_roots_);
    F_osv_diag_.resize(n_roots_);

    // initialize reblock array
    {
      std::size_t n_unocc = T2[0].trange().dim(0).extent();
      std::size_t n_occ = T2[0].trange().dim(2).extent();

      const TA::TiledRange1 occ_row = T2[0].trange().dim(3);

      std::vector<std::size_t> occ_blocks;
      for (std::size_t i = 0; i <= n_occ; ++i) {
        occ_blocks.push_back(i);
      }
      const TA::TiledRange1 occ_col =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());

      reblock_i_ = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<Array>>(world, occ_row, occ_col, 1.0);

      const TA::TiledRange1 uocc_row = T2[0].trange().dim(0);

      std::vector<std::size_t> uocc_blocks{0, n_unocc};
      const TA::TiledRange1 uocc_col =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());

      reblock_a_ = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<Array>>(world, uocc_row, uocc_col, 1.0);
    }

    for (std::size_t i = 0; i < n_roots_; i++) {
      auto T_reblock = detail::reblock_t2(T2[i], reblock_i_, reblock_a_);
      auto D = detail::construct_density(T_reblock);
      detail::construct_pno(D, F_uocc,
                            tpno_, tosv_,
                            pnos_[i], npnos_[i], F_pno_diag_[i],
                            osvs_[i], nosvs_[i], F_osv_diag_[i],
                            pno_canonical_);
    }
  }

  ~StateSpecificPNOEEPred() = default;

  virtual void operator()(const EigenVector<typename Tile::numeric_type> &e,
                          std::vector<::mpqc::cc::TPack<Array>> &guess) const {
    TA_ASSERT(e.size() == guess.size());
    TA_ASSERT(e.size() == n_roots_);

    // precondition for each root
    for (std::size_t j = 0; j < n_roots_; ++j) {
      compute(j, e[j], guess[j]);
    }
  }

  void compute(std::size_t i, const typename Tile::numeric_type &e,
               ::mpqc::cc::TPack<Array> &guess) const {
    // reblock
    guess[0] = detail::reblock_t1(guess[0], reblock_i_, reblock_a_);
    guess[1] = detail::reblock_t2(guess[1], reblock_i_, reblock_a_);

    // pno update
    guess[0] = detail::pno_jacobi_update_t1(guess[0], eps_o_, F_osv_diag_[i],
                                            osvs_[i], -e);
    guess[1] = detail::pno_jacobi_update_t2(guess[1], eps_o_, F_pno_diag_[i],
                                            pnos_[i], -e);

    // unblock
    guess[0] = detail::unblock_t1(guess[0], reblock_i_, reblock_a_);
    guess[1] = detail::unblock_t2(guess[1], reblock_i_, reblock_a_);

    // change the sign
    guess[0]("a,i") = -guess[0]("a,i");
    guess[1]("a,b,i,j") = -guess[1]("a,b,i,j");
  }

  /// return pnos for root i
  const std::vector<Matrix> &pnos(int i) const { return pnos_[i]; }

  /// return osvs for root i
  const std::vector<Matrix> &osvs(int i) const { return osvs_[i]; }

  /// return pno truncation threshold
  double tpno() const { return tpno_; }

  /// return osv truncation threshold
  double tosv() const { return tosv_; }

 private:
  double tpno_;
  double tosv_;
  bool pno_canonical_;
  std::size_t n_roots_;

  Array reblock_i_;
  Array reblock_a_;

  // diagonal of F_ij matrix
  Vector eps_o_;

  /// pnos for all roots
  std::vector<std::vector<Matrix>> pnos_;
  /// # of pnos for each pair
  std::vector<std::vector<int>> npnos_;
  /// diagonal of F_ab in pno basis for all roots
  std::vector<std::vector<Vector>> F_pno_diag_;
  /// osvs for all roots
  std::vector<std::vector<Matrix>> osvs_;
  /// # of osvs for each orbital
  std::vector<std::vector<int>> nosvs_;
  /// diagonal of F_ab in osv basis for all roots
  std::vector<std::vector<Vector>> F_osv_diag_;
};

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_PRECONDITIONER_H_
