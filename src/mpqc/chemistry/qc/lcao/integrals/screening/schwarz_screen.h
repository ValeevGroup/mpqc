

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCHWARZ_SCREEN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCHWARZ_SCREEN_H_

#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/math/groups/petite_list.h"
#include "mpqc/util/misc/exception.h"

#include <tiledarray.h>

#include <stdexcept>
#include <unordered_map>

#include <functional>

namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {
// Helper function to allow F norm screening
inline double l2Norm(Eigen::Map<const RowMatrixXd> const &map) {
  return map.norm();
}

// Helper function to allow inf norm screening
inline double inf_norm(Eigen::Map<const RowMatrixXd> const &map) {
  return map.lpNorm<Eigen::Infinity>();
}
}  // namespace detail

/// \brief Class which holds shell set information for screening.
class Qmatrix {
 private:
  using f2s_map = std::unordered_map<int64_t, int64_t>;
  std::array<f2s_map, 2> f2s_maps_;  // Function to shell maps for each basis

  // Matrix to hold the {shell,shell-pair} Schwarz bound factors
  // # of rows = # of shells in a basis
  // # of cols = {1,# of shells} if screening integrals over
  // {shells,shell-pairs}
  RowMatrixXd Q_;
  // Matrix to hold the cluster/cluster-pair Schwarz bound factors
  // # of rows = # of clusters in a basis
  // # of cols = {1,# of clusers} if screening integrals over
  // {clusters,cluster-pairs}
  RowMatrixXd Q_cluster_;
  Eigen::VectorXd max_elem_in_row_;  // Max element for each row of Q_
  double max_elem_in_Q_ = 0;  // Max val in Q_, should be set in constructor
  bool is_aux_;

  using norm_func_type = decltype(detail::inf_norm);

  // Function to return the shell index for a basis given the function
  // index(ind) The map passed in should correspond to the basis desired
  int64_t f2s(f2s_map const &, int64_t) const;

 public:
  Qmatrix() = default;
  Qmatrix(Qmatrix const &) = default;
  Qmatrix(Qmatrix &&) = default;
  Qmatrix &operator=(Qmatrix const &) = default;
  Qmatrix &operator=(Qmatrix &&) = default;

  RowMatrixXd const &Q() const { return Q_; }
  RowMatrixXd const &Qtile() const { return Q_cluster_; }

  /*! Constructs Qmatrix for the case of screening integrals over shells, e.g.
   * integrals of type \c (a|Op|b)
   *
   * \param world is a reference to the madness world to compute the Q matrix
   * construction tasks
   *
   * \param eng is an integral engine pool used to compute the integral
   * estimates
   *
   * \param bs is a Basis object
   *
   * \param norm_func is a function that computes the sqrt(norm(quartet)) of the
   * integrals with signature double (Eigen::Map<const RowMatrixXd> const &)
   */
  template <typename E>
  Qmatrix(madness::World &world, E const &eng, Basis const &bs,
          norm_func_type norm_op = detail::inf_norm)
      : is_aux_(true) {
    // Get bs shells
    auto shells = bs.flattened_shells();
    const auto nshells = shells.size();

    // Compute Q_
    Q_ = Eigen::VectorXd(nshells);

    auto Q_func = [=](Shell const *sh, double *Q_val) {
      // Compute the shell set (a unit | a unit)
      const auto &bufs = eng->local().compute(*sh, unit_shell, *sh, unit_shell);
      const auto nsh = sh->size();
      const auto bmap = Eigen::Map<const RowMatrixXd>(bufs[0], nsh, nsh);

      // Compute and write the sqrt(norm) of quartet (a unit | a unit)
      *Q_val = norm_op(bmap);
    };

    for (auto i = 0ul; i < nshells; ++i) {
      auto *Q_val_ptr = &Q_(i);
      world.taskq.add(Q_func, &shells[i], Q_val_ptr);
    }
    world.gop.fence();

    // Compute the function to shell map
    f2s_map map;
    auto func_ind = 0;
    for (auto i = 0ul; i < nshells; ++i) {
      map.emplace(func_ind, i);
      func_ind += shells[i].size();
    }
    f2s_maps_[0] = std::move(map);

    // Compute the max elem in Q
    max_elem_in_row_.resize(nshells);
    for (auto i = 0ul; i < nshells; ++i) {
      const auto max_elem = Q_(i);
      max_elem_in_Q_ = std::max(max_elem_in_Q_, max_elem);
      max_elem_in_row_[i] = max_elem;
    }

    // Build Q_cluster_
    auto const &cshells = bs.cluster_shells();
    const auto nclusters = cshells.size();
    Q_cluster_ = Eigen::VectorXd(nclusters);
    Q_cluster_.setZero();

    // Because we are using Qbra * Qket >= (bra|ket)^2 just sum these values
    // into
    // Q_cluster_
    auto first_shell_in_cluster = 0;
    for (auto c = 0; c < nclusters; ++c) {
      const auto nshells = cshells[c].size();
      const auto last = first_shell_in_cluster + nshells;
      for (auto s = first_shell_in_cluster; s < last; ++s) {
        Q_cluster_(c) += Q_(s);
      }
      first_shell_in_cluster += nshells;
    }
  }

  /*! Constructs Qmatrix for the case of screening integrals over shell-pairs,
   * e.g. integrals of type \c (ab|Op|cd)
   *
   * \param world is a reference to the madness world to compute the Q matrix
   * construction tasks
   *
   * \param eng is an integral engine pool used to compute the integral
   * estimates
   *
   * \param bs0 is a Basis object composed of shells that appear first in shell
   * pairs
   *
   * \param bs1 is a Basis object composed of shells that appear second in shell
   * pairs
   *
   * \param norm_func is a function that computes the sqrt(norm(quartet)) of the
   * integrals with signature double (Eigen::Map<const RowMatrixXd> const &)
   */
  template <typename E>
  Qmatrix(TA::World &world, E const &eng, Basis const &bs0, Basis const &bs1,
          norm_func_type norm_op = detail::inf_norm)
      : is_aux_(false) {
    // Get bs shells
    const auto shells0 = bs0.flattened_shells();
    const auto shells1 = bs1.flattened_shells();
    const auto nshells0 = shells0.size();
    const auto nshells1 = shells1.size();

    // Compute Q_
    Q_.resize(nshells0, nshells1);

    auto Q_func = [=](Shell const *sh0, Shell const *sh1, double *Q_val) {
      // Get shell sizes
      const auto nsh0 = sh0->size();
      const auto nsh1 = sh1->size();
      const auto n2 = nsh0 * nsh1;

      // Capture local engine
      auto &leng = eng->local();
      leng.set_precision(0.);
      auto const &result = leng.results();

      // Compute estimate
      leng.compute(*sh0, *sh1, *sh0, *sh1);
      const auto bmap = Eigen::Map<const RowMatrixXd>(result[0], n2, n2);
      *Q_val = norm_op(bmap);

      // Reset precision to something better than 0.
      leng.set_precision(std::numeric_limits<double>::epsilon());
    };

    for (auto i = 0ul; i < nshells0; ++i) {
      for (auto j = 0ul; j < nshells1; ++j) {
        auto *Q_val_ptr = &Q_(i, j);
        world.taskq.add(Q_func, &shells0[i], &shells1[j], Q_val_ptr);
      }
    }
    world.gop.fence();

    // Compute the function to shell map for bs0
    f2s_map map0;
    auto func_ind = 0;
    for (auto i = 0ul; i < nshells0; ++i) {
      map0.emplace(func_ind, i);
      func_ind += shells0[i].size();
    }
    f2s_maps_[0] = std::move(map0);

    // Compute the function to shell map for bs1
    f2s_map map1;
    func_ind = 0;
    for (auto i = 0ul; i < nshells1; ++i) {
      map1.emplace(func_ind, i);
      func_ind += shells1[i].size();
    }
    f2s_maps_[1] = std::move(map1);

    // Compute the max elem in Q
    max_elem_in_row_.resize(Q_.rows());
    for (auto i = 0ul; i < Q_.rows(); ++i) {
      const auto max_elem = Q_.row(i).cwiseAbs().maxCoeff();
      max_elem_in_Q_ = std::max(max_elem_in_Q_, max_elem);
      max_elem_in_row_[i] = max_elem;
    }

    // Make Q_cluster_
    auto const &cshells0 = bs0.cluster_shells();
    auto const &cshells1 = bs1.cluster_shells();

    auto const &nclusters0 = cshells0.size();
    auto const &nclusters1 = cshells1.size();

    Q_cluster_ = RowMatrixXd(nclusters0, nclusters1);
    Q_cluster_.setZero();

    auto first_shell_in_cluster0 = 0;
    for (auto c0 = 0; c0 < nclusters0; ++c0) {
      const auto nshells0 = cshells0[c0].size();
      const auto last0 = first_shell_in_cluster0 + nshells0;

      auto first_shell_in_cluster1 = 0;
      for (auto c1 = 0; c1 < nclusters1; ++c1) {
        const auto nshells1 = cshells1[c1].size();
        const auto last1 = first_shell_in_cluster1 + nshells1;

        auto val = 0.0;
        for (auto s0 = first_shell_in_cluster0; s0 < last0; ++s0) {
          for (auto s1 = first_shell_in_cluster1; s1 < last1; ++s1) {
            val += Q_(s0, s1);
          }
        }
        Q_cluster_(c0, c1) = val;

        first_shell_in_cluster1 += nshells1;
      }
      first_shell_in_cluster0 += nshells0;
    }
  }

  /// Will return the largest value in Q_
  double operator()() const;

  /// Will return the largest value in Q_.row(a).
  double operator()(int64_t a) const;

  /// Will return max in row but not do a f2s lookup
  double max_in_row(int64_t a) const;

  /// Will return the value Q_(a,b)
  double operator()(int64_t a, int64_t b) const;

  /// Return whether this Q is for an auxiliary basis, i.e. is a vector
  bool is_aux_Q() const;
};

/*! \brief Class for Schwarz based screening.
 *
 * We will assume that these screeners are replicated for the integrals they
 * need and so will not bother with serialization of them.
 *
 * SchwarzScreen has support for up to 4 different bases, currently there are
 * no optimizations made for all four bases being equal to each other.
 */
class SchwarzScreen : public Screener {
 private:
  std::shared_ptr<Qmatrix> Qbra_;  // Screening matrix for the bra
  std::shared_ptr<Qmatrix> Qket_;  // Screening matrix for the ket
  double thresh_;                  // Threshold used for screening
  double thresh2_;  // Threshold squared since estimates are not squared

  inline boost::optional<double> estimate(int64_t a) const;

  inline boost::optional<double> estimate(int64_t a, int64_t b) const;

  inline boost::optional<double> estimate(int64_t a, int64_t b,
                                          int64_t c) const;

  inline boost::optional<double> estimate(int64_t a, int64_t b, int64_t c,
                                          int64_t d) const;

  /* skip_ takes a list of indices and returns true when the integral estimate
   * is below the threshold and false when it is greater than or equal to the
   * threshold. See function estimate for the implementation.
   *
   * \param idx is a series of indices from which to generate an integral
   * estimate
   *
   * This is a one stop catch all for the overloads of the skip functions
   */
  template <typename... IDX>
  bool skip_(IDX... idx) const {
    auto est =
        estimate(std::forward<IDX>(idx)...);  // optional estimation for set idx
    if (est) {  // Check that we were able to estimate the integral set
      return est < thresh2_;  // If est below threshold then skip this set
    } else {         // We were unable to estimate this set for some reason
      return false;  // thus we should compute the integral
    }
  }

  /* skip_ takes norm of a density block and a list of indices.
   * This returns true when the product of density norm and integral estimate
   * is below the threshold and false when it is greater than or equal to the
   * threshold. See function estimate for the implementation.
   *
   * \param D is infinity norm of the corresponding density block
   * \param idx is a series of indices from which to generate an integral
   * estimate
   *
   * This is a one stop catch all for the overloads of the skip functions
   */
  template <typename... IDX>
  bool skip_(double D, IDX... idx) const {
    auto est =
        estimate(std::forward<IDX>(idx)...);  // optional estimation for set idx
    if (est && D != 0.0) {  // Check that we were able to estimate the integral
                            // set AND density norm is large than zero
      return (D * D * est.get()) <
             thresh2_;  // If D^2 * est below threshold then skip this set
    } else {            // We were unable to estimate this set for some reason
      return false;     // thus we should compute the integral
    }
  }

 public:
  SchwarzScreen() = default;
  SchwarzScreen(SchwarzScreen const &) = default;
  SchwarzScreen(SchwarzScreen &&) = default;
  SchwarzScreen &operator=(SchwarzScreen const &) = default;
  SchwarzScreen &operator=(SchwarzScreen &&) = default;

  virtual ~SchwarzScreen() noexcept = default;

  /*! \brief Constructor which requires two Q matrices
   *
   * \param Qbra is the Qmatrix that describes the bra in screened integrals,
   * e.g.
   * it represents (cd| in (ab|cd), or (X| in (X|cd)
   *
   * \param Qket is the Qmatrix that describes the ket in screened integrals,
   * e.g.
   * it represents |cd) in (ab|cd), or |X) in (ab|X)
   *
   * \param thresh is the screening threshold used when evaluating whether or
   * not the skip function returns true for Qbra(a,b) * Qket(c,d) <= thresh_.
   * The norm types of Qbra and Qket are defined by their constructors.
   */
  SchwarzScreen(std::shared_ptr<Qmatrix> Qbra, std::shared_ptr<Qmatrix> Qket,
                double thresh);

  /// Returns the screening matrix for the Bra indices
  inline Qmatrix const &Qbra() const { return *Qbra_; }

  /// Returns the screening matrix for the Ket indices
  inline Qmatrix const &Qket() const { return *Qket_; }

  /// Reports the threshold being used for skipping integrals
  double skip_threshold() const;

  // Overrides of the skips following.  See function skip_ for implementation
  bool skip(int64_t) override;
  bool skip(int64_t) const override;

  bool skip(int64_t, int64_t) override;
  bool skip(int64_t, int64_t) const override;

  bool skip(int64_t, int64_t, int64_t) override;
  bool skip(int64_t, int64_t, int64_t) const override;

  bool skip(int64_t, int64_t, int64_t, int64_t) override;
  bool skip(int64_t, int64_t, int64_t, int64_t) const override;

  bool skip(int64_t, int64_t, int64_t, int64_t, double) override;
  bool skip(int64_t, int64_t, int64_t, int64_t, double) const override;

  /*! \brief returns an estimate of shape norms for the given basis vector, in
   * presence of symmetry described by a math::PetiteList object.
   *
   * This function will only compute the estimate for tiles which are
   * considered local by the pmap, the user will be responsible for using the
   * world based constructor for the TA::Shape.
   */
  TA::Tensor<float> norm_estimate(madness::World &world,
                                  std::vector<gaussian::Basis> const &bs_array,
                                  TA::Pmap const &pmap,
                                  bool replicate = false) const override;
};

/*! \brief Creates a Schwarz Screener
 *
 * \param world is a reference to a madness world
 *
 * \param eng is a reference to a ShrPool<E>
 *
 * \param bs_array is a reference to a vector of basis sets, if the length is 3
 * then DF integrals are assumed and the first basis is assumed to be the
 * auxiliary basis if the length is 4 then it is assumed that four center
 * screening is desired.  There is no requirement that any basis sets be the
 * same.
 *
 * \param thresh is the Schwarz Screening threshold
 *
 * \param norm_func is the function pointer that returns a norm, it should have
 * a signature double (Eigen::Map<const RowMatrixXd> const &)
 *
 */
template <typename E>
SchwarzScreen create_schwarz_screener(
    TA::World &world, ShrPool<E> const &eng, std::vector<Basis> const &bs_array,
    double thresh, decltype(detail::inf_norm) norm_func = detail::l2Norm) {
  if (bs_array.size() == 3) {  // One index must be auxiliary assume first
    auto Q_bra =
        std::make_shared<Qmatrix>(Qmatrix(world, eng, bs_array[0], norm_func));
    auto Q_ket = std::make_shared<Qmatrix>(
        Qmatrix(world, eng, bs_array[1], bs_array[2], norm_func));
    return SchwarzScreen(std::move(Q_bra), std::move(Q_ket), thresh);
  } else if (bs_array.size() == 4) {  // Four center estimator
    auto Q_bra = std::make_shared<Qmatrix>(
        Qmatrix(world, eng, bs_array[0], bs_array[1], norm_func));
    auto Q_ket = std::make_shared<Qmatrix>(
        Qmatrix(world, eng, bs_array[2], bs_array[3], norm_func));
    return SchwarzScreen(std::move(Q_bra), std::move(Q_ket), thresh);
  } else {
    throw InputError(
        "The bs_array passed into create_schwarz_screener must be of length 3 "
        "or 4",
        __FILE__, __LINE__);
    return SchwarzScreen();
  }
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCHWARZ_SCREEN_H_
