

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCHWARZ_SCREEN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCHWARZ_SCREEN_H_

#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/math/groups/petite_list.h"

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
  return map.diagonal().norm();
}

// Helper function to allow inf norm screening
inline double inf_norm(Eigen::Map<const RowMatrixXd> const &map) {
  return map.diagonal().lpNorm<Eigen::Infinity>();
}
}  // namespace detail

/*! \brief Class which holds shell set information for screening.  */
class Qmatrix {
 private:
  using f2s_map = std::unordered_map<int64_t, int64_t>;
  std::array<f2s_map, 2> f2s_maps_;  // Function to shell maps for each basis

  RowMatrixXd Q_;                    // Matrix to hold the integral estimates
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

  RowMatrixXd const &Q() const { return Q_; };

  /*! The basis constructor that takes a single argument is for the use case
   * where the Qmatrix is suppose to represent an Auxilary basisset.
   *
   * \param world is a reference to the madness world to compute the Q matrix
   * construction tasks
   *
   * \param eng is an integral engine pool used to compute the integral
   * estimates
   *
   * \param bs is a basis set which we will compute estimates via eng(shell,
   * unit, shell, unit);
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
  }

  /*! The double basis constructor that takes two basis sets is for the use case
   * where the Qmatrix is suppose to represent a basis with two center products
   *
   * \param world is a reference to the madness world to compute the Q matrix
   * construction tasks
   *
   * \param eng is an integral engine pool used to compute the integral
   * estimates
   *
   * \param bs0 is a basis set for the left center
   *
   * \param bs1 is a basis set for the right center
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
  std::shared_ptr<Qmatrix> Qab_;  // Screening mat for basis ab
  std::shared_ptr<Qmatrix> Qcd_;  // Screening mat for basis cd
  double thresh_;                 // Threshold used for screening
  double thresh2_;  // Threshold squared since estimates are not squared

  inline Qmatrix const &Qab() const { return *Qab_; }
  inline Qmatrix const &Qcd() const { return *Qcd_; }

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

 public:
  SchwarzScreen() = default;
  SchwarzScreen(SchwarzScreen const &) = default;
  SchwarzScreen(SchwarzScreen &&) = default;
  SchwarzScreen &operator=(SchwarzScreen const &) = default;
  SchwarzScreen &operator=(SchwarzScreen &&) = default;

  virtual ~SchwarzScreen() noexcept = default;

  /*! \brief Constructor which requires two Q matrices
   *
   * \param Qab is the Qmatrix needed for mu and nu in (mu nu | rho sigma) and
   * also the vector needed for X in (X | rho sigma), a current limitation of
   * Schwarz screening is that all auxillary basis Qs must be placed in Qab_.
   *
   * \param Qcd is the Qmatrix needed for rho and sigma in (mu nu | rho sigma)
   * and (X | rho sigma)
   *
   * \param thresh is the screening threshold used when evaluating whether or
   * not the skip function returns true for Qab(a,b) * Qcd(c,d) <= thresh_.
   * How Qab and Qcd are determined (infinity norm of a shell set or F norm) is
   * to be part of the Q matrix construction.
   */
  SchwarzScreen(std::shared_ptr<Qmatrix> Qab, std::shared_ptr<Qmatrix> Qcd,
                double thresh);

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

  /*! \brief returns an estimate of shape norms for the given basis vector.
   *
   * This will use the outer product of Qab * Qcd to determine the
   * shape. Will replicate the result on ever node. 
   */
  TA::Tensor<float> norm_estimate(
      madness::World &world, std::vector<gaussian::Basis> const &bs_array,
      const math::PetiteList &plist =
          math::SymmPetiteList<math::PetiteList::Symmetry::e>()) const override;
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
    double thresh, decltype(detail::inf_norm) norm_func = detail::inf_norm) {
  if (bs_array.size() == 3) {  // One index must be auxiliary assume first
    auto Q_ab =
        std::make_shared<Qmatrix>(Qmatrix(world, eng, bs_array[0], norm_func));
    auto Q_cd = std::make_shared<Qmatrix>(
        Qmatrix(world, eng, bs_array[1], bs_array[2], norm_func));
    return SchwarzScreen(std::move(Q_ab), std::move(Q_cd), thresh);
  } else if (bs_array.size() == 4) {  // Four center estimator
    auto Q_ab = std::make_shared<Qmatrix>(
        Qmatrix(world, eng, bs_array[0], bs_array[1], norm_func));
    auto Q_cd = std::make_shared<Qmatrix>(
        Qmatrix(world, eng, bs_array[2], bs_array[3], norm_func));
    return SchwarzScreen(std::move(Q_ab), std::move(Q_cd), thresh);
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
