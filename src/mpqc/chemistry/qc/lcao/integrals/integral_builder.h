
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_BUILDER_H_

#include <array>
#include <functional>
#include <memory>

#include <tiledarray.h>

#include "mpqc/util/misc/pool.h"

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integral_kernels.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

/*! \brief Builds integrals from an  array of bases and an integral engine pool.
 *
 * \param Op is a function or functor that takes a TA::Tensor && and returns a
 *tile.
 * The simplest type of Op is simply the following:
 * ```
 * auto op = [](TA::Tensor<double> && t){ return std::move(t) };
 * ```
 */
template <typename Tile, typename Engine = libint2::Engine>
class IntegralBuilder
    : public std::enable_shared_from_this<IntegralBuilder<Tile, Engine>> {
 public:
  using Op = std::function<Tile(TA::TensorD &&)>;

 private:
  std::shared_ptr<BasisVector> bases_;
  ShrPool<Engine> engines_;
  std::shared_ptr<Screener> screen_;
  Op op_;

 public:
  /*! \brief Constructor which copies all shared_ptr members
   *
   * \param world Should be the same world as the one for the array which will
   * hold the tiles.
   * \param shr_epool is a shared pointer to an IntegralTSPool
   * \param shr_bases is a shared pointer to an array of Basis
   * \param screen is a shared pointer to a Screener type
   * \param op should be a thread safe function or functor that takes a
   *  rvalue of a TA::TensorD and returns a valid TA::Array tile.
   */
  IntegralBuilder(ShrPool<Engine> shr_epool, std::shared_ptr<BasisVector> shr_bases,
                  std::shared_ptr<Screener> screen, Op op)
      : bases_(std::move(shr_bases)),
        engines_(std::move(shr_epool)),
        screen_(std::move(screen)),
        op_(std::move(op)) {
    std::size_t N = bases_->size();
    TA_ASSERT((N == 2) || (N == 3) || (N == 4));
  }

  virtual ~IntegralBuilder() = default;

  Tile operator()(std::vector<std::size_t> const &idx, TA::Range range) {
    return op_(integrals(idx, std::move(range)));
  }

  TA::TensorD integrals(std::vector<std::size_t> const &idx, TA::Range range) {
    auto size = bases_->size();

    if (size == 2) {
      // Get integral shells
      detail::VecArray<2> shellvec_ptrs;
      for (auto i = 0ul; i < size; ++i) {
        auto const &basis_i = bases_->operator[](i);
        shellvec_ptrs[i] = &basis_i.cluster_shells()[idx[i]];
      }

      // Compute integrals over the selected shells.
      return detail::integral_kernel(engines_->local(), std::move(range),
                                     shellvec_ptrs, *screen_);
    } else if (size == 3) {
      // Get integral shells
      detail::VecArray<3> shellvec_ptrs;
      for (auto i = 0ul; i < size; ++i) {
        auto const &basis_i = bases_->operator[](i);
        shellvec_ptrs[i] = &basis_i.cluster_shells()[idx[i]];
      }

      // Compute integrals over the selected shells.
      return detail::integral_kernel(engines_->local(), std::move(range),
                                     shellvec_ptrs, *screen_);
    } else if (size == 4) {
      // Get integral shells
      detail::VecArray<4> shellvec_ptrs;
      for (auto i = 0ul; i < size; ++i) {
        auto const &basis_i = bases_->operator[](i);
        shellvec_ptrs[i] = &basis_i.cluster_shells()[idx[i]];
      }

      // Compute integrals over the selected shells.
      return detail::integral_kernel(engines_->local(), std::move(range),
                                     shellvec_ptrs, *screen_);
    } else {
      throw std::runtime_error(
          "Invalid Size of Basis Sets!! Must be 2 or 3 or 4!! \n");
    }
  }

  Tile op(TA::TensorD &&tensor) { return op_(std::move(tensor)); }
};

template <typename Tile, typename Engine = libint2::Engine>
class DirectIntegralBuilder : public IntegralBuilder<Tile, Engine> {
 public:
  using Op = typename IntegralBuilder<Tile,Engine>::Op;

  DirectIntegralBuilder(madness::World &world, ShrPool<Engine> shr_epool,
                        std::shared_ptr<BasisVector> shr_bases,
                        std::shared_ptr<Screener> screen, Op op)
      : IntegralBuilder<Tile, Engine>(shr_epool, shr_bases, screen, op),
        id_(world.register_ptr(this)) {}

  madness::uniqueidT id() const { return id_; }

  ~DirectIntegralBuilder() {
    if (madness::initialized()) {
      madness::World *world = madness::World::world_from_id(id_.get_world_id());
      world->unregister_ptr(this);
    }
  }

  using IntegralBuilder<Tile,Engine>::operator();
  using IntegralBuilder<Tile,Engine>::integrals;
  using IntegralBuilder<Tile,Engine>::op;

 private:
  madness::uniqueidT id_;
};

/*!
 * \brief Function to make detection of template parameters easier, see
 * IntegralBuilder for details.
 */
template <typename Tile, typename Engine>
std::shared_ptr<IntegralBuilder<Tile, Engine>> make_integral_builder(
    ShrPool<Engine> shr_epool, std::shared_ptr<BasisVector> shr_bases,
    std::shared_ptr<Screener> shr_screen,
    std::function<Tile(TA::TensorD &&)> op) {
  return std::make_shared<IntegralBuilder<Tile, Engine>>(
      std::move(shr_epool), std::move(shr_bases), std::move(shr_screen),
      std::move(op));
}

/*!
 * \brief Function to make detection of template parameters easier, see
 * DirectIntegralBuilder for details.
 */
template <typename Tile, typename Engine>
std::shared_ptr<DirectIntegralBuilder<Tile, Engine>> make_direct_integral_builder(
    madness::World &world, ShrPool<Engine> shr_epool,
    std::shared_ptr<BasisVector> shr_bases, std::shared_ptr<Screener> shr_screen,
    std::function<Tile(TA::TensorD &&)> op) {
  return std::make_shared<DirectIntegralBuilder<Tile, Engine>>(
      world, std::move(shr_epool), std::move(shr_bases), std::move(shr_screen),
      std::move(op));
}

/*!
 * \brief Function to make detection of template parameters easier for complex
 * tensor types, see DirectIntegralBuilder for details.
 */
template <typename Tile, typename Engine>
std::shared_ptr<DirectIntegralBuilder<Tile, Engine>> make_direct_integral_builder(
    madness::World &world, ShrPool<Engine> shr_epool,
    std::shared_ptr<BasisVector> shr_bases, std::shared_ptr<Screener> shr_screen,
    std::function<Tile(TA::TensorZ &&)> op) {
  return std::make_shared<DirectIntegralBuilder<Tile, Engine>>(
      world, std::move(shr_epool), std::move(shr_bases), std::move(shr_screen),
      std::move(op));
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_BUILDER_H_
