
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_BUILDER_H_

#include <array>
#include <functional>
#include <memory>

#include <tiledarray.h>

#include "mpqc/util/misc/pool.h"

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/expression/formula.h"
#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integral_kernels.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/math/groups/petite_list.h"
#include "mpqc/util/core/exception.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

/*! \brief Builds integrals from an  array of bases and an integral engine pool.
 *
 * \param Op is a function or functor that takes a TA::Tensor && and returns a
 * tile. The simplest type of Op is simply the following:
 * ```
 * auto op = [](TA::Tensor<double> && t){ return std::move(t) };
 * ```
 */
template <typename Tile, typename Engine = libint2::Engine>
class IntegralBuilder {
 public:
  using Op = std::function<Tile(TA::TensorD &&)>;

 private:
  std::shared_ptr<BasisVector> bases_;
  ShrPool<Engine> engines_;
  std::shared_ptr<Screener> screen_;
  Op op_;
  std::shared_ptr<const math::PetiteList> plist_;

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
   * \param plist the PetiteList object describing the symmetry properties of
   * the set of AO integrals
   */
  IntegralBuilder(ShrPool<Engine> shr_epool,
                  std::shared_ptr<BasisVector> shr_bases,
                  std::shared_ptr<Screener> screen, Op op,
                  std::shared_ptr<const math::PetiteList> plist =
                      math::PetiteList::make_trivial())
      : bases_(std::move(shr_bases)),
        engines_(std::move(shr_epool)),
        screen_(std::move(screen)),
        op_(std::move(op)),
        plist_(std::move(plist)) {
    std::size_t N = bases_->size();
    TA_ASSERT((N == 2) || (N == 3) || (N == 4));
  }

  virtual ~IntegralBuilder() {}

  Tile operator()(std::vector<std::size_t> const &idx, TA::Range range) {
    return op_(integrals(idx, std::move(range)));
  }

  TA::TensorD integrals(std::vector<std::size_t> const &idx, TA::Range range) {
    auto size = bases_->size();

    switch (size) {
      case 2:
        // Get integral shells
        detail::VecArray<2> shellvec_ptrs2;
        for (auto i = 0ul; i < size; ++i) {
          auto const &basis_i = bases_->operator[](i);
          shellvec_ptrs2[i] = &basis_i.cluster_shells()[idx[i]];
        }
        // Compute integrals over the selected shells.
        return detail::integral_kernel(engines_->local(), std::move(range),
                                       shellvec_ptrs2, *screen_, *plist_);
      case 3:
        // Get integral shells
        detail::VecArray<3> shellvec_ptrs3;
        for (auto i = 0ul; i < size; ++i) {
          auto const &basis_i = bases_->operator[](i);
          shellvec_ptrs3[i] = &basis_i.cluster_shells()[idx[i]];
        }
        // Compute integrals over the selected shells.
        return detail::integral_kernel(engines_->local(), std::move(range),
                                       shellvec_ptrs3, *screen_, *plist_);
      case 4:
        // Get integral shells
        detail::VecArray<4> shellvec_ptrs4;
        for (auto i = 0ul; i < size; ++i) {
          auto const &basis_i = bases_->operator[](i);
          shellvec_ptrs4[i] = &basis_i.cluster_shells()[idx[i]];
        }

        // Compute integrals over the selected shells.
        return detail::integral_kernel(engines_->local(), std::move(range),
                                       shellvec_ptrs4, *screen_, *plist_);
      default:
        throw std::runtime_error(
            "Invalid Size of Basis Sets!! Must be 2 or 3 or 4!! \n");
    }
  }

  /*!
   * @brief This computes an \c std::array of integrals for operators that have
   * \c nopers components, e.g. multipole moments, geometrical derivatives,
   * etc. Note that only one tile is computed for each component, and all
   * \c nopers tiles share the same index and range.
   * @tparam libint2_oper libint2 operator type
   * @param idx index of the target tile
   * @param range range of the target tile
   * @return
   */
  template <libint2::Operator libint2_oper>
  std::array<TA::TensorD, libint2::operator_traits<libint2_oper>::nopers>
  integrals(std::vector<std::size_t> const &idx, TA::Range range) {
    auto size = bases_->size();

    switch (size) {
      case 2:
        // get integral shells
        detail::VecArray<2> shellvec_ptrs2;
        for (auto i = 0ul; i < size; ++i) {
          auto const &basis_i = bases_->operator[](i);
          shellvec_ptrs2[i] = &basis_i.cluster_shells()[idx[i]];
        }
        // Compute integrals over the selected shells.
        return detail::integral_kernel<libint2_oper>(
            engines_->local(), std::move(range), shellvec_ptrs2, *screen_,
            *plist_);

      default:
        throw FeatureNotImplemented("Only support basis set size = 2!",
                                    __FILE__, __LINE__);
    }
  };

  Tile op(TA::TensorD &&tensor) { return op_(std::move(tensor)); }
};

template <typename Tile, typename Engine = libint2::Engine>
class DirectIntegralBuilder
    : public IntegralBuilder<Tile, Engine>,
      public std::enable_shared_from_this<DirectIntegralBuilder<Tile, Engine>> {
 public:
  using Op = typename IntegralBuilder<Tile, Engine>::Op;

  DirectIntegralBuilder(madness::World &world, ShrPool<Engine> shr_epool,
                        std::shared_ptr<BasisVector> shr_bases,
                        std::shared_ptr<Screener> screen, Op op,
                        std::shared_ptr<const math::PetiteList> plist)
      : IntegralBuilder<Tile, Engine>(shr_epool, shr_bases, screen, op, plist),
        id_(world.register_ptr(this)) {}

  DirectIntegralBuilder(const DirectIntegralBuilder &builder) = delete;
  DirectIntegralBuilder &operator=(const DirectIntegralBuilder &builder) =
      delete;

  madness::uniqueidT id() const { return id_; }

  ~DirectIntegralBuilder() {
    if (madness::initialized()) {
      madness::World *world = madness::World::world_from_id(id_.get_world_id());
      world->unregister_ptr(this);
    }
  }

  using IntegralBuilder<Tile, Engine>::operator();
  using IntegralBuilder<Tile, Engine>::integrals;
  using IntegralBuilder<Tile, Engine>::op;

 private:
  madness::uniqueidT id_;
};

template <typename Tile>
struct DFTaskGemm {
  typedef Tile result_type;
  typedef Tile first_argument_type;
  typedef Tile second_argument_type;

  TA::Range range;
  TA::math::GemmHelper gemm_helper;

  DFTaskGemm(const TA::Range &range, const TA::math::GemmHelper &helper)
      : range(range), gemm_helper(helper) {}

  result_type operator()() const { return Tile(range, 0.0); }

  const result_type &operator()(const result_type &result) const {
    return result;
  }

  void add(result_type &result, const result_type &arg) const {
    TA::math::inplace_vector_op_serial(
        [](TA::detail::numeric_t<Tile> &l,
           const TA::detail::numeric_t<Tile> r) { l += r; },
        result.range().volume(), result.data(), arg.data());
  }

  void operator()(result_type &result, const result_type &arg) const {
    add(result, arg);
  }

  void operator()(result_type &result, const first_argument_type &first,
                  const second_argument_type &second) const {
    Tile tmp = first.gemm(second, 1.0, gemm_helper);
    add(result, tmp);
  }
};

template <typename Tile, typename Policy>
class DirectDFIntegralBuilder : public std::enable_shared_from_this<
                                    DirectDFIntegralBuilder<Tile, Policy>> {
 public:
  // constructor
  DirectDFIntegralBuilder() = default;
  DirectDFIntegralBuilder(
      const TA::DistArray<Tile, Policy> &left,
      const TA::DistArray<Tile, Policy> &right,
      Formula::Notation notation = Formula::Notation::Chemical)
      : bra_(left),
        ket_(right),
        id_(left.world().register_ptr(this)),
        notation_(notation) {
    df_lobound_ = bra_.trange().data().front().tiles_range().first;
    df_upbound_ = bra_.trange().data().front().tiles_range().second;

    TA_ASSERT(df_lobound_ == ket_.trange().data().front().tiles_range().first);
    TA_ASSERT(df_upbound_ == ket_.trange().data().front().tiles_range().second);
  }

  DirectDFIntegralBuilder(DirectDFIntegralBuilder &&) = delete;
  DirectDFIntegralBuilder(const DirectDFIntegralBuilder &) = delete;
  DirectDFIntegralBuilder &operator=(const DirectDFIntegralBuilder &) = delete;

  ~DirectDFIntegralBuilder() {
    if (madness::initialized()) {
      madness::World *world = madness::World::world_from_id(id_.get_world_id());
      world->unregister_ptr(this);
    }
  }

  madness::uniqueidT id() const { return id_; }

  // compute Tile for particular block
  madness::Future<Tile> operator()(const std::vector<std::size_t> &idx,
                                   const TA::Range &range) const {
    TA_ASSERT(idx.size() == 4);

    auto &world = bra_.world();
    // create tile
    //    madness::Future<Tile> result(Tile(range, 0.0));
    //        Tile result(range, 0.0);

    std::vector<std::size_t> bra_idx(3);
    std::vector<std::size_t> ket_idx(3);
    bra_idx[1] = idx[0];
    ket_idx[2] = idx[3];
    if (notation_ == Formula::Notation::Physical) {
      bra_idx[2] = idx[2];
      ket_idx[1] = idx[1];
    } else {
      bra_idx[2] = idx[1];
      ket_idx[1] = idx[2];
    }

    TA::math::GemmHelper gemm_helper(madness::cblas::Trans,
                                     madness::cblas::NoTrans, 4, 3, 3);

    TA::Range this_range;
    if (notation_ == Formula::Notation::Physical) {
      this_range = TA::Range(TA::Permutation({0, 2, 1, 3}), range);
    } else {
      this_range = range;
    }

    DFTaskGemm<Tile> task_gemm(this_range, gemm_helper);

    TA::detail::ReducePairTask<decltype(task_gemm)> reduce_pair_task(world,
                                                                     task_gemm);
    // loop over density fitting space
    for (std::size_t i = df_lobound_; i < df_upbound_; ++i) {
      bra_idx[0] = i;
      ket_idx[0] = i;
      madness::Future<Tile> future_bra_tile = bra_.find(bra_idx);
      madness::Future<Tile> future_ket_tile = ket_.find(ket_idx);

      reduce_pair_task.add(future_bra_tile, future_ket_tile);
    }

    // permute function
    auto permute = [](const Tile &tile) {
      return tile.permute(TA::Permutation({0, 2, 1, 3}));
    };

    if (notation_ == Formula::Notation::Chemical) {
      return reduce_pair_task.submit();
    } else {
      return world.taskq.add(permute, reduce_pair_task.submit());
    }
  }

 private:
  // left hand size three center integral
  TA::DistArray<Tile, Policy> bra_;
  // right hand size three center integral
  TA::DistArray<Tile, Policy> ket_;
  // low bound for density fitting dimension, should be zero
  std::size_t df_lobound_;
  // up bound for density fitting dimension, should be the max
  std::size_t df_upbound_;
  // madness id for serailization
  madness::uniqueidT id_;
  // notation
  Formula::Notation notation_;
};

/*!
 * \brief Function to make detection of template parameters easier, see
 * IntegralBuilder for details.
 */
template <typename Tile, typename Engine>
std::shared_ptr<IntegralBuilder<Tile, Engine>> make_integral_builder(
    ShrPool<Engine> shr_epool, std::shared_ptr<BasisVector> shr_bases,
    std::shared_ptr<Screener> shr_screen,
    std::function<Tile(TA::TensorD &&)> op,
    std::shared_ptr<const math::PetiteList> plist =
        math::PetiteList::make_trivial()) {
  return std::make_shared<IntegralBuilder<Tile, Engine>>(
      std::move(shr_epool), std::move(shr_bases), std::move(shr_screen),
      std::move(op), std::move(plist));
}

/*!
 * \brief Function to make detection of template parameters easier, see
 * DirectIntegralBuilder for details.
 */
template <typename Tile, typename Engine>
std::shared_ptr<DirectIntegralBuilder<Tile, Engine>>
make_direct_integral_builder(madness::World &world, ShrPool<Engine> shr_epool,
                             std::shared_ptr<BasisVector> shr_bases,
                             std::shared_ptr<Screener> shr_screen,
                             std::function<Tile(TA::TensorD &&)> op,
                             std::shared_ptr<const math::PetiteList> plist =
                                 math::PetiteList::make_trivial()) {
  return std::make_shared<DirectIntegralBuilder<Tile, Engine>>(
      world, std::move(shr_epool), std::move(shr_bases), std::move(shr_screen),
      std::move(op), std::move(plist));
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_BUILDER_H_
