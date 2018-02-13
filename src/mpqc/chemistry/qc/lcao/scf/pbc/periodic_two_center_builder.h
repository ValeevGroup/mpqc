#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_TWO_CENTER_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_TWO_CENTER_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"

#include <mutex>

namespace mpqc {
namespace scf {

/// PeriodicTwoCenterBuilder computes 2-center integral in zRHF that takes
/// advantage of parallelized summation over unit cells for coulomb interaction
template <typename Tile, typename Policy>
class PeriodicTwoCenterBuilder
    : public madness::WorldObject<PeriodicTwoCenterBuilder<Tile, Policy>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using WorldObject_ =
      madness::WorldObject<PeriodicTwoCenterBuilder<Tile, Policy>>;
  using PeriodicTwoCenterBuilder_ = PeriodicTwoCenterBuilder<Tile, Policy>;
  using Factory = ::mpqc::lcao::gaussian::PeriodicAOFactory<Tile, Policy>;
  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using BasisVector = std::vector<Basis>;
  using Shell = typename ::mpqc::lcao::gaussian::Shell;
  using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::vector<std::vector<size_t>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;
  using OperType = ::mpqc::Operator::Type;

  PeriodicTwoCenterBuilder(madness::World &world,
                           std::shared_ptr<const Basis> basis,
                           std::shared_ptr<const UnitCell> unitcell,
                           Vector3d &dcell, Vector3i &R_max, Vector3i &RJ_max,
                           int64_t R_size, int64_t RJ_size,
                           shellpair_list_t &sig_shellpair_list)
      : WorldObject_(world),
        basis0_(basis),
        unitcell_(unitcell),
        dcell_(dcell),
        R_max_(R_max),
        RJ_max_(RJ_max),
        R_size_(R_size),
        RJ_size_(RJ_size),
        sig_shellpair_list_(sig_shellpair_list) {
    assert(basis0_ != nullptr && "No basis is provided");
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    // initialize compound bases
    init();
  }

  PeriodicTwoCenterBuilder(Factory &ao_factory)
      : WorldObject_(ao_factory.world()),
        basis0_(ao_factory.basis_registry()->retrieve(OrbitalIndex(L"λ"))),
        unitcell_(std::make_shared<const UnitCell>(ao_factory.unitcell())),
        dcell_(ao_factory.unitcell().dcell()),
        R_max_(ao_factory.R_max()),
        RJ_max_(ao_factory.RJ_max()),
        R_size_(ao_factory.R_size()),
        RJ_size_(ao_factory.RJ_size()),
        sig_shellpair_list_(ao_factory.significant_shell_pairs()) {
    assert(basis0_ != nullptr && "No basis is provided");
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    // initialize compound bases
    init();
  }

  ~PeriodicTwoCenterBuilder() { engines_.resize(0); }

  array_type eval(OperType oper, double target_precision =
                                     std::numeric_limits<double>::epsilon()) {
    using ::mpqc::detail::direct_vector;
    using ::mpqc::detail::shift_mol_origin;
    using ::mpqc::lcao::gaussian::make_engine_pool;
    using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
    using ::mpqc::lcao::gaussian::detail::to_libint2_operator_params;

    // initialize engines
    const auto basis0 = *basis0_;
    const auto basisR = *basisR_;
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      if (oper == OperType::Nuclear) {
        auto shifted_mol = shift_mol_origin(*unitcell_, vec_RJ);
        engines_.emplace_back(make_engine_pool(
            to_libint2_operator(oper),
            utility::make_array_of_refs(basis0, basisR), libint2::BraKet::x_x,
            to_libint2_operator_params(oper, *shifted_mol)));
      } else {
        throw ProgrammingError("Operator type not supported", __FILE__,
                               __LINE__);
      }
    }

    // prepare input data
    auto &compute_world = this->get_world();
    const auto me = compute_world.rank();
    const auto nproc = compute_world.nproc();
    target_precision_ = target_precision;

    // # of tiles per basis
    auto ntiles0 = basis0_->nclusters();
    auto ntilesR = basisR_->nclusters();

    for (auto tile0 = 0ul, tile01 = 0ul; tile0 != ntiles0; ++tile0) {
      for (auto tileR = 0ul; tileR != ntilesR; ++tileR) {
        for (auto RJ = 0; RJ != RJ_size_; ++RJ, ++tile01) {
          if (tile01 % nproc == me)
            WorldObject_::task(me, &PeriodicTwoCenterBuilder_::compute_task_ab,
                               RJ, std::array<size_t, 2>{{tile0, tileR}});
        }
      }
    }

    compute_world.gop.fence();

    // collect local tiles
    if (compute_world.size() > 1) {
      for (const auto &local_tile : local_result_tiles_) {
        const auto tile_ord = local_tile.first;
        const auto proc = result_pmap_->owner(tile_ord);
        WorldObject_::task(proc,
                           &PeriodicTwoCenterBuilder_::accumulate_global_task,
                           local_tile.second, tile_ord);
      }
      local_result_tiles_.clear();
      compute_world.gop.fence();

      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> global_tile_norms;
        for (const auto &global_tile : global_result_tiles_) {
          const auto tile_ord = global_tile.first;
          const auto i = tile_ord / ntilesR;
          const auto j = tile_ord % ntilesR;
          const auto norm = global_tile.second.norm();
          global_tile_norms.push_back(
              std::make_pair(std::array<size_t, 2>{{i, j}}, norm));
        }
        shape =
            decltype(shape)(compute_world, global_tile_norms, result_trange_);
      }

      array_type result(compute_world, result_trange_, shape, result_pmap_);
      for (const auto &global_tile : global_result_tiles_) {
        if (!result.shape().is_zero(global_tile.first))
          result.set(global_tile.first, global_tile.second);
      }
      result.fill_local(0.0, true);
      global_result_tiles_.clear();

      return result;
    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> local_tile_norms;
        for (const auto &local_tile : local_result_tiles_) {
          const auto tile_ord = local_tile.first;
          const auto i = tile_ord / ntilesR;
          const auto j = tile_ord % ntilesR;
          const auto norm = local_tile.second.norm();
          local_tile_norms.push_back(
              std::make_pair(std::array<size_t, 2>{{i, j}}, norm));
        }
        shape =
            decltype(shape)(compute_world, local_tile_norms, result_trange_);
      }

      array_type result(compute_world, result_trange_, shape, result_pmap_);
      for (const auto &local_tile : local_result_tiles_) {
        if (!result.shape().is_zero(local_tile.first))
          result.set(local_tile.first, local_tile.second);
      }
      result.fill_local(0.0, true);
      local_result_tiles_.clear();

      return result;
    }
  }

 private:
  // set by ctor
  std::shared_ptr<const Basis> basis0_;
  std::shared_ptr<const UnitCell> unitcell_;
  const Vector3d dcell_;
  const Vector3i R_max_;
  const Vector3i RJ_max_;
  const int64_t R_size_;
  const int64_t RJ_size_;
  const shellpair_list_t sig_shellpair_list_;

  // mutated by compute_ functions
  mutable std::shared_ptr<Basis> basisR_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> local_result_tiles_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> global_result_tiles_;
  mutable TA::TiledRange result_trange_;
  mutable std::shared_ptr<TA::Pmap> result_pmap_;
  mutable double target_precision_ = 0.0;
  mutable std::vector<Engine> engines_;
  mutable std::unordered_map<size_t, size_t> basis0_shell_offset_map_;
  mutable std::unordered_map<size_t, size_t> basisR_shell_offset_map_;

  void init() {
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::detail::compute_shell_offset;

    // make compound basis set for ket
    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    basisR_ = shift_basis_origin(*basis0_, zero_shift_base, R_max_, dcell_);

    const auto basis0 = *basis0_;
    const auto basisR = *basisR_;

    auto &world = this->get_world();
    // make trange and pmap for the result
    {
      auto result_bases = ::mpqc::lcao::gaussian::BasisVector{{basis0, basisR}};
      result_trange_ =
          ::mpqc::lcao::gaussian::detail::create_trange(result_bases);
      result_pmap_ =
          Policy::default_pmap(world, result_trange_.tiles_range().volume());
    }

    // compute basis shell offsets and cluster shell offsets
    basis0_shell_offset_map_ = compute_shell_offset(basis0);
    basisR_shell_offset_map_ = compute_shell_offset(basisR);
  }

  void accumulate_global_task(Tile arg_tile, long tile_ord) {
    // if reducer does not exist, create entry and store F, else accumulate F to
    // the existing contents
    typename decltype(global_result_tiles_)::accessor acc;
    // try inserting, otherwise, accumulate
    if (!global_result_tiles_.insert(
            acc, std::make_pair(tile_ord, arg_tile))) {  // CRITICAL SECTION
      // NB can't do acc->second += fock_matrix_tile to avoid spawning TBB
      // tasks from critical section
      const auto size = arg_tile.range().volume();
      TA::math::inplace_vector_op_serial(
          [](TA::detail::numeric_t<Tile> &l,
             const TA::detail::numeric_t<Tile> r) { l += r; },
          size, acc->second.data(), arg_tile.data());
    }
    acc.release();  // END OF CRITICAL SECTION
  }

  void accumulate_local_task(Tile arg_tile, long tile_ord) {
    // if reducer does not exist, create entry and store F, else accumulate F to
    // the existing contents
    typename decltype(local_result_tiles_)::accessor acc;
    // try inserting, otherwise, accumulate
    if (!local_result_tiles_.insert(
            acc, std::make_pair(tile_ord, arg_tile))) {  // CRITICAL SECTION
      // NB can't do acc->second += target_tile to avoid spawning TBB
      // tasks from critical section
      const auto size = arg_tile.range().volume();
      TA::math::inplace_vector_op_serial(
          [](TA::detail::numeric_t<Tile> &l,
             const TA::detail::numeric_t<Tile> r) { l += r; },
          size, acc->second.data(), arg_tile.data());
    }
    acc.release();  // END OF CRITICAL SECTION
  }

  void compute_task_ab(size_t RJ, std::array<size_t, 2> tile_idx) {
    const auto tile0 = tile_idx[0];
    const auto tileR = tile_idx[1];

    // get reference to basis sets
    const auto &basis0 = basis0_;
    const auto &basisR = basisR_;

    // shell clusters for this tile
    const auto &cluster0 = basis0->cluster_shells()[tile0];
    const auto &clusterR = basisR->cluster_shells()[tileR];

    // number of shells in each cluster
    const auto nshells0 = cluster0.size();
    const auto nshellsR = clusterR.size();

    // 1-d tile ranges
    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);
    const auto ntilesR = tr1.tile_extent();
    const auto &rng0 = tr0.tile(tile0);
    const auto &rngR = tr1.tile(tileR);
    const auto rngR_size = rngR.second - rngR.first;

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng0, rngR});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);
    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();

    using ::mpqc::lcao::gaussian::detail::compute_func_offset_list;

    {
      // index of first shell in this cluster
      const auto sh0_offset = basis0_shell_offset_map_[tile0];
      const auto shR_offset = basisR_shell_offset_map_[tileR];

      // index of last shell in this cluster
      const auto sh0_max = sh0_offset + nshells0;
      const auto shR_max = shR_offset + nshellsR;

      // determine if this task contains significant shell-pairs
      auto is_significant = false;
      {
        auto sh0_in_basis = sh0_offset;
        for (; sh0_in_basis != sh0_max; ++sh0_in_basis) {
          for (const auto shR_in_basis : sig_shellpair_list_[sh0_in_basis]) {
            if (shR_in_basis >= shR_offset && shR_in_basis < shR_max) {
              is_significant = true;
              break;
            }
          }
          if (is_significant) break;
        }
      }

      if (is_significant) {
        auto engine = engines_[RJ]->local();
        const auto engine_precision = target_precision_;
        engine.set_precision(engine_precision);
        const auto &computed_shell_sets = engine.results();

        // compute offset list of cluster1 and cluster3
        auto offset_list = compute_func_offset_list(clusterR, rngR.first);

        // this is the index of the first basis functions for each shell *in
        // this shell cluster*
        auto cf0_offset = 0;
        // this is the index of the first basis functions for each shell *in the
        // basis set*
        auto bf0_offset = rng0.first;

        size_t cfR_offset, bfR_offset;
        // loop over all shell sets
        for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
          const auto &shell0 = cluster0[sh0];
          const auto nf0 = shell0.size();

          const auto sh0_in_basis = sh0 + sh0_offset;
          for (const auto &shR_in_basis : sig_shellpair_list_[sh0_in_basis]) {
            if (shR_in_basis < shR_offset || shR_in_basis >= shR_max) continue;

            const auto shR = shR_in_basis - shR_offset;
            std::tie(cfR_offset, bfR_offset) = offset_list[shR];

            const auto &shellR = clusterR[shR];
            const auto nfR = shellR.size();

            // TODO add screening for 2-center ints
            engine.compute1(shell0, shellR);

            const auto *ints = computed_shell_sets[0];
            if (ints != nullptr) {
              for (auto f0 = 0, f0R = 0; f0 != nf0; ++f0) {
                const auto cf0 = f0 + cf0_offset;
                for (auto fR = 0; fR != nfR; ++fR, ++f0R) {
                  const auto cfR = fR + cfR_offset;
                  const auto cf0R = cf0 * rngR_size + cfR;

                  const auto value = ints[f0R];
                  result_ptr[cf0R] += value;
                }
              }
            }
          }

          cf0_offset += nf0;
          bf0_offset += nf0;
        }
      }
    }

    // accumulate the local contributions
    {
      auto tile0R = tile0 * ntilesR + tileR;
      PeriodicTwoCenterBuilder_::accumulate_local_task(result_tile, tile0R);
    }
  }
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_TWO_CENTER_BUILDER_H_
