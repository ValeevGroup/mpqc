#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_TWO_CENTER_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_TWO_CENTER_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
class PeriodicTwoCenterBuilder
    : public madness::WorldObject<PeriodicTwoCenterBuilder<Tile, Policy>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using WorldObject_ =
      madness::WorldObject<PeriodicTwoCenterBuilder<Tile, Policy>>;
  using PeriodicTwoCenterBuilder_ = PeriodicTwoCenterBuilder<Tile, Policy>;
  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using BasisVector = std::vector<Basis>;
  using Shell = typename ::mpqc::lcao::gaussian::Shell;
  using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;
  using OperType = ::mpqc::Operator::Type;

  PeriodicTwoCenterBuilder(madness::World &world,
                           std::shared_ptr<const Basis> basis,
                           std::shared_ptr<const UnitCell> unitcell,
                           Vector3d &dcell, Vector3i &R_max, Vector3i &RJ_max,
                           int64_t R_size, int64_t RJ_size)
      : WorldObject_(world),
        basis0_(basis),
        unitcell_(unitcell),
        dcell_(dcell),
        R_max_(R_max),
        RJ_max_(RJ_max),
        R_size_(R_size),
        RJ_size_(RJ_size) {
    assert(basis0_ != nullptr && "No basis is provided");
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    // initialize compound bases
    init();
  }

  ~PeriodicTwoCenterBuilder() { engines_.resize(0); }

  array_type eval(OperType oper, double target_precision =
                                     std::numeric_limits<double>::epsilon()) {
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::detail::shift_mol_origin;
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
    target_precision_ = target_precision;

    // compute significant shell pair list
    //    ExEnv::out0() << "Computing shell pair list...\n" << std::endl;
    //    sig_shellpair_list_ = parallel_compute_shellpair_list(basis0, basisR);

    // # of tiles per basis
    auto ntiles0 = basis0_->nclusters();
    auto ntilesR = basisR_->nclusters();
    const auto ntile_tasks =
        static_cast<uint64_t>(ntiles0 * ntilesR * RJ_size_);
    auto pmap = std::make_shared<const TA::detail::BlockedPmap>(compute_world,
                                                                ntile_tasks);

    for (auto tile0 = 0ul, tile01 = 0ul; tile0 != ntiles0; ++tile0) {
      for (auto tileR = 0ul; tileR != ntilesR; ++tileR) {
        for (auto RJ = 0; RJ != RJ_size_; ++RJ, ++tile01) {
          if (pmap->is_local(tile01))
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

  // mutated by compute_ functions
  mutable std::shared_ptr<Basis> basisR_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> local_result_tiles_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> global_result_tiles_;
  mutable TA::TiledRange result_trange_;
  mutable std::shared_ptr<TA::Pmap> result_pmap_;
  mutable double target_precision_ = 0.0;
  mutable std::vector<Engine> engines_;
  mutable shellpair_list_t sig_shellpair_list_;

  void init() {
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
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

    // 1-d tile ranges
    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);
    const auto ntiles0 = tr0.tile_extent();
    const auto ntilesR = tr1.tile_extent();
    const auto &rng0 = tr0.tile(tile0);
    const auto &rngR = tr1.tile(tileR);
    const auto rng0_size = rng0.second - rng0.first;
    const auto rngR_size = rngR.second - rngR.first;

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng0, rngR});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);
    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();

    {
      const auto engine_precision = target_precision_;
      const auto &basis0 = basis0_;
      const auto &basisR = basisR_;

      // shell clusters for this tile
      const auto &cluster0 = basis0->cluster_shells()[tile0];
      const auto &clusterR = basisR->cluster_shells()[tileR];

      // number of shells in each cluster
      const auto nshells0 = cluster0.size();
      const auto nshellsR = clusterR.size();

      auto engine = engines_[RJ]->local();
      engine.set_precision(engine_precision);
      const auto &computed_shell_sets = engine.results();

      // make non-negligible shell pair list
      auto sig_shellpair_list = compute_shellpair_list(cluster0, clusterR);

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

        for (const auto &shR : sig_shellpair_list[sh0]) {
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

    // accumulate the local contributions
    {
      auto tile0R = tile0 * ntilesR + tileR;
      PeriodicTwoCenterBuilder_::accumulate_local_task(result_tile, tile0R);
    }
  }

  /*!
   * \brief This computes non-negligible shell pair list; ; shells \c i and \c j
   * form a non-negligible pair if they share a center or the Frobenius norm of
   * their overlap isgreater than threshold
   * \param basis1 a basis
   * \param basis2 a basis
   * \param threshold
   *
   * \return a list of pairs with
   * key: shell index
   * mapped value: a vector of shell indices
   */
  shellpair_list_t parallel_compute_shellpair_list(
      const Basis &basis1, const Basis &basis2,
      double threshold = 1e-12) const {
    using ::mpqc::lcao::gaussian::make_engine_pool;
    using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
    // initialize engine
    auto engine_pool = make_engine_pool(
        libint2::Operator::overlap, utility::make_array_of_refs(basis1, basis2),
        libint2::BraKet::x_x);

    auto &world = this->get_world();

    shellpair_list_t result;
    const auto bs1_equiv_bs2 = false;

    const auto &shv1 = basis1.flattened_shells();
    const auto &shv2 = basis2.flattened_shells();
    const auto nsh1 = shv1.size();
    const auto nsh2 = shv2.size();

    for (auto s1 = 0, s12 = 0; s1 != nsh1; ++s1) {
      result.insert(std::make_pair(s1, std::vector<size_t>()));
      auto n1 = shv1[s1].size();

      auto compute = [&](int s1) {

        const auto engine_precision = target_precision_;
        auto engine = engine_pool->local();
        engine.set_precision(engine_precision);
        const auto &buf = engine.results();

        auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for (auto s2 = 0l; s2 <= s2_max; ++s2, ++s12) {
          auto on_same_center = (shv1[s1].O == shv2[s2].O);
          bool significant = on_same_center;
          if (!on_same_center) {
            auto n2 = shv2[s2].size();
            engine.compute1(shv1[s1], shv2[s2]);
            Eigen::Map<const RowMatrixXd> buf_mat(buf[0], n1, n2);
            auto norm = buf_mat.norm();
            significant = (norm >= threshold);
          }

          if (significant) result[s1].emplace_back(s2);
        }

      };

      // world.taskq.add(compute, s1);
      compute(s1);
    }
    world.gop.fence();

    // resort shell list in increasing order
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      auto &list = result[s1];
      std::sort(list.begin(), list.end());
    }

    return result;
  }

  /*!
   * \brief This computes non-negligible shell pair list; ; shells \c i and \c j
   * form a non-negligible pair if they share a center or the Frobenius norm of
   * their overlap isgreater than threshold
   * \param shv1 a cluster (a.k.a. std::vector<Shell>)
   * \param shv2 a cluster (a.k.a. std::vector<Shell>)
   * \param threshold
   *
   * \return a list of pairs with
   * key: shell index
   * mapped value: a vector of shell indices
   */
  shellpair_list_t compute_shellpair_list(
      const ShellVec &shv1,
      const ShellVec &_shv2 = std::vector<Shell>({Shell()}),
      double threshold = 1e-12) const {
    const ShellVec &shv2 =
        ((_shv2.size() == 1 && _shv2[0] == Shell()) ? shv1 : _shv2);
    const auto nsh1 = shv1.size();
    const auto nsh2 = shv2.size();
    const auto shv1_equiv_shv2 = (&shv1 == &shv2);

    // determine max # of primitives in a shell cluster
    auto max_nprim = [](const ShellVec &shv) {
      size_t n = 0;
      for (auto shell : shv) n = std::max(shell.nprim(), n);
      return n;
    };
    const auto max_nprim_1 = max_nprim(shv1);
    const auto max_nprim_2 = max_nprim(shv2);

    // determine max angular momentum of a shell cluster
    auto max_l = [](const ShellVec &shv) {
      int l = 0;
      for (auto shell : shv)
        for (auto c : shell.contr) l = std::max(c.l, l);
      return l;
    };
    const auto max_l_1 = max_l(shv1);
    const auto max_l_2 = max_l(shv2);

    // initialize libint2 engine
    auto engine = libint2::Engine(libint2::Operator::overlap,
                                  std::max(max_nprim_1, max_nprim_2),
                                  std::max(max_l_1, max_l_2), 0);
    const auto &buf = engine.results();
    shellpair_list_t result;

    // compute non-negligible shell-pair list
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      result.insert(std::make_pair(s1, std::vector<size_t>()));
      auto n1 = shv1[s1].size();

      auto s2_max = shv1_equiv_shv2 ? s1 : nsh2 - 1;
      for (auto s2 = 0l; s2 <= s2_max; ++s2, ++s12) {
        auto on_same_center = (shv1[s1].O == shv2[s2].O);
        bool significant = on_same_center;
        if (!on_same_center) {
          auto n2 = shv2[s2].size();
          engine.compute(shv1[s1], shv2[s2]);
          Eigen::Map<const RowMatrixXd> buf_mat(buf[0], n1, n2);
          auto norm = buf_mat.norm();
          significant = (norm >= threshold);
        }

        if (significant) result[s1].emplace_back(s2);
      }
    }

    // resort shell list in increasing order
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      auto &list = result[s1];
      std::sort(list.begin(), list.end());
    }

    return result;
  }

  /*!
   * \brief This computes basis function offsets for every shell in a cluster
   * \param cluster a cluster (a.k.a. std::vector<Shell>)
   * \param bf_first basis function index of the first function in this \c
   * cluster
   *
   * \return a list of <key, mapped value> pairs with
   * key: shell index
   * mapped value: {cluster function offset, basis function offset} tuple
   */
  func_offset_list compute_func_offset_list(const ShellVec &cluster,
                                            const size_t bf_first) const {
    func_offset_list result;

    auto cf_offset = 0;
    auto bf_offset = bf_first;

    const auto nshell = cluster.size();
    for (auto s = 0; s != nshell; ++s) {
      const auto &shell = cluster[s];
      const auto nf = shell.size();
      result.insert(std::make_pair(s, std::make_tuple(cf_offset, bf_offset)));
      bf_offset += nf;
      cf_offset += nf;
    }

    return result;
  }
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_TWO_CENTER_BUILDER_H_
