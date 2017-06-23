#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_THREE_CENTER_CONTRACTION_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_THREE_CENTER_CONTRACTION_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"

#include <mutex>

namespace mpqc {
namespace scf {

/// PeriodicThreeCenterContractionBuilder is an integral-direct implementation
/// of two types of contractions (X|μν) D_μν and (X|μν) X that appear in
/// periodic RI-J method. This builder takes advantage of shell-level screening
/// and parallelized summation over unit cells for coulomb interaction
template <typename Tile, typename Policy>
class PeriodicThreeCenterContractionBuilder
    : public madness::WorldObject<
          PeriodicThreeCenterContractionBuilder<Tile, Policy>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using const_data_ptr = typename Tile::allocator_type::const_pointer;

  using WorldObject_ =
      madness::WorldObject<PeriodicThreeCenterContractionBuilder<Tile, Policy>>;
  using PeriodicThreeCenterContractionBuilder_ =
      PeriodicThreeCenterContractionBuilder<Tile, Policy>;

  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using BasisVector = std::vector<Basis>;
  using Shell = typename ::mpqc::lcao::gaussian::Shell;
  using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;

  PeriodicThreeCenterContractionBuilder(
      madness::World &world, std::shared_ptr<const Basis> basis,
      std::shared_ptr<const Basis> aux_basis, Vector3d &dcell, Vector3i &R_max,
      Vector3i &RJ_max, Vector3i &RD_max, int64_t R_size, int64_t RJ_size,
      int64_t RD_size, std::string screen = "schwarz",
      double screen_threshold = 1.0e-20, double shell_pair_threshold = 1.0e-12)
      : WorldObject_(world),
        basis0_(basis),
        aux_basis_(aux_basis),
        screen_(screen),
        screen_threshold_(screen_threshold),
        dcell_(dcell),
        R_max_(R_max),
        RJ_max_(RJ_max),
        RD_max_(RD_max),
        R_size_(R_size),
        RJ_size_(RJ_size),
        RD_size_(RD_size),
        shell_pair_threshold_(shell_pair_threshold) {
        RD_size_(RD_size) {
    assert(basis0_ != nullptr && "No basis is provided");
    assert(aux_basis_ != nullptr && "No auxiliary basis is provided");
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    // initialize compound bases, engines, and screeners
    init();
  }

  ~PeriodicThreeCenterContractionBuilder() {}

  /*!
   * \brief This computes either (X|μν) D_μν or (X|μν) X contraction.
   * \tparam target_rank is 1 for (X|μν) D_μν and is 2 for (X|μν) X
   * \param M is density matrix or X
   * \param target_precision
   * \return result array should have same rank as \tparam target_rank
   */
  template <size_t target_rank>
  array_type contract_with(array_type const &M, double target_precision) {
    auto arg_trange = M.trange();
    auto arg_elements_range = arg_trange.elements_range();
    auto arg_tiles_range = arg_trange.tiles_range();
    auto arg_rank = arg_trange.rank();
    auto ntiles_obs = basis0_->nclusters();
    auto nbf_obs = basis0_->nfunctions();
    auto ntiles_aux = aux_basis_->nclusters();
    auto nbf_aux = aux_basis_->nfunctions();
    if (target_rank == 1) {
      // validate preconditions
      assert(arg_rank == 2);
      assert(arg_elements_range.extent(0) == nbf_obs &&
             arg_elements_range.extent(1) == (nbf_obs * RD_size_));
      assert(arg_tiles_range.extent(0) == ntiles_obs &&
             arg_tiles_range.extent(1) == (ntiles_obs * RD_size_));
      // compute contraction
      return compute_contr_Xmn_mn(M, target_precision);
    } else if (target_rank == 2) {
      // validate preconditions
      assert(arg_rank == 1);
      assert(arg_elements_range.extent(0) == nbf_aux);
      assert(arg_tiles_range.extent(0) == ntiles_aux);
      // compute contraction
      return compute_contr_Xmn_X(M, target_precision);
    }
  }

  /*!
   * \brief This is the implementation of (X|μν) D_μν contraction
   * \param D density matrix
   * \param target_precision
   * \return
   */
  array_type compute_contr_Xmn_mn(array_type const &D,
                                  double target_precision) const {
    // Copy D and make it replicated.
    array_type D_repl;
    D_repl("i,j") = D("i,j");
    D_repl.make_replicated();
    arg_pmap_repl_ = D_repl.pmap();
    arg_trange_ = D_repl.trange();

    // prepare input data
    auto &compute_world = this->get_world();
    const auto me = compute_world.rank();
    const auto nproc = compute_world.nproc();
    target_precision_ = target_precision;

    // make trange and pmap for the result
    result_trange_ = TA::TiledRange({trange1_aux_});
    result_pmap_ = Policy::default_pmap(compute_world,
                                        result_trange_.tiles_range().volume());

    // # of tiles per basis
    auto ntiles0 = basis0_->nclusters();
    auto ntilesRD = basisRD_->nclusters();
    auto ntiles_aux = aux_basis_->nclusters();

    // make shell block norm of D
    auto shblk_norm_D = compute_shellblock_norm(*basis0_, *basisRD_, D_repl);
    shblk_norm_D.make_replicated();  // make sure it is replicated

    // initialize engines
    {
      using ::mpqc::lcao::gaussian::make_engine_pool;
      auto oper_type = libint2::Operator::coulomb;
      assert(RJ_size_ > 0 && RJ_size_ % 2 == 1);
      auto ref_uc = (RJ_size_ - 1) / 2;
      const auto aux_basis_RJ = *(aux_basis_RJ_[ref_uc]);
      const auto basis0 = *basis0_;
      const auto basisRD = *basisRD_;
      engines_ = make_engine_pool(
          oper_type, utility::make_array_of_refs(aux_basis_RJ, basis0, basisRD),
          libint2::BraKet::xs_xx);
    }

    auto empty = TA::Future<Tile>(Tile());
    for (auto tile_aux = 0ul, tile012 = 0ul; tile_aux != ntiles_aux;
         ++tile_aux) {
      for (auto tile0 = 0ul; tile0 != ntiles0; ++tile0) {
        for (auto tileRD = 0ul; tileRD != ntilesRD; ++tileRD) {
          auto D_0RD = (D_repl.is_zero({tile0, tileRD}))
                          ? empty
                          : D_repl.find({tile0, tileRD});
          auto norm_D_0RD = (shblk_norm_D.is_zero({tile0, tileRD}))
                               ? empty
                               : shblk_norm_D.find({tile0, tileRD});
          if (D_0RD.get().data() == nullptr || norm_D_0RD.get().data() == nullptr)
            continue;
          for (auto RJ = 0; RJ != RJ_size_; ++RJ, ++tile012) {
            if (tile012 % nproc == me)
              WorldObject_::task(
                  me,
                  &PeriodicThreeCenterContractionBuilder_::compute_task_Xmn_mn,
                  D_0RD, norm_D_0RD, RJ,
                  std::array<size_t, 3>{{tile_aux, tile0, tileRD}});
          }
        }
      }
    }

    compute_world.gop.fence();

    // clean up
    engines_.reset();

    // collect local tiles
    if (arg_pmap_repl_->is_replicated() && compute_world.size() > 1) {
      for (const auto &local_tile : local_contr_tiles_) {
        const auto tile_ord = local_tile.first;
        const auto proc = result_pmap_->owner(tile_ord);
        WorldObject_::task(
            proc,
            &PeriodicThreeCenterContractionBuilder_::accumulate_global_task,
            local_tile.second, tile_ord);
      }
      local_contr_tiles_.clear();
      compute_world.gop.fence();

      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 1>, double>> global_tile_norms;
        for (const auto &global_tile : global_contr_tiles_) {
          const auto tile_ord = global_tile.first;
          const auto norm = global_tile.second.norm();
          global_tile_norms.push_back(
              std::make_pair(std::array<size_t, 1>{{tile_ord}}, norm));
        }
        shape =
            decltype(shape)(compute_world, global_tile_norms, result_trange_);
      }

      array_type result(compute_world, result_trange_, shape, result_pmap_);
      for (const auto &global_tile : global_contr_tiles_) {
        if (!result.shape().is_zero(global_tile.first))
          result.set(global_tile.first, global_tile.second);
      }
      result.fill_local(0.0, true);
      global_contr_tiles_.clear();

      return result;
    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 1>, double>> local_tile_norms;
        for (const auto &local_tile : local_contr_tiles_) {
          const auto tile_ord = local_tile.first;
          const auto norm = local_tile.second.norm();
          local_tile_norms.push_back(
              std::make_pair(std::array<size_t, 1>{{tile_ord}}, norm));
        }
        shape =
            decltype(shape)(compute_world, local_tile_norms, result_trange_);
      }

      array_type result(compute_world, result_trange_, shape, result_pmap_);
      for (const auto &local_tile : local_contr_tiles_) {
        if (!result.shape().is_zero(local_tile.first))
          result.set(local_tile.first, local_tile.second);
      }
      result.fill_local(0.0, true);
      local_contr_tiles_.clear();

      return result;
    }
  }

  /*!
   * \brief This is the implementation of (X|μν) X contraction
   * \param X is array X
   * \param target_precision
   * \return
   */
  array_type compute_contr_Xmn_X(array_type const &X,
                                 double target_precision) const {
    // copy X and make it replicated
    array_type X_repl;
    X_repl("X") = X("X");
    X_repl.make_replicated();
    arg_pmap_repl_ = X_repl.pmap();

    // prepare input data
    auto &compute_world = this->get_world();
    const auto me = compute_world.rank();
    const auto nproc = compute_world.nproc();
    target_precision_ = target_precision;

    {
      // make trange and pmap for the result
      const auto basis0 = *basis0_;
      const auto basisR = *basisR_;
      auto result_bases = ::mpqc::lcao::gaussian::BasisVector{{basis0, basisR}};
      result_trange_ =
          ::mpqc::lcao::gaussian::detail::create_trange(result_bases);
      result_pmap_ = Policy::default_pmap(
          compute_world, result_trange_.tiles_range().volume());

      // initialize engines
      using ::mpqc::lcao::gaussian::make_engine_pool;
      auto oper_type = libint2::Operator::coulomb;
      assert(RJ_size_ > 0 && RJ_size_ % 2 == 1);
      auto ref_uc = (RJ_size_ - 1) / 2;
      const auto aux_basis_RJ = *(aux_basis_RJ_[ref_uc]);
      engines_ = make_engine_pool(
          oper_type, utility::make_array_of_refs(aux_basis_RJ, basis0, basisR),
          libint2::BraKet::xs_xx);
    }

    // # of tiles per basis
    auto ntiles0 = basis0_->nclusters();
    auto ntilesR = basisR_->nclusters();
    auto ntiles_aux = aux_basis_->nclusters();

    auto empty = TA::Future<Tile>(Tile());
    for (auto tile_aux = 0ul, tile012 = 0ul; tile_aux != ntiles_aux;
         ++tile_aux) {
      auto X_aux =
          (X_repl.is_zero({tile_aux})) ? empty : X_repl.find({tile_aux});
      if (X_aux.get().data() == nullptr) continue;

      for (auto tile0 = 0ul; tile0 != ntiles0; ++tile0) {
        for (auto tileR = 0ul; tileR != ntilesR; ++tileR) {
          for (auto RJ = 0; RJ != RJ_size_; ++RJ, ++tile012) {
            if (tile012 % nproc == me)
              WorldObject_::task(
                  me,
                  &PeriodicThreeCenterContractionBuilder_::compute_task_Xmn_X,
                  X_aux, RJ, std::array<size_t, 3>{{tile_aux, tile0, tileR}});
          }
        }
      }
    }

    compute_world.gop.fence();

    // clean up
    engines_.reset();

    // collect local tiles
    if (arg_pmap_repl_->is_replicated() && compute_world.size() > 1) {
      for (const auto &local_tile : local_contr_tiles_) {
        const auto tile_ord = local_tile.first;
        const auto proc = result_pmap_->owner(tile_ord);
        WorldObject_::task(
            proc,
            &PeriodicThreeCenterContractionBuilder_::accumulate_global_task,
            local_tile.second, tile_ord);
      }
      local_contr_tiles_.clear();
      compute_world.gop.fence();

      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> global_tile_norms;
        for (const auto &global_tile : global_contr_tiles_) {
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
      for (const auto &global_tile : global_contr_tiles_) {
        if (!result.shape().is_zero(global_tile.first))
          result.set(global_tile.first, global_tile.second);
      }
      result.fill_local(0.0, true);
      global_contr_tiles_.clear();

      return result;
    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> local_tile_norms;
        for (const auto &local_tile : local_contr_tiles_) {
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
      for (const auto &local_tile : local_contr_tiles_) {
        if (!result.shape().is_zero(local_tile.first))
          result.set(local_tile.first, local_tile.second);
      }
      result.fill_local(0.0, true);
      local_contr_tiles_.clear();

      return result;
    }
  }

 private:
  // set by ctor
  std::shared_ptr<const Basis> basis0_;
  std::shared_ptr<const Basis> aux_basis_;
  const std::string screen_;
  const double screen_threshold_;
  const double shell_pair_threshold_;
  const Vector3d dcell_;
  const Vector3i R_max_;
  const Vector3i RJ_max_;
  const Vector3i RD_max_;
  const int64_t R_size_;
  const int64_t RJ_size_;
  const int64_t RD_size_;

  // mutated by compute_ functions
  mutable std::shared_ptr<lcao::Screener> p_screener_R_;
  mutable std::shared_ptr<lcao::Screener> p_screener_RD_;
  mutable std::shared_ptr<Basis> basisR_;
  mutable std::shared_ptr<Basis> basisRD_;
  mutable std::vector<std::shared_ptr<Basis>> aux_basis_RJ_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> local_contr_tiles_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> global_contr_tiles_;
  mutable TA::TiledRange arg_trange_;
  mutable TA::TiledRange1 trange1_aux_;
  mutable TA::TiledRange result_trange_;
  mutable std::shared_ptr<TA::Pmap> result_pmap_;
  mutable std::shared_ptr<TA::Pmap> dist_pmap_D_;
  mutable std::shared_ptr<TA::Pmap> arg_pmap_repl_;
  mutable double target_precision_ = 0.0;
  mutable Engine engines_;
  mutable array_type shblk_norm_D_;
  mutable shellpair_list_t sig_shellpair_list_R_;
  mutable shellpair_list_t sig_shellpair_list_RD_;
  mutable std::unordered_map<size_t, size_t> basis0_shell_offset_map_;
  mutable std::unordered_map<size_t, size_t> basisR_shell_offset_map_;
  mutable std::unordered_map<size_t, size_t> basisRD_shell_offset_map_;

  void init() {
    trange1_aux_ = aux_basis_->create_trange1();

    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::make_engine_pool;
    using ::mpqc::lcao::gaussian::detail::make_screener;

    // make compound basis set for ket1 in (bra | ket0 ket1)
    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    basisR_ = shift_basis_origin(*basis0_, zero_shift_base, R_max_, dcell_);
    basisRD_ = shift_basis_origin(*basis0_, zero_shift_base, RD_max_, dcell_);

    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      // make compound basis set for bra (aux)
      aux_basis_RJ_.emplace_back(shift_basis_origin(*aux_basis_, vec_RJ));
    }

    const auto basis0 = *basis0_;
    const auto basisR = *basisR_;
    const auto basisRD = *basisRD_;
    auto oper_type = libint2::Operator::coulomb;
    auto &world = this->get_world();
    // initialize screener
    if (screen_ == "schwarz") {
      const auto tmp_basis =
          *(shift_basis_origin(*aux_basis_, zero_shift_base, RJ_max_, dcell_));
      // screener involving (X | μ0 νR)
      {
        auto tmp_eng = make_engine_pool(
            oper_type, utility::make_array_of_refs(tmp_basis, basis0, basisR),
            libint2::BraKet::xx_xx);
        auto bases =
            ::mpqc::lcao::gaussian::BasisVector{{tmp_basis, basis0, basisR}};
        p_screener_R_ =
            make_screener(world, tmp_eng, bases, screen_, screen_threshold_);
      }

      // screener involving (X | μ0 ν_Rd)
      if (R_max_ == RD_max_) {
        p_screener_RD_ = p_screener_R_;
      } else {
        auto tmp_eng = make_engine_pool(
            oper_type, utility::make_array_of_refs(tmp_basis, basis0, basisRD),
            libint2::BraKet::xx_xx);
        auto bases =
            ::mpqc::lcao::gaussian::BasisVector{{tmp_basis, basis0, basisRD}};
        p_screener_RD_ =
            make_screener(world, tmp_eng, bases, screen_, screen_threshold_);
      }
    } else {
      throw InputError("Wrong screening method", __FILE__, __LINE__, "screen");
    }

    // compute significant shell pair list
    {
      sig_shellpair_list_R_ = parallel_compute_shellpair_list(
          basis0, basisR, shell_pair_threshold_);
      basis0_shell_offset_map_ = compute_shell_offset(basis0);
      basisR_shell_offset_map_ = compute_shell_offset(basisR);

      if (R_max_ == RD_max_) {
        sig_shellpair_list_RD_ = sig_shellpair_list_R_;
        basisRD_shell_offset_map_ = basisR_shell_offset_map_;
      } else {
        sig_shellpair_list_RD_ = parallel_compute_shellpair_list(
              basis0, basisRD, shell_pair_threshold_);
        basisRD_shell_offset_map_ = compute_shell_offset(basisRD);
      }
    }
  }

  void accumulate_global_task(Tile arg_tile, long tile_ord) {
    // if reducer does not exist, create entry and store F, else accumulate F to
    // the existing contents
    typename decltype(global_contr_tiles_)::accessor acc;
    // try inserting, otherwise, accumulate
    if (!global_contr_tiles_.insert(
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
    typename decltype(local_contr_tiles_)::accessor acc;
    // try inserting, otherwise, accumulate
    if (!local_contr_tiles_.insert(
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

  void compute_task_Xmn_mn(Tile D_0RD, Tile norm_D_0RD, size_t RJ,
                           std::array<size_t, 3> tile_idx) {
    const auto tile_aux = tile_idx[0];
    const auto tile0 = tile_idx[1];
    const auto tileRD = tile_idx[2];

    // get reference to basis sets
    const auto &basisRJ_aux = aux_basis_RJ_[RJ];
    const auto &basis0 = basis0_;
    const auto &basisRD = basisRD_;

    // shell clusters for this tile
    const auto &clusterRJ_aux = basisRJ_aux->cluster_shells()[tile_aux];
    const auto &cluster0 = basis0->cluster_shells()[tile0];
    const auto &clusterRD = basisRD->cluster_shells()[tileRD];

    // number of shells in each cluster
    const auto nshellsRJ_aux = clusterRJ_aux.size();
    const auto nshells0 = cluster0.size();
    const auto nshellsRD = clusterRD.size();

    // 1-d tile ranges
    const auto &tr0 = arg_trange_.dim(0);
    const auto &tr1 = arg_trange_.dim(1);
    const auto &rng_aux = trange1_aux_.tile(tile_aux);
    const auto &rng0 = tr0.tile(tile0);
    const auto &rngRD = tr1.tile(tileRD);
    const auto rngRD_size = rngRD.second - rngRD.first;

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng_aux});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);

    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();
    const auto *D_0RD_ptr = D_0RD.data();
    const auto *norm_D_0RD_ptr = norm_D_0RD.data();
    assert(D_0RD_ptr != nullptr);
    assert(norm_D_0RD_ptr != nullptr);

    const auto nbf_aux_per_uc = basisRJ_aux->nfunctions();

    // compute contributions to all result matrices
    {
      // index of first shell in this cluster
      const auto sh0_offset = basis0_shell_offset_map_[tile0];
      const auto shRD_offset = basisRD_shell_offset_map_[tileRD];

      // index of last shell in this cluster
      const auto sh0_max = sh0_offset + nshells0;
      const auto shRD_max = shRD_offset + nshellsRD;

      // determine if this task contains significant shell-pairs
      auto is_significant = false;
      {
        auto sh0_in_basis = sh0_offset;
        for (; sh0_in_basis != sh0_max; ++sh0_in_basis) {
          for (const auto shRD_in_basis : sig_shellpair_list_RD_[sh0_in_basis]) {
            if (shRD_in_basis >= shRD_offset && shRD_in_basis < shRD_max) {
              is_significant = true;
              break;
            }
          }
          if (is_significant) break;
        }
      }

      if (is_significant) {
        auto &screen = *(p_screener_RD_);
        auto engine = engines_->local();
        const auto engine_precision = target_precision_;
        engine.set_precision(engine_precision);
        const auto &computed_shell_sets = engine.results();

        // compute offset list of clusterRD (ket1)
        auto offset_list_ket1 = compute_func_offset_list(clusterRD, rngRD.first);

        // this is the index of the first basis functions for each shell *in
        // this shell cluster*
        auto cf_aux_offset = 0;
        // this is the index of the first basis functions for each shell *in the
        // basis set*
        auto bf_aux_offset = rng_aux.first;

        size_t cf_RD_offset, bf_RD_offset;

        // loop over all shell sets
        for (auto sh_aux = 0; sh_aux != nshellsRJ_aux; ++sh_aux) {
          const auto &shell_aux = clusterRJ_aux[sh_aux];
          const auto nf_aux = shell_aux.size();

          const auto bf_aux_in_screener = bf_aux_offset + nbf_aux_per_uc * RJ;

          auto cf0_offset = 0;
          auto bf0_offset = rng0.first;
          for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
            const auto &shell0 = cluster0[sh0];
            const auto nf0 = shell0.size();

            const auto sh0_in_basis = sh0 + sh0_offset;
            for (const auto &shRD_in_basis : sig_shellpair_list_RD_[sh0_in_basis]) {
              if (shRD_in_basis < shRD_offset || shRD_in_basis >= shRD_max)
                continue;

              const auto shRD = shRD_in_basis - shRD_offset;
              std::tie(cf_RD_offset, bf_RD_offset) = offset_list_ket1[shRD];

              const auto &shellRD = clusterRD[shRD];
              const auto nfRD = shellRD.size();

              const auto sh0RD =
                  sh0 * nshellsRD + shRD;  // index of {sh0, shRD} in norm_D0RD
              const double Dnorm0RD = norm_D_0RD_ptr[sh0RD];

              if (screen.skip(bf_aux_in_screener, bf0_offset, bf_RD_offset,
                              Dnorm0RD))
                continue;

              // compute shell set
              // TODO call 3-body version of compute2 to avoid extra copies
              engine.compute(shell_aux, shell0, shellRD);
              const auto *eri_aux_0RD = computed_shell_sets[0];

              if (eri_aux_0RD != nullptr) {
                for (auto f_aux = 0, f_aux_0RD = 0; f_aux != nf_aux; ++f_aux) {
                  const auto cf_aux = f_aux + cf_aux_offset;
                  for (auto f0 = 0; f0 != nf0; ++f0) {
                    const auto cf0 = f0 + cf0_offset;
                    for (auto fRD = 0; fRD != nfRD; ++fRD, ++f_aux_0RD) {
                      const auto cfRD = fRD + cf_RD_offset;
                      const auto cf0RD = cf0 * rngRD_size + cfRD;

                      const auto value = eri_aux_0RD[f_aux_0RD];

                      result_ptr[cf_aux] += D_0RD_ptr[cf0RD] * value;
                    }
                  }
                }
              }
            }

            cf0_offset += nf0;
            bf0_offset += nf0;
          }

          cf_aux_offset += nf_aux;
          bf_aux_offset += nf_aux;
        }
      }
    }

    // accumulate the local contributions
    {
      PeriodicThreeCenterContractionBuilder_::accumulate_local_task(result_tile,
                                                                    tile_aux);
    }
  }

  void compute_task_Xmn_X(Tile X_aux, size_t RJ,
                          std::array<size_t, 3> tile_idx) {
    const auto tile_aux = tile_idx[0];
    const auto tile0 = tile_idx[1];
    const auto tileR = tile_idx[2];

    // get reference to basis sets
    const auto &basisRJ_aux = aux_basis_RJ_[RJ];
    const auto &basis0 = basis0_;
    const auto &basisR = basisR_;

    // shell clusters for this tile
    const auto &clusterRJ_aux = basisRJ_aux->cluster_shells()[tile_aux];
    const auto &cluster0 = basis0->cluster_shells()[tile0];
    const auto &clusterR = basisR->cluster_shells()[tileR];

    // number of shells in each cluster
    const auto nshellsRJ_aux = clusterRJ_aux.size();
    const auto nshells0 = cluster0.size();
    const auto nshellsR = clusterR.size();

    // 1-d tile ranges
    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);
    const auto ntilesR = tr1.tile_extent();
    const auto &rng_aux = trange1_aux_.tile(tile_aux);
    const auto &rng0 = tr0.tile(tile0);
    const auto &rngR = tr1.tile(tileR);
    const auto rngR_size = rngR.second - rngR.first;

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng0, rngR});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);

    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();
    const auto *X_aux_ptr = X_aux.data();
    assert(X_aux_ptr != nullptr);

    const auto nbf_aux_per_uc = basisRJ_aux->nfunctions();

    // TODO: check the sparsity of X. Will screening (X|mu) M_X be useful?
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
          for (const auto shR_in_basis : sig_shellpair_list_R_[sh0_in_basis]) {
            if (shR_in_basis >= shR_offset && shR_in_basis < shR_max) {
              is_significant = true;
              break;
            }
          }
          if (is_significant) break;
        }
      }

      if (is_significant) {
        auto &screen = *(p_screener_R_);
        auto engine = engines_->local();
        const auto engine_precision = target_precision_;
        engine.set_precision(engine_precision);
        const auto &computed_shell_sets = engine.results();

        // compute offset list of clusterR (ket1)
        auto offset_list_ket1 = compute_func_offset_list(clusterR, rngR.first);

        // this is the index of the first basis functions for each shell *in
        // this shell cluster*
        auto cf_aux_offset = 0;
        // this is the index of the first basis functions for each shell *in the
        // basis set*
        auto bf_aux_offset = rng_aux.first;

        size_t cf_R_offset, bf_R_offset;

        // loop over all shell sets
        for (auto sh_aux = 0; sh_aux != nshellsRJ_aux; ++sh_aux) {
          const auto &shell_aux = clusterRJ_aux[sh_aux];
          const auto nf_aux = shell_aux.size();

          const auto bf_aux_in_screener = bf_aux_offset + nbf_aux_per_uc * RJ;

          auto cf0_offset = 0;
          auto bf0_offset = rng0.first;
          for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
            const auto &shell0 = cluster0[sh0];
            const auto nf0 = shell0.size();

            const auto sh0_in_basis = sh0 + sh0_offset;
            for (const auto &shR_in_basis : sig_shellpair_list_R_[sh0_in_basis]) {
              if (shR_in_basis < shR_offset || shR_in_basis >= shR_max)
                continue;

              const auto shR = shR_in_basis - shR_offset;
              std::tie(cf_R_offset, bf_R_offset) = offset_list_ket1[shR];

              const auto &shellR = clusterR[shR];
              const auto nfR = shellR.size();

              if (screen.skip(bf_aux_in_screener, bf0_offset, bf_R_offset))
                continue;

              // compute shell set
              // TODO call 3-body version of compute2 to avoid extra copies
              engine.compute(shell_aux, shell0, shellR);

              const auto *eri_aux_0R = computed_shell_sets[0];

              if (eri_aux_0R != nullptr) {
                for (auto f_aux = 0, f_aux_0R = 0; f_aux != nf_aux; ++f_aux) {
                  const auto cf_aux = f_aux + cf_aux_offset;
                  for (auto f0 = 0; f0 != nf0; ++f0) {
                    const auto cf0 = f0 + cf0_offset;
                    for (auto fR = 0; fR != nfR; ++fR, ++f_aux_0R) {
                      const auto cfR = fR + cf_R_offset;
                      const auto cf0R = cf0 * rngR_size + cfR;

                      const auto value = eri_aux_0R[f_aux_0R];

                      result_ptr[cf0R] += X_aux_ptr[cf_aux] * value;
                    }
                  }
                }
              }
            }

            cf0_offset += nf0;
            bf0_offset += nf0;
          }

          cf_aux_offset += nf_aux;
          bf_aux_offset += nf_aux;
        }
      }
    }

    // accumulate the local contributions
    {
      auto tile0R = tile0 * ntilesR + tileR;
      PeriodicThreeCenterContractionBuilder_::accumulate_local_task(result_tile,
                                                                    tile0R);
    }
  }

  /*!
   * \brief This computes shell-block norm of density matrix \c D
   * \param bs Basis
   * \param D density matrix
   * \return
   */
  array_type compute_shellblock_norm(const Basis &bs0, const Basis &bs1,
                                     const array_type &D) const {
    auto &world = this->get_world();
    // make trange1
    auto make_shblk_trange1 = [](const Basis &bs) {
      const auto &shells_Vec = bs.cluster_shells();
      auto blocking = std::vector<int64_t>{0};
      for (const auto &shells : shells_Vec) {
        const auto nshell = shells.size();
        auto next = blocking.back() + nshell;
        blocking.emplace_back(next);
      }
      return TA::TiledRange1(blocking.begin(), blocking.end());
    };

    const auto tr0 = make_shblk_trange1(bs0);
    const auto tr1 = make_shblk_trange1(bs1);

    auto eig_D = ::mpqc::array_ops::array_to_eigen(D);
    // compute shell block norms
    const auto shells0 = bs0.flattened_shells();
    const auto shells1 = bs1.flattened_shells();
    const auto nshell0 = shells0.size();
    const auto nshell1 = shells1.size();
    RowMatrixXd norm_D(nshell0, nshell1);
    for (auto sh0 = 0, sh0_first = 0; sh0 != nshell0; ++sh0) {
      const auto sh0_size = shells0[sh0].size();
      for (auto sh1 = 0, sh1_first = 0; sh1 != nshell1; ++sh1) {
        const auto sh1_size = shells1[sh1].size();

        norm_D(sh0, sh1) = eig_D.block(sh0_first, sh1_first, sh0_size, sh1_size)
                               .template lpNorm<Eigen::Infinity>();

        sh1_first += sh1_size;
      }

      sh0_first += sh0_size;
    }

    return array_ops::eigen_to_array<Tile, Policy>(world, norm_D, tr0, tr1);
  }

  /*!
   * \brief This computes non-negligible shell pair list; ; shells \c i and \c j
   * form a non-negligible pair if they share a center or the Frobenius norm of
   * their overlap is greater than threshold
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
    std::mutex mx;
    shellpair_list_t result;

    const auto &shv1 = basis1.flattened_shells();
    const auto &shv2 = basis2.flattened_shells();
    const auto nsh1 = shv1.size();
    const auto nsh2 = shv2.size();

    auto compute = [&](int64_t input_s1) {

      auto n1 = shv1[input_s1].size();
      const auto engine_precision = target_precision_;
      auto engine = engine_pool->local();
      engine.set_precision(engine_precision);
      const auto &buf = engine.results();

      for (auto s2 = 0l; s2 != nsh2; ++s2) {
        auto on_same_center = (shv1[input_s1].O == shv2[s2].O);
        bool significant = on_same_center;
        if (!on_same_center) {
          auto n2 = shv2[s2].size();
          engine.compute1(shv1[input_s1], shv2[s2]);
          Eigen::Map<const RowMatrixXd> buf_mat(buf[0], n1, n2);
          auto norm = buf_mat.norm();
          significant = (norm >= threshold);
        }

        if (significant) {
          mx.lock();
          result[input_s1].emplace_back(s2);
          mx.unlock();
        }
      }
    };

    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      result.insert(std::make_pair(s1, std::vector<size_t>()));
      world.taskq.add(compute, s1);
    }
    world.gop.fence();

    engine_pool.reset();

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

  /*!
   * \brief This computes shell offsets for every cluster in a basis
   * \param basis
   * \return a list of <key, mapped value> pairs with
   * key: cluster index
   * mapped value: index of first shell in a cluster
   */
  std::unordered_map<size_t, size_t> compute_shell_offset(
      const Basis &basis) const {
    std::unordered_map<size_t, size_t> result;

    auto shell_offset = 0;
    const auto &cluster_shells = basis.cluster_shells();
    const auto nclusters = cluster_shells.size();
    for (auto c = 0; c != nclusters; ++c) {
      const auto nshells = cluster_shells[c].size();
      result.insert(std::make_pair(c, shell_offset));
      shell_offset += nshells;
    }

    return result;
  }
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_THREE_CENTER_CONTRACTION_BUILDER_H_
