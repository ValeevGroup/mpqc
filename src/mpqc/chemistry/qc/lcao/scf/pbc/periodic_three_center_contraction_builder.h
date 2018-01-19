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
  using Factory = ::mpqc::lcao::gaussian::PeriodicAOFactory<Tile, Policy>;

  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using BasisVector = std::vector<Basis>;
  using Shell = typename ::mpqc::lcao::gaussian::Shell;
  using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::vector<std::vector<size_t>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;

  PeriodicThreeCenterContractionBuilder(
      madness::World &world, std::shared_ptr<const Basis> basis,
      std::shared_ptr<const Basis> aux_basis, Vector3d &dcell, Vector3i &R_max,
      Vector3i &RJ_max, Vector3i &RD_max, int64_t R_size, int64_t RJ_size,
      int64_t RD_size, shellpair_list_t &sig_shellpair_list,
      std::string screen = "schwarz", double screen_threshold = 1.0e-20,
      double density_threshold = Policy::shape_type::threshold())
      : WorldObject_(world),
        basis0_(basis),
        aux_basis_(aux_basis),
        screen_(screen),
        screen_threshold_(screen_threshold),
        density_threshold_(density_threshold),
        dcell_(dcell),
        R_max_(R_max),
        RJ_max_(RJ_max),
        RD_max_(RD_max),
        R_size_(R_size),
        RJ_size_(RJ_size),
        RD_size_(RD_size),
        sig_shellpair_list_(sig_shellpair_list) {
    assert(basis0_ != nullptr && "No basis is provided");
    assert(aux_basis_ != nullptr && "No auxiliary basis is provided");
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    // initialize compound bases, engines, and screeners
    init();
  }

  PeriodicThreeCenterContractionBuilder(Factory &ao_factory)
      : WorldObject_(ao_factory.world()),
        basis0_(ao_factory.basis_registry()->retrieve(OrbitalIndex(L"λ"))),
        aux_basis_(ao_factory.basis_registry()->retrieve(OrbitalIndex(L"Κ"))),
        screen_(ao_factory.screen()),
        screen_threshold_(ao_factory.screen_threshold()),
        density_threshold_(ao_factory.density_threshold()),
        dcell_(ao_factory.unitcell().dcell()),
        R_max_(ao_factory.R_max()),
        RJ_max_(ao_factory.RJ_max()),
        RD_max_(ao_factory.RD_max()),
        R_size_(ao_factory.R_size()),
        RJ_size_(ao_factory.RJ_size()),
        RD_size_(ao_factory.RD_size()),
        sig_shellpair_list_(ao_factory.significant_shell_pairs()) {
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
    auto &compute_world = this->get_world();
    const auto me = compute_world.rank();
    const auto nproc = compute_world.nproc();
    target_precision_ = target_precision;

    // Copy D and make it replicated.
    array_type D_repl;
    D_repl("i,j") = D("i,j");
    D_repl.make_replicated();
    compute_world.gop.fence();

    // make trange and pmap for the result
    result_trange_ = TA::TiledRange({trange1_aux_});
    result_pmap_ = Policy::default_pmap(compute_world,
                                        result_trange_.tiles_range().volume());

    // # of tiles per basis
    auto ntiles = basis0_->nclusters();

    // make shell block norm of D
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::detail::compute_shellblock_norm;
    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    auto basis1 =
        shift_basis_origin(*basis0_, zero_shift_base, RD_max_, dcell_);
    auto shblk_norm_D = compute_shellblock_norm(*basis0_, *basis1, D_repl);
    shblk_norm_D.make_replicated();  // make sure it is replicated
    compute_world.gop.fence();

    // initialize engines
    {
      using ::mpqc::lcao::gaussian::make_engine_pool;
      auto oper_type = libint2::Operator::coulomb;
      const auto aux_basis = *aux_basis_;
      const auto basis0 = *basis0_;
      const auto basisR = *basisR_;
      engines_ = make_engine_pool(
          oper_type, utility::make_array_of_refs(aux_basis, basis0, basisR),
          libint2::BraKet::xs_xx);
    }

    using ::mpqc::detail::direct_3D_idx;
    using ::mpqc::detail::direct_ord_idx;
    using ::mpqc::detail::is_in_lattice_range;

    const auto &Dnorm = D_repl.shape().data();
    auto task_id = 0ul;
    for (auto R1_ord = 0; R1_ord != R_size_; ++R1_ord) {
      const auto R1_3D = direct_3D_idx(R1_ord, R_max_);
      if (!is_in_lattice_range(R1_3D, RD_max_)) {
        continue;
      }

      const auto uc_ord_D01 = direct_ord_idx(R1_3D, RD_max_);
      for (auto tile0 = 0ul; tile0 != ntiles; ++tile0) {
        for (auto tile1 = 0ul; tile1 != ntiles; ++tile1) {
          const std::array<long, 2> idx_D01{
              {long(tile0), long(tile1 + uc_ord_D01 * ntiles_per_uc_)}};
          if (Dnorm(idx_D01) < density_threshold_ || D_repl.is_zero(idx_D01) ||
              shblk_norm_D.is_zero(idx_D01)) {
            continue;
          }

          auto D01 = D_repl.find(idx_D01);
          auto norm_D01 = shblk_norm_D.find(idx_D01);
          for (auto RJ_ord = 0; RJ_ord != RJ_size_; ++RJ_ord) {
            for (auto tile_aux = 0ul; tile_aux != ntiles; ++tile_aux) {
              if (task_id % nproc == me) {
                WorldObject_::task(
                    me,
                    &PeriodicThreeCenterContractionBuilder_::
                        compute_task_Xmn_mn,
                    D01, norm_D01,
                    std::array<size_t, 2>{{size_t(R1_ord), size_t(RJ_ord)}},
                    std::array<size_t, 3>{{tile_aux, tile0, tile1}});
              }
              task_id++;
            }
          }
        }
      }
    }

    compute_world.gop.fence();

    // clean up
    engines_.reset();

    // collect local tiles
    if (compute_world.size() > 1) {
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
    auto &compute_world = this->get_world();
    const auto me = compute_world.rank();
    const auto nproc = compute_world.nproc();
    target_precision_ = target_precision;

    // copy X and make it replicated
    array_type X_repl;
    X_repl("X") = X("X");
    X_repl.make_replicated();
    compute_world.gop.fence();

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
      const auto aux_basis = *aux_basis_;
      engines_ = make_engine_pool(
          oper_type, utility::make_array_of_refs(aux_basis, basis0, basisR),
          libint2::BraKet::xs_xx);
    }

    // # of tiles per basis
    auto ntiles_aux = aux_basis_->nclusters();
    auto ntiles0 = basis0_->nclusters();
    auto ntiles1 = basisR_->nclusters();

    auto task_id = 0ul;
    for (auto tile_aux = 0ul; tile_aux != ntiles_aux; ++tile_aux) {
      if (X_repl.is_zero({tile_aux})) {
        continue;
      }

      auto X_aux = X_repl.find({tile_aux});
      for (auto tile0 = 0ul; tile0 != ntiles0; ++tile0) {
        for (auto tileR = 0ul; tileR != ntiles1; ++tileR) {
          for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
            if (task_id % nproc == me) {
              WorldObject_::task(
                  me,
                  &PeriodicThreeCenterContractionBuilder_::compute_task_Xmn_X,
                  X_aux, RJ, std::array<size_t, 3>{{tile_aux, tile0, tileR}});
            }
            task_id++;
          }
        }
      }
    }

    compute_world.gop.fence();

    // clean up
    engines_.reset();

    // collect local tiles
    if (compute_world.size() > 1) {
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
          const auto i = tile_ord / ntiles1;
          const auto j = tile_ord % ntiles1;
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
          const auto i = tile_ord / ntiles1;
          const auto j = tile_ord % ntiles1;
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
  const double density_threshold_;
  const Vector3d dcell_;
  const Vector3i R_max_;
  const Vector3i RJ_max_;
  const Vector3i RD_max_;
  const int64_t R_size_;
  const int64_t RJ_size_;
  const int64_t RD_size_;
  const shellpair_list_t &sig_shellpair_list_;

  // mutable by init function
  mutable size_t ntiles_per_uc_;
  mutable TA::TiledRange1 trange1_aux_;
  mutable TA::TiledRange trange_eri3_;
  mutable std::shared_ptr<lcao::Screener> p_screener_;
  mutable std::shared_ptr<Basis> basisR_;
  mutable std::vector<std::shared_ptr<Basis>> aux_basis_RJ_;
  mutable std::unordered_map<size_t, size_t> basis0_shell_offset_map_;
  mutable std::unordered_map<size_t, size_t> basisR_shell_offset_map_;

  // mutated by compute_ functions
  mutable madness::ConcurrentHashMap<std::size_t, Tile> local_contr_tiles_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile> global_contr_tiles_;
  mutable TA::TiledRange result_trange_;
  mutable std::shared_ptr<TA::Pmap> result_pmap_;
  mutable double target_precision_ = 0.0;
  mutable Engine engines_;
  mutable array_type shblk_norm_D_;

  void init() {
    ntiles_per_uc_ = basis0_->nclusters();
    assert(ntiles_per_uc_ == aux_basis_->nclusters());

    trange1_aux_ = aux_basis_->create_trange1();
    auto &world = this->get_world();

    using ::mpqc::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::detail::compute_shell_offset;
    using ::mpqc::lcao::gaussian::make_engine_pool;
    using ::mpqc::lcao::gaussian::detail::make_screener;

    // make compound basis set for ν_R in (X_Rj | μ_0 ν_R)
    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    basisR_ = shift_basis_origin(*basis0_, zero_shift_base, R_max_, dcell_);

    // make compound basis set for X in (X_Rj | μ_0 ν_R)
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      aux_basis_RJ_.emplace_back(shift_basis_origin(*aux_basis_, vec_RJ));
    }

    const auto basis0 = *basis0_;
    const auto basisR = *basisR_;
    const auto basis_aux = *aux_basis_;
    // initialize screener for (X_Rj | μ_0 ν_R)
    if (screen_ == "schwarz") {
      const auto basisX =
          *(shift_basis_origin(basis_aux, zero_shift_base, RJ_max_, dcell_));
      auto screen_engine =
          make_engine_pool(libint2::Operator::coulomb,
                           utility::make_array_of_refs(basisX, basis0, basisR),
                           libint2::BraKet::xx_xx);
      auto bases =
          ::mpqc::lcao::gaussian::BasisVector{{basisX, basis0, basisR}};
      p_screener_ = make_screener(world, screen_engine, bases, screen_,
                                  screen_threshold_);
    } else {
      throw InputError("Wrong screening method", __FILE__, __LINE__, "screen");
    }

    // compute the index of the first shell in each cluster
    basis0_shell_offset_map_ = compute_shell_offset(basis0);
    basisR_shell_offset_map_ = compute_shell_offset(basisR);

    // create a TiledRange for three-center ERIs
    trange_eri3_ = ::mpqc::lcao::gaussian::detail::create_trange(
        BasisVector{{basis_aux, basis0, basisR}});
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

  void compute_task_Xmn_mn(Tile D01, Tile norm_D01,
                           std::array<size_t, 2> lattice_ord_idx,
                           std::array<size_t, 3> tile_idx) {
    const auto tile_aux = tile_idx[0];
    const auto tile0 = tile_idx[1];
    const auto tile1 = tile_idx[2];

    const auto R1_ord = lattice_ord_idx[0];
    const auto RJ_ord = lattice_ord_idx[1];

    // translate tile index by unit cell index
    const auto tile1_R1 = tile1 + R1_ord * ntiles_per_uc_;

    // get reference to basis sets
    const auto &basis_aux = aux_basis_RJ_[RJ_ord];
    const auto &basis0 = basis0_;
    const auto &basis1 = basisR_;

    // shell clusters for this tile
    const auto &cluster_aux = basis_aux->cluster_shells()[tile_aux];
    const auto &cluster0 = basis0->cluster_shells()[tile0];
    const auto &cluster1 = basis1->cluster_shells()[tile1_R1];

    // number of shells in each cluster
    const auto nshells_aux = cluster_aux.size();
    const auto nshells0 = cluster0.size();
    const auto nshells1 = cluster1.size();

    // 1-d tile ranges
    const auto &rng_aux = trange_eri3_.dim(0).tile(tile_aux);
    const auto &rng0 = trange_eri3_.dim(1).tile(tile0);
    const auto &rng1 = trange_eri3_.dim(2).tile(tile1_R1);
    const auto rng1_size = rng1.second - rng1.first;

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng_aux});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);

    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();
    const auto *D01_ptr = D01.data();
    const auto *norm_D01_ptr = norm_D01.data();
    assert(D01_ptr != nullptr);
    assert(norm_D01_ptr != nullptr);

    const auto nbf_aux_per_uc = basis_aux->nfunctions();

    using ::mpqc::lcao::gaussian::detail::compute_func_offset_list;

    // compute contributions to all result matrices
    {
      // index of first shell in this cluster
      const auto sh0_offset = basis0_shell_offset_map_[tile0];
      const auto sh1_offset = basisR_shell_offset_map_[tile1_R1];

      // index of last shell in this cluster
      const auto sh0_max = sh0_offset + nshells0;
      const auto sh1_max = sh1_offset + nshells1;

      // determine if this task contains significant shell-pairs
      auto is_significant = false;
      {
        for (auto sh0_in_basis = sh0_offset; sh0_in_basis != sh0_max;
             ++sh0_in_basis) {
          for (const auto sh1_in_basis : sig_shellpair_list_[sh0_in_basis]) {
            if (sh1_in_basis >= sh1_offset && sh1_in_basis < sh1_max) {
              is_significant = true;
              break;
            }
          }
          if (is_significant) break;
        }
      }

      if (is_significant) {
        auto &screen = *(p_screener_);
        auto engine = engines_->local();
        const auto engine_precision = target_precision_;
        engine.set_precision(engine_precision);
        const auto &computed_shell_sets = engine.results();

        // compute offset list of cluster1 (ket1)
        auto offset_list_c1 = compute_func_offset_list(cluster1, rng1.first);
        // this is the index of the first basis functions for each shell *in
        // this shell cluster*
        auto cf0_offset = 0;
        // this is the index of the first basis functions for each shell *in the
        // basis set*
        auto bf0_offset = rng0.first;

        size_t cf1_offset, bf1_offset;
        // loop over all shell sets
        for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
          const auto &shell0 = cluster0[sh0];
          const auto nf0 = shell0.size();

          const auto sh0_in_basis = sh0 + sh0_offset;
          for (const auto &sh1_in_basis : sig_shellpair_list_[sh0_in_basis]) {
            if (sh1_in_basis < sh1_offset || sh1_in_basis >= sh1_max) continue;

            const auto sh1 = sh1_in_basis - sh1_offset;
            std::tie(cf1_offset, bf1_offset) = offset_list_c1[sh1];

            const auto &shell1 = cluster1[sh1];
            const auto nf1 = shell1.size();

            const auto sh01 =
                sh0 * nshells1 + sh1;  // index of {sh0, sh1} in norm_D01
            const double Dnorm01 = norm_D01_ptr[sh01];

            auto cf_aux_offset = 0;
            auto bf_aux_offset = rng_aux.first;
            for (auto sh_aux = 0; sh_aux != nshells_aux; ++sh_aux) {
              const auto &shell_aux = cluster_aux[sh_aux];
              const auto nf_aux = shell_aux.size();

              const auto bf_aux_in_screener =
                  bf_aux_offset + nbf_aux_per_uc * RJ_ord;

              if (screen.skip(bf_aux_in_screener, bf0_offset, bf1_offset,
                              Dnorm01))
                continue;

              // compute shell set
              // TODO call 3-body version of compute2 to avoid extra copies
              engine.compute(shell_aux, shell0, shell1);
              const auto *eri_aux_01 = computed_shell_sets[0];

              if (eri_aux_01 != nullptr) {
                for (auto f_aux = 0, f_aux_01 = 0; f_aux != nf_aux; ++f_aux) {
                  const auto cf_aux = f_aux + cf_aux_offset;
                  for (auto f0 = 0; f0 != nf0; ++f0) {
                    const auto cf0 = f0 + cf0_offset;
                    for (auto f1 = 0; f1 != nf1; ++f1, ++f_aux_01) {
                      const auto cf1 = f1 + cf1_offset;
                      const auto cf01 = cf0 * rng1_size + cf1;

                      const auto value = eri_aux_01[f_aux_01];

                      result_ptr[cf_aux] += D01_ptr[cf01] * value;
                    }
                  }
                }
              }
              cf_aux_offset += nf_aux;
              bf_aux_offset += nf_aux;
            }
          }
          cf0_offset += nf0;
          bf0_offset += nf0;
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
    const auto tile1 = tile_idx[2];

    // get reference to basis sets
    const auto &basis_aux = aux_basis_RJ_[RJ];
    const auto &basis0 = basis0_;
    const auto &basis1 = basisR_;

    // shell clusters for this tile
    const auto &cluster_aux = basis_aux->cluster_shells()[tile_aux];
    const auto &cluster0 = basis0->cluster_shells()[tile0];
    const auto &cluster1 = basis1->cluster_shells()[tile1];

    // number of shells in each cluster
    const auto nshells_aux = cluster_aux.size();
    const auto nshells0 = cluster0.size();
    const auto nshells1 = cluster1.size();

    // 1-d tile ranges
    const auto &rng_aux = trange_eri3_.dim(0).tile(tile_aux);
    const auto &rng0 = trange_eri3_.dim(1).tile(tile0);
    const auto &rng1 = trange_eri3_.dim(2).tile(tile1);
    const auto rng1_size = rng1.second - rng1.first;
    const auto ntiles1 = trange_eri3_.dim(2).tile_extent();

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng0, rng1});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);

    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();
    const auto *X_aux_ptr = X_aux.data();
    assert(X_aux_ptr != nullptr);

    const auto nbf_aux_per_uc = basis_aux->nfunctions();

    using ::mpqc::lcao::gaussian::detail::compute_func_offset_list;

    // TODO: check the sparsity of X. Will screening (X|mu) M_X be useful?
    {
      // index of first shell in this cluster
      const auto sh0_offset = basis0_shell_offset_map_[tile0];
      const auto sh1_offset = basisR_shell_offset_map_[tile1];

      // index of last shell in this cluster
      const auto sh0_max = sh0_offset + nshells0;
      const auto sh1_max = sh1_offset + nshells1;

      // determine if this task contains significant shell-pairs
      auto is_significant = false;
      {
        for (auto sh0_in_basis = sh0_offset; sh0_in_basis != sh0_max;
             ++sh0_in_basis) {
          for (const auto sh1_in_basis : sig_shellpair_list_[sh0_in_basis]) {
            if (sh1_in_basis >= sh1_offset && sh1_in_basis < sh1_max) {
              is_significant = true;
              break;
            }
          }
          if (is_significant) break;
        }
      }

      if (is_significant) {
        auto &screen = *(p_screener_);
        auto engine = engines_->local();
        const auto engine_precision = target_precision_;
        engine.set_precision(engine_precision);
        const auto &computed_shell_sets = engine.results();

        // compute offset list of clusterR (ket1)
        auto offset_list_c1 = compute_func_offset_list(cluster1, rng1.first);

        // this is the index of the first basis functions for each shell *in
        // this shell cluster*
        auto cf0_offset = 0;
        // this is the index of the first basis functions for each shell *in the
        // basis set*
        auto bf0_offset = rng0.first;

        size_t cf1_offset, bf1_offset;
        // loop over all shell sets
        for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
          const auto &shell0 = cluster0[sh0];
          const auto nf0 = shell0.size();

          const auto sh0_in_basis = sh0 + sh0_offset;
          for (const auto &sh1_in_basis : sig_shellpair_list_[sh0_in_basis]) {
            if (sh1_in_basis < sh1_offset || sh1_in_basis >= sh1_max) continue;

            const auto sh1 = sh1_in_basis - sh1_offset;
            std::tie(cf1_offset, bf1_offset) = offset_list_c1[sh1];

            const auto &shell1 = cluster1[sh1];
            const auto nf1 = shell1.size();

            auto cf_aux_offset = 0;
            auto bf_aux_offset = rng_aux.first;
            for (auto sh_aux = 0; sh_aux != nshells_aux; ++sh_aux) {
              const auto &shell_aux = cluster_aux[sh_aux];
              const auto nf_aux = shell_aux.size();

              const auto bf_aux_in_screener =
                  bf_aux_offset + nbf_aux_per_uc * RJ;

              if (screen.skip(bf_aux_in_screener, bf0_offset, bf1_offset))
                continue;

              // compute shell set
              // TODO call 3-body version of compute2 to avoid extra copies
              engine.compute(shell_aux, shell0, shell1);
              const auto *eri_aux_01 = computed_shell_sets[0];

              if (eri_aux_01 != nullptr) {
                for (auto f_aux = 0, f_aux_01 = 0; f_aux != nf_aux; ++f_aux) {
                  const auto cf_aux = f_aux + cf_aux_offset;
                  for (auto f0 = 0; f0 != nf0; ++f0) {
                    const auto cf0 = f0 + cf0_offset;
                    for (auto f1 = 0; f1 != nf1; ++f1, ++f_aux_01) {
                      const auto cf1 = f1 + cf1_offset;
                      const auto cf01 = cf0 * rng1_size + cf1;

                      const auto value = eri_aux_01[f_aux_01];

                      result_ptr[cf01] += X_aux_ptr[cf_aux] * value;
                    }
                  }
                }
              }
              cf_aux_offset += nf_aux;
              bf_aux_offset += nf_aux;
            }
          }
          cf0_offset += nf0;
          bf0_offset += nf0;
        }
      }
    }

    // accumulate the local contributions
    {
      auto tile01 = tile0 * ntiles1 + tile1;
      PeriodicThreeCenterContractionBuilder_::accumulate_local_task(result_tile,
                                                                    tile01);
    }
  }

};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_THREE_CENTER_CONTRACTION_BUILDER_H_
