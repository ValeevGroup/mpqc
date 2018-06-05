//
// Created by Drew Lewis on 05/24/2017
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_FOUR_CENTER_EXACT_K_DIAGONAL_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_FOUR_CENTER_EXACT_K_DIAGONAL_BUILDER_H_

#include <cassert>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/basis/util.h"
#include "mpqc/chemistry/qc/lcao/factory/factory_utility.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace scf {

/// ExactKDiagonalBuilder is an integral-direct implementation of FockBuilder
/// in a Gaussian AO basis that uses 4-center integrals to compute the exact K
/// diagonal and optimally takes / advantage of the permutational symmetry and
/// shell-level screening.
template <typename Tile, typename Policy>
class ExactKDiagonalBuilder
    : public FockBuilder<Tile, Policy>,
      public madness::WorldObject<ExactKDiagonalBuilder<Tile, Policy>> {
 public:
  using array_type = typename FockBuilder<Tile, Policy>::array_type;
  using const_data_ptr = typename Tile::allocator_type::const_pointer;

  using WorldObject_ =
      madness::WorldObject<ExactKDiagonalBuilder<Tile, Policy>>;
  using ExactKDiagonalBuilder_ = ExactKDiagonalBuilder<Tile, Policy>;

  using Basis = ::mpqc::lcao::gaussian::Basis;
  using Shell = ::mpqc::lcao::gaussian::Shell;
  using ShellVec = ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;

  ExactKDiagonalBuilder(madness::World& world,
                        std::shared_ptr<const Basis> bra_basis,
                        std::shared_ptr<const Basis> ket_basis,
                        std::shared_ptr<const Basis> density_basis,
                        std::string screen = "schwarz",
                        double screen_threshold = 1.0e-12)
      : WorldObject_(world),
        bra_basis_(std::move(bra_basis)),
        ket_basis_(std::move(ket_basis)),
        density_basis_(std::move(density_basis)),
        screen_(screen),
        screen_threshold_(screen_threshold) {
    // same basis on each center only
    assert(bra_basis_ == ket_basis_ && bra_basis_ == density_basis_ &&
           "not yet implemented");
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();
  }

  virtual ~ExactKDiagonalBuilder() {}

  array_type operator()(array_type const& D, array_type const&,
                        double target_precision) override {
    // validate preconditions
    auto ntiles_D = density_basis_->nclusters();
    {
      auto trange_D = D.trange();
      auto elements_range_D = trange_D.elements_range();
      auto tiles_range_D = trange_D.tiles_range();
      auto nbf_D = density_basis_->nfunctions();
      assert(elements_range_D.extent(0) == nbf_D &&
             elements_range_D.extent(1) == nbf_D);
      assert(tiles_range_D.extent(0) == ntiles_D &&
             tiles_range_D.extent(1) == ntiles_D);
    }

    if (bra_basis_ == density_basis_ && ket_basis_ == density_basis_)
      return compute_K_aaaa(D, target_precision);
    assert(false && "feature not implemented");
    return array_type{};
  }

  array_type compute_K_aaaa(array_type const& D, double target_precision) {
    dist_pmap_D_ = D.pmap();

    // Copy D and make it replicated.
    array_type D_repl;
    D_repl("i,j") = D("i,j");
    D_repl.make_replicated();

    // prepare input data
    auto& compute_world = this->get_world();
    const auto me = compute_world.rank();
    const auto nproc = compute_world.nproc();
    target_precision_ = target_precision;

    auto ntiles = bra_basis_->nclusters();
    trange_D_ = D_repl.trange();
    pmap_D_ = D_repl.pmap();

    // make the engine pool
    auto oper_type = libint2::Operator::coulomb;
    const auto& basis = *bra_basis_;  // all basis sets are the same
    engines_ = ::mpqc::lcao::gaussian::make_engine_pool(
        oper_type, utility::make_array_of_refs(basis, basis, basis, basis),
        libint2::BraKet::xx_xx);

    // make screener
    auto bases =
        ::mpqc::lcao::gaussian::BasisVector{{basis, basis, basis, basis}};
    p_screener_ = ::mpqc::lcao::gaussian::detail::make_screener(
        compute_world, engines_, bases, screen_, screen_threshold_);

    num_ints_computed_ = 0;

    // make shell block norm of D
    using ::mpqc::lcao::gaussian::detail::compute_shellblock_norm;
    auto shblk_norm_D = compute_shellblock_norm(basis, basis, D);
    shblk_norm_D.make_replicated();  // make sure it is replicated

    // Define this so I don't have to keep removing them.
    bool compute_K_ = true;

    auto empty = TA::Future<Tile>(Tile());
    // todo screen loop with schwarz
    for (auto tile0 = 0ul, tile0123 = 0ul; tile0 != ntiles; ++tile0) {
      for (auto tile1 = 0ul; tile1 <= tile0; ++tile1) {
        for (auto tile2 = 0ul; tile2 <= tile0; ++tile2) {
          // if tile0==tile2 there will be shell blocks such that shell0 >
          // shell2, hence need shell3<=shell2 -> tile3<=tile2
          for (auto tile3 = 0ul; tile3 <= tile2; ++tile3, ++tile0123) {
            bool contains_K_diag = tile0 == tile2 || tile0 == tile3 ||
                                   tile1 == tile2 || tile1 == tile3;
            if (false == contains_K_diag) {
              continue;
            }

            // TODO screen D blocks using schwarz estimate for this Coulomb
            // operator tile
            auto D02 = (!compute_K_ || D_repl.is_zero({tile0, tile2}))
                           ? empty
                           : D_repl.find({tile0, tile2});
            auto D03 = (!compute_K_ || D_repl.is_zero({tile0, tile3}))
                           ? empty
                           : D_repl.find({tile0, tile3});
            auto D12 = (!compute_K_ || D_repl.is_zero({tile1, tile2}))
                           ? empty
                           : D_repl.find({tile1, tile2});
            auto D13 = (!compute_K_ || D_repl.is_zero({tile1, tile3}))
                           ? empty
                           : D_repl.find({tile1, tile3});

            // shell block norms of D
            auto norm_D02 =
                (!compute_K_ || shblk_norm_D.is_zero({tile0, tile2}))
                    ? empty
                    : shblk_norm_D.find({tile0, tile2});
            auto norm_D03 =
                (!compute_K_ || shblk_norm_D.is_zero({tile0, tile3}))
                    ? empty
                    : shblk_norm_D.find({tile0, tile3});
            auto norm_D12 =
                (!compute_K_ || shblk_norm_D.is_zero({tile1, tile2}))
                    ? empty
                    : shblk_norm_D.find({tile1, tile2});
            auto norm_D13 =
                (!compute_K_ || shblk_norm_D.is_zero({tile1, tile3}))
                    ? empty
                    : shblk_norm_D.find({tile1, tile3});

            // clang-format off
            // Using lambda as a task argument fails
            // because madness cannot find type trait (constness) of
            // lambda.
            // To be fixed ...
//            auto task_func = [&, this]() {
//              this->compute_task(
//                  D01, D23, D02, D03, D12, D13,
//                  std::array<size_t, 4>{{tile0, tile1, tile2, tile3}},
//                  std::array<Tile, 6>{{norm_D01, norm_D23, norm_D02, norm_D03,
//                                       norm_D12, norm_D13}});
//            };
//            if (pmap->is_local(tile0123)) WorldObject_::task(me, task_func);
            // clang-format on

            if (tile0123 % nproc == me)
              WorldObject_::task(
                  me, &ExactKDiagonalBuilder::compute_task, D02, D03, D12, D13,
                  std::array<size_t, 4>{{tile0, tile1, tile2, tile3}},
                  std::array<Tile, 4>{
                      {norm_D02, norm_D03, norm_D12, norm_D13}});
          }
        }
      }
    }

    // fence ensures everyone is done
    compute_world.gop.fence();

    // cleanup
    engines_.reset();

    ExEnv::out0() << "\nIntegrals per node:" << std::endl;
    for (auto i = 0; i < compute_world.nproc(); ++i) {
      if (me == i) {
        ExEnv::outn() << indent << "Integrals on node(" << i
                      << "): " << num_ints_computed_ << std::endl;
      }
      compute_world.gop.fence();
    }
    ExEnv::out0() << std::endl;

    if (pmap_D_->is_replicated() && compute_world.size() > 1) {
      // Each process has its own copy of G which is treated differently.
      // Reduce all G's to a dist array
      for (const auto& local_tile : local_fock_tiles_) {
        const auto ij = local_tile.first;
        const auto proc01 = dist_pmap_D_->owner(ij);
        WorldObject_::task(proc01, &ExactKDiagonalBuilder_::accumulate_array,
                           local_tile.second, ij);
      }
      local_fock_tiles_.clear();

      compute_world.gop.fence();

      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> global_tile_norms;
        for (const auto& global_tile : global_fock_tiles_) {
          const auto ij = global_tile.first;
          const auto i = ij / ntiles;
          const auto j = ij % ntiles;
          const auto ij_norm = global_tile.second.norm();
          global_tile_norms.push_back(
              std::make_pair(std::array<size_t, 2>{{i, j}}, ij_norm));
        }
        shape = decltype(shape)(compute_world, global_tile_norms, trange_D_);
      }

      array_type G_dist(compute_world, trange_D_, shape, dist_pmap_D_);
      for (const auto& global_tile : global_fock_tiles_) {
        if (!G_dist.shape().is_zero(global_tile.first))
          G_dist.set(global_tile.first, global_tile.second);
      }
      G_dist.fill_local(0.0, true);
      global_fock_tiles_.clear();

      // symmetrize to account for permutation symmetry use
      G_dist("i,j") = 0.5 * (G_dist("i,j") + G_dist("j,i"));
      return G_dist;

    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> local_tile_norms;
        for (const auto& local_tile_iter : local_fock_tiles_) {
          const auto ij = local_tile_iter.first;
          const auto i = ij / ntiles;
          const auto j = ij % ntiles;
          const auto ij_norm = local_tile_iter.second.norm();
          local_tile_norms.push_back(
              std::make_pair(std::array<size_t, 2>{{i, j}}, ij_norm));
        }
        shape = decltype(shape)(compute_world, local_tile_norms, trange_D_);
      }

      array_type G(compute_world, trange_D_, shape, pmap_D_);

      // copy results of local reduction tasks into the local copy of G
      for (const auto& local_tile : local_fock_tiles_) {
        // if this tile was not truncated away
        if (!G.shape().is_zero(local_tile.first))
          G.set(local_tile.first, local_tile.second);
      }
      // set the remaining local tiles to 0 (this should only be needed for
      // dense policy)
      G.fill_local(0.0, true);

      local_fock_tiles_.clear();

      // symmetrize to account for permutation symmetry use
      G("i,j") = 0.5 * (G("i,j") + G("j,i"));

      return G;
    }
  }

  void register_fock(const array_type& fock,
                     FormulaRegistry<array_type>& registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  inline void print_iter(std::string const& leader) override {}

 private:
  // set by ctor
  std::shared_ptr<const Basis> bra_basis_;
  std::shared_ptr<const Basis> ket_basis_;
  std::shared_ptr<const Basis> density_basis_;
  const std::string screen_;
  const double screen_threshold_;

  // mutated by compute_ functions
  std::shared_ptr<lcao::Screener> p_screener_;
  madness::ConcurrentHashMap<std::size_t, Tile> local_fock_tiles_;
  madness::ConcurrentHashMap<std::size_t, Tile> global_fock_tiles_;
  TA::TiledRange trange_D_;
  std::shared_ptr<TA::Pmap> pmap_D_;
  std::shared_ptr<TA::Pmap> dist_pmap_D_;
  double target_precision_ = 0.0;
  ::mpqc::lcao::gaussian::ShrPool<libint2::Engine> engines_;
  std::atomic<size_t> num_ints_computed_{0};

  void accumulate_array(Tile arg_tile, long tile01) {
    typename decltype(global_fock_tiles_)::accessor acc;
    if (!global_fock_tiles_.insert(
            acc, std::make_pair(tile01, arg_tile))) {  // CRITICAL SECTION
      const auto size = arg_tile.range().volume();
      TA::math::inplace_vector_op_serial(
          [](TA::detail::numeric_t<Tile>& l,
             const TA::detail::numeric_t<Tile> r) { l += r; },
          size, acc->second.data(), arg_tile.data());
    }
    acc.release();  // END OF CRITICAL SECTION
  }

  void accumulate_task(Tile fock_matrix_tile, long tile0, long tile1) {
    const auto ntiles = trange_D_.dim(0).tile_extent();
    const auto tile01 = tile0 * ntiles + tile1;
    assert(pmap_D_->is_local(tile01));
    // if reducer does not exist, create entry and store F, else accumulate F to
    // the existing contents
    typename decltype(local_fock_tiles_)::accessor acc;
    // try inserting, otherwise, accumulate
    if (!local_fock_tiles_.insert(
            acc,
            std::make_pair(tile01, fock_matrix_tile))) {  // CRITICAL SECTION
      // NB can't do acc->second += fock_matrix_tile to avoid spawning TBB
      // tasks from critical section
      const auto size = fock_matrix_tile.range().volume();
      TA::math::inplace_vector_op_serial(
          [](TA::detail::numeric_t<Tile>& l,
             const TA::detail::numeric_t<Tile> r) { l += r; },
          size, acc->second.data(), fock_matrix_tile.data());
    }
    acc.release();  // END OF CRITICAL SECTION
  }

  void compute_task(Tile D02, Tile D03, Tile D12, Tile D13,
                    std::array<size_t, 4> tile_idx,
                    std::array<Tile, 4> norm_D) {
    const auto tile0 = tile_idx[0];
    const auto tile1 = tile_idx[1];
    const auto tile2 = tile_idx[2];
    const auto tile3 = tile_idx[3];

    // 1-d tile ranges
    const auto& trange1 = trange_D_.dim(0);
    const auto& rng0 = trange1.tile(tile0);
    const auto& rng1 = trange1.tile(tile1);
    const auto& rng2 = trange1.tile(tile2);
    const auto& rng3 = trange1.tile(tile3);
    const auto rng2_size = rng2.second - rng2.first;
    const auto rng3_size = rng3.second - rng3.first;

    // 2-d tile ranges describing the Fock contribution blocks produced by this
    auto rng02 = TA::Range({rng0, rng2});
    auto rng03 = TA::Range({rng0, rng3});
    auto rng12 = TA::Range({rng1, rng2});
    auto rng13 = TA::Range({rng1, rng3});

    // initialize contribution to the Fock matrices
    auto F02 = Tile(std::move(rng02), 0.0);
    auto F03 = Tile(std::move(rng03), 0.0);
    auto F12 = Tile(std::move(rng12), 0.0);
    auto F13 = Tile(std::move(rng13), 0.0);

    // grab ptrs to tile data to make addressing more efficient
    auto* F02_ptr = F02.data();
    auto* F03_ptr = F03.data();
    auto* F12_ptr = F12.data();
    auto* F13_ptr = F13.data();
    const auto* D02_ptr = D02.data();
    const auto* D03_ptr = D03.data();
    const auto* D12_ptr = D12.data();
    const auto* D13_ptr = D13.data();
    const auto* norm_D02_ptr = norm_D[0].data();
    const auto* norm_D03_ptr = norm_D[1].data();
    const auto* norm_D12_ptr = norm_D[2].data();
    const auto* norm_D13_ptr = norm_D[3].data();

    // compute contributions to all Fock matrices
    {
      auto& screen = *p_screener_;
      auto engine = engines_->local();
      const auto engine_precision = target_precision_;
      engine.set_precision(engine_precision);
      const auto& computed_shell_sets = engine.results();

      const auto& basis = bra_basis_;
      // shell clusters for this tile
      const auto& cluster0 = basis->cluster_shells()[tile0];
      const auto& cluster1 = basis->cluster_shells()[tile1];
      const auto& cluster2 = basis->cluster_shells()[tile2];
      const auto& cluster3 = basis->cluster_shells()[tile3];

      // make unique shell pair list
      const auto same_c0c1 = tile0 == tile1;
      const auto same_c2c3 = tile2 == tile3;
      const auto same_c0c2 = tile0 == tile2;
      shellpair_list_t bra_shellpair_list, ket_shellpair_list;
      if (same_c0c1 && same_c2c3 && same_c0c2) {
        bra_shellpair_list = compute_shellpair_list(cluster0);
        ket_shellpair_list = bra_shellpair_list;
      } else if (same_c0c1 && same_c2c3 && !same_c0c2) {
        bra_shellpair_list = compute_shellpair_list(cluster0);
        ket_shellpair_list = compute_shellpair_list(cluster2);
      } else if (same_c0c1 && !same_c2c3) {
        bra_shellpair_list = compute_shellpair_list(cluster0);
        ket_shellpair_list = compute_shellpair_list(cluster2, cluster3);
      } else if (!same_c0c1 && same_c2c3) {
        bra_shellpair_list = compute_shellpair_list(cluster0, cluster1);
        ket_shellpair_list = compute_shellpair_list(cluster2);
      } else {
        bra_shellpair_list = compute_shellpair_list(cluster0, cluster1);
        ket_shellpair_list = compute_shellpair_list(cluster2, cluster3);
      }

      // number of shells in each cluster
      const auto nshells0 = cluster0.size();
      const auto nshells2 = cluster2.size();
      const auto nshells3 = cluster3.size();

      // compute offset list of cluster1 and cluster3
      using ::mpqc::lcao::gaussian::detail::compute_func_offset_list;
      auto offset_list_c1 = compute_func_offset_list(cluster1, rng1.first);
      auto offset_list_c3 = compute_func_offset_list(cluster3, rng3.first);

      // this is the index of the first basis functions for each shell *in this
      // shell cluster*
      auto cf0_offset = 0;
      // this is the index of the first basis functions for each shell *in the
      // basis set*
      auto bf0_offset = rng0.first;

      size_t cf1_offset, bf1_offset, cf3_offset, bf3_offset;

      // loop over unique shell sets
      // N.B. skip nonunique shell sets that did not get eliminated by unique
      // cluster set iteration
      for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
        const auto& shell0 = cluster0[sh0];
        const auto nf0 = shell0.size();

        for (const auto& sh1 : bra_shellpair_list[sh0]) {
          std::tie(cf1_offset, bf1_offset) = offset_list_c1[sh1];
          // skip if shell set is nonunique
          if (bf0_offset < bf1_offset)
            break;  // assuming basis functions increase monotonically in the
                    // basis

          const auto& shell1 = cluster1[sh1];
          const auto nf1 = shell1.size();

          const auto multiplicity01 = bf0_offset == bf1_offset ? 1.0 : 2.0;

          auto cf2_offset = 0;
          auto bf2_offset = rng2.first;

          for (auto sh2 = 0; sh2 != nshells2; ++sh2) {
            // skip if shell set is nonunique
            if (bf0_offset < bf2_offset) break;

            const auto& shell2 = cluster2[sh2];
            const auto nf2 = shell2.size();

            const auto sh02 =
                sh0 * nshells2 + sh2;  // index of {sh0, sh2} in norm_D02
            const auto sh12 =
                sh1 * nshells2 + sh2;  // index of {sh1, sh2} in norm_D12
            const auto Dnorm12 = (tile0 == tile3) ? norm_D12_ptr[sh12] : 0.0;
            const auto Dnorm02 = (tile1 == tile3) ? norm_D02_ptr[sh02] : 0.0;
            const auto Dnorm012 = std::max({Dnorm02, Dnorm12});

            for (const auto& sh3 : ket_shellpair_list[sh2]) {
              std::tie(cf3_offset, bf3_offset) = offset_list_c3[sh3];
              // skip if shell set is nonunique
              if (bf2_offset < bf3_offset ||
                  (bf0_offset == bf2_offset && bf1_offset < bf3_offset))
                break;

              const auto& shell3 = cluster3[sh3];
              const auto nf3 = shell3.size();

              const auto sh03 =
                  sh0 * nshells3 + sh3;  // index of {sh0, sh3} in norm_D03
              const auto sh13 =
                  sh1 * nshells3 + sh3;  // index of {sh1, sh3} in norm_D13
              const auto Dnorm03 = (tile1 == tile2) ? norm_D03_ptr[sh03] : 0.0;
              const auto Dnorm13 = (tile0 == tile2) ? norm_D13_ptr[sh13] : 0.0;
              const auto Dnorm0123 = std::max({Dnorm03, Dnorm13, Dnorm012});

              if (screen.skip(bf0_offset, bf1_offset, bf2_offset, bf3_offset,
                              Dnorm0123))
                continue;

              num_ints_computed_ += nf0 * nf1 * nf2 * nf3;

              const auto multiplicity23 = bf2_offset == bf3_offset ? 1.0 : 2.0;
              const auto multiplicity0213 =
                  (bf0_offset == bf2_offset && bf1_offset == bf3_offset) ? 1.0
                                                                         : 2.0;
              const auto multiplicity =
                  multiplicity01 * multiplicity23 * multiplicity0213;

              // compute shell set
              engine.compute2<libint2::Operator::coulomb,
                              libint2::BraKet::xx_xx, 0>(shell0, shell1, shell2,
                                                         shell3);
              const auto* eri_0123 = computed_shell_sets[0];

              if (eri_0123 !=
                  nullptr) {  // if the shell set is not screened out

                for (auto f0 = 0, f0123 = 0; f0 != nf0; ++f0) {
                  const auto cf0 = f0 + cf0_offset;  // basis function index in
                                                     // the tile (i.e. shell
                                                     // cluster)
                  for (auto f1 = 0; f1 != nf1; ++f1) {
                    const auto cf1 = f1 + cf1_offset;
                    for (auto f2 = 0; f2 != nf2; ++f2) {
                      const auto cf2 = f2 + cf2_offset;
                      const auto cf02 =
                          cf0 * rng2_size +
                          cf2;  // index of {cf0,cf2} in D02 or F02
                      const auto cf12 =
                          cf1 * rng2_size +
                          cf2;  // index of {cf1,cf2} in D12 or F12
                      for (auto f3 = 0; f3 != nf3; ++f3, ++f0123) {
                        const auto cf3 = f3 + cf3_offset;
                        const auto cf03 =
                            cf0 * rng3_size +
                            cf3;  // index of {cf0,cf3} in D03 or F03
                        const auto cf13 =
                            cf1 * rng3_size +
                            cf3;  // index of {cf1,cf3} in D13 or F13

                        const auto value = eri_0123[f0123];

                        const auto value_scaled_by_multiplicity =
                            value * multiplicity;

                        F02_ptr[cf02] -= (D13_ptr != nullptr)
                                             ? 0.25 * D13_ptr[cf13] *
                                                   value_scaled_by_multiplicity
                                             : 0.0;
                        F13_ptr[cf13] -= (D02_ptr != nullptr)
                                             ? 0.25 * D02_ptr[cf02] *
                                                   value_scaled_by_multiplicity
                                             : 0.0;
                        F03_ptr[cf03] -= (D12_ptr != nullptr)
                                             ? 0.25 * D12_ptr[cf12] *
                                                   value_scaled_by_multiplicity
                                             : 0.0;
                        F12_ptr[cf12] -= (D03_ptr != nullptr)
                                             ? 0.25 * D03_ptr[cf03] *
                                                   value_scaled_by_multiplicity
                                             : 0.0;
                      }
                    }
                  }
                }
              }
            }

            cf2_offset += nf2;
            bf2_offset += nf2;
          }
        }

        cf0_offset += nf0;
        bf0_offset += nf0;
      }
    }

    // accumulate the contributions by submitting tasks to the owners of their
    // tiles
    ExactKDiagonalBuilder_::accumulate_task(F02, tile_idx[0], tile_idx[2]);
    ExactKDiagonalBuilder_::accumulate_task(F03, tile_idx[0], tile_idx[3]);
    ExactKDiagonalBuilder_::accumulate_task(F12, tile_idx[1], tile_idx[2]);
    ExactKDiagonalBuilder_::accumulate_task(F13, tile_idx[1], tile_idx[3]);
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
      const ShellVec& shv1,
      const ShellVec& _shv2 = std::vector<Shell>({Shell()}),
      double threshold = 1e-12) {
    const ShellVec& shv2 =
        ((_shv2.size() == 1 && _shv2[0] == Shell()) ? shv1 : _shv2);
    const auto nsh1 = shv1.size();
    const auto nsh2 = shv2.size();
    const auto shv1_equiv_shv2 = (&shv1 == &shv2);

    // determine max # of primitives in a shell cluster
    auto max_nprim = [](const ShellVec& shv) {
      size_t n = 0;
      for (auto shell : shv) n = std::max(shell.nprim(), n);
      return n;
    };
    const auto max_nprim_1 = max_nprim(shv1);
    const auto max_nprim_2 = max_nprim(shv2);

    // determine max angular momentum of a shell cluster
    auto max_l = [](const ShellVec& shv) {
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
    const auto& buf = engine.results();
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
      auto& list = result[s1];
      std::sort(list.begin(), list.end());
    }

    return result;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_FOUR_CENTER_EXACT_K_DIAGONAL_BUILDER_H_
