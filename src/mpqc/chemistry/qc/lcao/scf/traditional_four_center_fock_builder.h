
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_

#include <cassert>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/factory/factory_utility.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace scf {

/// ReferenceFourCenterFockBuilder is a reference implementation
/// that uses 4-center (stored or direct) integrals.

/// @warning This is a very inefficient builder when used with direct integrals:
/// the number of integrals it evaluates is roughly 6 to 16 times greater than
/// optimal (this depends on whether the integral arrays account for any
/// permutational symmetry). It is only useful for reference computation
template <typename Tile, typename Policy, typename Integral>
class ReferenceFourCenterFockBuilder : public FockBuilder<Tile, Policy> {
 public:
  using array_type = typename FockBuilder<Tile, Policy>::array_type;

  ReferenceFourCenterFockBuilder(Integral const &eri4_J, Integral const &eri4_K)
      : eri4_J_(eri4_J), eri4_K_(eri4_K) {}

  array_type operator()(array_type const &D, array_type const &,
                        double) override {
    const auto make_J = eri4_J_.is_initialized();
    const auto make_K = eri4_K_.is_initialized();
    assert(make_J || make_K);

		// Make J
    array_type J;
    if (make_J) {
      J("mu, nu") = eri4_J_("mu, nu, rho, sig") * D("rho, sig");
      // symmetrize to account for petite list
      J("mu, nu") = 0.5 * (J("mu, nu") + J("nu, mu"));
    }

    // Make K
    array_type K;
    if (make_K) {
      K("mu, nu") = eri4_K_("mu, rho, nu, sig") * D("rho, sig");
      // symmetrize to account for petite list
      K("mu, nu") = 0.5 * (K("mu, nu") + K("nu, mu"));
    }

    // Make and return G
    array_type G;
    if (make_J && make_K)
      G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");
    else if (make_J)
      G("mu, nu") = 2 * J("mu, nu");
    else if (make_K)
      G("mu, nu") = -K("mu, nu");

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  inline void print_iter(std::string const &leader) override {}

 private:
  Integral eri4_J_;
  Integral eri4_K_;
};

/// FourCenterFockBuilder is an integral-direct implementation of FockBuilder
/// in a Gaussian AO basis that uses 4-center integrals and optimally takes
/// advantage of the permutational symmetry and shell-level screening.
template <typename Tile, typename Policy>
class FourCenterFockBuilder
    : public FockBuilder<Tile, Policy>,
      public madness::WorldObject<FourCenterFockBuilder<Tile, Policy>> {
 public:
  using array_type = typename FockBuilder<Tile, Policy>::array_type;

  using WorldObject_ =
      madness::WorldObject<FourCenterFockBuilder<Tile, Policy>>;
  using FourCenterFockBuilder_ = FourCenterFockBuilder<Tile, Policy>;

  using Basis = ::mpqc::lcao::gaussian::Basis;
	using Shell = typename ::mpqc::lcao::gaussian::Shell;
	using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
	using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

  FourCenterFockBuilder(madness::World &world,
                        std::shared_ptr<const Basis> bra_basis,
                        std::shared_ptr<const Basis> ket_basis,
                        std::shared_ptr<const Basis> density_basis,
                        bool compute_J, bool compute_K,
                        std::string screen = "schwarz",
                        double screen_threshold = 1.0e-10)
      : WorldObject_(world),
        compute_J_(compute_J),
        compute_K_(compute_K),
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

  virtual ~FourCenterFockBuilder() {}

  array_type operator()(array_type const &D, array_type const &,
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
      return compute_JK_aaaa(D, target_precision);
    assert(false && "feature not implemented");
  }

  array_type compute_JK_aaaa(array_type const &D,
                             double target_precision) const {
		// prepare input data
    auto &compute_world = this->get_world();
		auto t_jk_0 = mpqc::fenced_now(compute_world);
		const auto me = compute_world.rank();
    target_precision_ = target_precision;

		auto t_prepare_0 = mpqc::fenced_now(compute_world);
		auto nsh = bra_basis_->nshells();
    auto ntiles = bra_basis_->nclusters();
    trange_D_ = D.trange();
    pmap_D_ = D.pmap();
    auto trange1 = trange_D_.dim(0);
    auto trange = TA::TiledRange({trange1, trange1, trange1, trange1});
    const auto ntile_tasks =
        static_cast<uint64_t>(ntiles * (ntiles + 1)) *
        static_cast<uint64_t>((ntiles + 2) * (3 * ntiles + 1)) / 24;
    auto pmap = std::make_shared<const TA::detail::BlockedPmap>(
        compute_world, ntile_tasks);

		auto t_prepare_1 = mpqc::fenced_now(compute_world);
		auto dur_prepare = mpqc::duration_in_s(t_prepare_0, t_prepare_1);

		auto t_eng_0 = mpqc::fenced_now(compute_world);
		// make the engine pool
    auto oper_type = libint2::Operator::coulomb;
    const auto& basis = *bra_basis_;  // all basis sets are the same
    engines_ = ::mpqc::lcao::gaussian::make_engine_pool(
        oper_type,
        utility::make_array_of_refs(basis, basis, basis, basis),
        libint2::BraKet::xx_xx);
		auto t_eng_1 = mpqc::fenced_now(compute_world);
		auto dur_eng = mpqc::duration_in_s(t_eng_0, t_eng_1);

		auto t_screen_0 = mpqc::fenced_now(compute_world);
		// make screener
    auto bases = ::mpqc::lcao::gaussian::BasisVector{{basis, basis, basis, basis}};
    p_screener_ = ::mpqc::lcao::gaussian::detail::make_screener(compute_world, engines_, bases, screen_, screen_threshold_);
		auto t_screen_1 = mpqc::fenced_now(compute_world);
		auto dur_screen = mpqc::duration_in_s(t_screen_0, t_screen_1);

		auto t_shblknormD_0 = mpqc::fenced_now(compute_world);
		// make shell block norm of D
		shblk_norm_D_ = compute_shellblock_norm(basis, D);
		auto t_shblknormD_1 = mpqc::fenced_now(compute_world);
		auto dur_shblknormD = mpqc::duration_in_s(t_shblknormD_0, t_shblknormD_1);

		num_ints_computed_ = 0;
		dur_int_compute_ = 0.0;
		dur_contr_ = 0.0;
		dur_preLoopShell_ = 0.0;
		dur_skip_ = 0.0;
		dur_fock_tile_ = 0.0;
		dur_accumul_ = 0.0;
		dur_preLoop_ = 0.0;
		dur_preLoopInts_s0_ = 0.0;
		dur_preLoopInts_s1_ = 0.0;
		dur_preLoopInts_s2_ = 0.0;
		dur_preLoopInts_s3_ = 0.0;
		dur_s3_ = 0.0;

		auto t_task_0 = mpqc::fenced_now(compute_world);
    // compute shell-level schwarz matrices Q
    //auto Q = make_schwarz_Q_shells(compute_world, bra_basis_, bra_basis_);

    auto empty = TA::Future<Tile>(Tile());
		double dur_findD = 0.0;
    // todo screen loop with schwarz
    for (auto tile0 = 0ul, tile0123 = 0ul; tile0 != ntiles; ++tile0) {
      for (auto tile1 = 0ul; tile1 <= tile0; ++tile1) {
        for (auto tile2 = 0ul; tile2 <= tile0; ++tile2) {
          // if tile0==tile2 there will be shell blocks such that shell0 >
          // shell2, hence need shell3<=shell2 -> tile3<=tile2
          for (auto tile3 = 0ul; tile3 <= tile2; ++tile3, ++tile0123) {
            // TODO screen D blocks using schwarz estimate for this Coulomb
            // operator tile
						auto t_findD_0 = mpqc::fenced_now(compute_world);
						auto D01 =
                (!compute_J_ || D.is_zero({tile0, tile1})) ? empty : D.find({tile0, tile1});
            auto D23 =
                (!compute_J_ || D.is_zero({tile2, tile3})) ? empty : D.find({tile2, tile3});
            auto D02 =
                (!compute_K_ || D.is_zero({tile0, tile2})) ? empty : D.find({tile0, tile2});
            auto D03 =
                (!compute_K_ || D.is_zero({tile0, tile3})) ? empty : D.find({tile0, tile3});
            auto D12 =
                (!compute_K_ || D.is_zero({tile1, tile2})) ? empty : D.find({tile1, tile2});
            auto D13 =
                (!compute_K_ || D.is_zero({tile1, tile3})) ? empty : D.find({tile1, tile3});
						auto t_findD_1 = mpqc::fenced_now(compute_world);
						dur_findD += mpqc::duration_in_s(t_findD_0, t_findD_1);

						if (pmap->is_local(tile0123))
							WorldObject_::task(
									me, &FourCenterFockBuilder_::compute_task, D01, D23, D02, D03, D12, D13,
									std::array<size_t, 4>{{tile0, tile1, tile2, tile3}});

          }
        }
      }
    }

    // fence ensures everyone is done
    compute_world.gop.fence();

    // cleanup
    engines_.reset();
		auto t_task_1 = mpqc::fenced_now(compute_world);
		auto dur_task = mpqc::duration_in_s(t_task_0, t_task_1);

		auto t_shape_0 = mpqc::fenced_now(compute_world);
		typename Policy::shape_type shape;
    // compute the shape, if sparse
    if (!decltype(shape)::is_dense()) {
      // extract local contribution to the shape of G, construct global shape
      std::vector<std::pair<std::array<size_t, 2>, double>> local_tile_norms;
      for (const auto &local_tile_iter : local_fock_tiles_) {
        const auto ij = local_tile_iter.first;
        const auto i = ij / ntiles;
        const auto j = ij % ntiles;
        const auto ij_norm = local_tile_iter.second.norm();
        local_tile_norms.push_back(std::make_pair(std::array<size_t,2>{{i, j}}, ij_norm));
      }
      shape = decltype(shape)(compute_world, local_tile_norms, trange_D_);
    }
		auto t_shape_1 = mpqc::fenced_now(compute_world);
		auto dur_shape = mpqc::duration_in_s(t_shape_0, t_shape_1);

		auto t_g_0 = mpqc::fenced_now(compute_world);
		array_type G(compute_world, trange_D_, shape, pmap_D_);

    // copy results of local reduction tasks into G
    for (const auto &local_tile : local_fock_tiles_) {
      // if this tile was not truncated away
      if (!G.shape().is_zero(local_tile.first))
        G.set(local_tile.first, local_tile.second);
    }
    // set the remaining local tiles to 0 (this should only be needed for dense policy)
    G.fill_local(0.0, true);
    local_fock_tiles_.clear();

    // symmetrize to account for permutation symmetry use
    G("i,j") = 0.5 * (G("i,j") + G("j,i"));
		auto t_g_1 = mpqc::fenced_now(compute_world);
		auto dur_g = mpqc::duration_in_s(t_g_0, t_g_1);

		auto t_jk_1 = mpqc::fenced_now(compute_world);
		auto dur_jk = mpqc::duration_in_s(t_jk_0, t_jk_1);

		ExEnv::out0() << "\n# of integrals = " << num_ints_computed_ << std::endl;

		ExEnv::out0() << "\ntime for compute_JK_aaaa = " << dur_jk << std::endl;
		ExEnv::out0() << "\ttime for prepare input = " << dur_prepare << std::endl;
		ExEnv::out0() << "\ttime for init engine = " << dur_eng << std::endl;
		ExEnv::out0() << "\ttime for init screener = " << dur_screen << std::endl;
		ExEnv::out0() << "\ttime for computing shell block norm of D = " << dur_shblknormD << std::endl;
		ExEnv::out0() << "\ttime for computing tasks = " << dur_task << std::endl;
		ExEnv::out0() << "\t\ttime for finding D tiles = " << dur_findD << std::endl;
		ExEnv::out0() << "\t\ttime for pre Shell loop = " << dur_preLoopShell_ << std::endl;
		ExEnv::out0() << "\t\ttime for 6 fock tiles = " << dur_fock_tile_ << std::endl;

		ExEnv::out0() << "\t\t\ttime for pre loop = " << dur_preLoop_ << std::endl;
		ExEnv::out0() << "\t\t\ttime for pre Ints loop shell0 = " << dur_preLoopInts_s0_ << std::endl;
		ExEnv::out0() << "\t\t\ttime for pre Ints loop shell1 = " << dur_preLoopInts_s1_ << std::endl;
		ExEnv::out0() << "\t\t\ttime for pre Ints loop shell2 = " << dur_preLoopInts_s2_ << std::endl;
		ExEnv::out0() << "\t\t\ttime for pre Ints loop shell3 = " << dur_preLoopInts_s3_ << std::endl;
//		ExEnv::out0() << "\t\t\ttime for loop shell3 = " << dur_s3_ << std::endl;

		ExEnv::out0() << "\t\t\t\ttime for skipping = " << dur_skip_ << std::endl;
		ExEnv::out0() << "\t\t\t\ttime for integrals = " << dur_int_compute_ << std::endl;
		ExEnv::out0() << "\t\t\t\ttime for contraction = " << dur_contr_ << std::endl;

		ExEnv::out0() << "\t\ttime for accumulate = " << dur_accumul_ << std::endl;
		ExEnv::out0() << "\ttime for shape = " << dur_shape << std::endl;
		ExEnv::out0() << "\ttime for g = " << dur_g << std::endl;

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  inline void print_iter(std::string const &leader) override {}

 private:
  // set by ctor
  const bool compute_J_;
  const bool compute_K_;
  std::shared_ptr<const Basis> bra_basis_;
  std::shared_ptr<const Basis> ket_basis_;
  std::shared_ptr<const Basis> density_basis_;
  const std::string screen_;
  const double screen_threshold_;

  // mutated by compute_ functions
  mutable std::shared_ptr<lcao::Screener> p_screener_;
  mutable madness::ConcurrentHashMap<std::size_t, Tile>
      local_fock_tiles_;
  mutable TA::TiledRange trange_D_;
  mutable std::shared_ptr<TA::Pmap> pmap_D_;
  mutable double target_precision_ = 0.0;
  mutable ::mpqc::lcao::gaussian::ShrPool<libint2::Engine> engines_;
	mutable array_type shblk_norm_D_;
	mutable std::atomic<size_t> num_ints_computed_{0};
	mutable double dur_int_compute_;
	mutable double dur_contr_;
	mutable double dur_preLoopShell_;
	mutable double dur_skip_;
	mutable double dur_fock_tile_;
	mutable double dur_accumul_;
	mutable double dur_preLoop_;
	mutable double dur_preLoopInts_s0_;
	mutable double dur_preLoopInts_s1_;
	mutable double dur_preLoopInts_s2_;
	mutable double dur_preLoopInts_s3_;
	mutable double dur_s3_;


  void accumulate_task(Tile fock_matrix_tile, long tile0, long tile1) {
    const auto ntiles = trange_D_.dim(0).tile_extent();
    const auto tile01 = tile0 * ntiles + tile1;
    assert(pmap_D_->is_local(tile01));
    // if reducer does not exist, create entry and store F, else accumulate F to the existing contents
    typename decltype(local_fock_tiles_)::accessor acc;
    // try inserting, otherwise, accumulate
    if (!local_fock_tiles_.insert(acc, std::make_pair(tile01, fock_matrix_tile))) {  // CRITICAL SECTION
      // NB can't do acc->second += fock_matrix_tile to avoid spawning TBB
      // tasks from critical section
      const auto size = fock_matrix_tile.range().volume();
      TA::math::inplace_vector_op_serial(
          [](TA::detail::numeric_t<Tile>& l,
             const TA::detail::numeric_t<Tile> r) { l += r; },
          size, acc->second.data(), fock_matrix_tile.data());
    }
    acc.release();   // END OF CRITICAL SECTION
  }

  void compute_task(Tile D01, Tile D23, Tile D02, Tile D03, Tile D12, Tile D13,
                    std::array<size_t, 4> tile_idx) {
		const auto t_preLoopShell_0 = mpqc::now();

		const auto tile0 = tile_idx[0];
		const auto tile1 = tile_idx[1];
		const auto tile2 = tile_idx[2];
		const auto tile3 = tile_idx[3];

		// 1-d tile ranges
    const auto &trange1 = trange_D_.dim(0);
    const auto ntiles = trange1.tile_extent();
		const auto &rng0 = trange1.tile(tile0);
		const auto &rng1 = trange1.tile(tile1);
		const auto &rng2 = trange1.tile(tile2);
		const auto &rng3 = trange1.tile(tile3);
    const auto rng0_size = rng0.second - rng0.first;
    const auto rng1_size = rng1.second - rng1.first;
    const auto rng2_size = rng2.second - rng2.first;
    const auto rng3_size = rng3.second - rng3.first;

    // 2-d tile ranges describing the Fock contribution blocks produced by this
    auto rng01 = compute_J_ ? TA::Range({rng0, rng1}) : TA::Range();
    auto rng23 = compute_J_ ? TA::Range({rng2, rng3}) : TA::Range();
    auto rng02 = compute_K_ ? TA::Range({rng0, rng2}) : TA::Range();
    auto rng03 = compute_K_ ? TA::Range({rng0, rng3}) : TA::Range();
    auto rng12 = compute_K_ ? TA::Range({rng1, rng2}) : TA::Range();
    auto rng13 = compute_K_ ? TA::Range({rng1, rng3}) : TA::Range();

    // initialize contribution to the Fock matrices
    auto F01 = compute_J_ ? Tile(std::move(rng01), 0.0) : Tile();
    auto F23 = compute_J_ ? Tile(std::move(rng23), 0.0) : Tile();
    auto F02 = compute_K_ ? Tile(std::move(rng02), 0.0) : Tile();
    auto F03 = compute_K_ ? Tile(std::move(rng03), 0.0) : Tile();
    auto F12 = compute_K_ ? Tile(std::move(rng12), 0.0) : Tile();
    auto F13 = compute_K_ ? Tile(std::move(rng13), 0.0) : Tile();

		// find shell block norm of D
		auto norm_D01 =
				(!compute_J_ || shblk_norm_D_.is_zero({tile0, tile1})) ? Tile() : shblk_norm_D_.find({tile0, tile1});
		auto norm_D23 =
				(!compute_J_ || shblk_norm_D_.is_zero({tile2, tile3})) ? Tile() : shblk_norm_D_.find({tile2, tile3});
		auto norm_D02 =
				(!compute_K_ || shblk_norm_D_.is_zero({tile0, tile2})) ? Tile() : shblk_norm_D_.find({tile0, tile2});
		auto norm_D03 =
				(!compute_K_ || shblk_norm_D_.is_zero({tile0, tile3})) ? Tile() : shblk_norm_D_.find({tile0, tile3});
		auto norm_D12 =
				(!compute_K_ || shblk_norm_D_.is_zero({tile1, tile2})) ? Tile() : shblk_norm_D_.find({tile1, tile2});
		auto norm_D13 =
				(!compute_K_ || shblk_norm_D_.is_zero({tile1, tile3})) ? Tile() : shblk_norm_D_.find({tile1, tile3});

    // grab ptrs to tile data to make addressing more efficient
    auto* F01_ptr = compute_J_ ? F01.data() : nullptr;
    auto* F23_ptr = compute_J_ ? F23.data() : nullptr;
    auto* F02_ptr = compute_K_ ? F02.data() : nullptr;
    auto* F03_ptr = compute_K_ ? F03.data() : nullptr;
    auto* F12_ptr = compute_K_ ? F12.data() : nullptr;
    auto* F13_ptr = compute_K_ ? F13.data() : nullptr;
    const auto* D01_ptr = compute_J_ ? D01.data() : nullptr;
    const auto* D23_ptr = compute_J_ ? D23.data() : nullptr;
    const auto* D02_ptr = compute_K_ ? D02.data() : nullptr;
    const auto* D03_ptr = compute_K_ ? D03.data() : nullptr;
    const auto* D12_ptr = compute_K_ ? D12.data() : nullptr;
    const auto* D13_ptr = compute_K_ ? D13.data() : nullptr;
		const auto* norm_D01_ptr = compute_J_ ? norm_D01.data() : nullptr;
		const auto* norm_D23_ptr = compute_J_ ? norm_D23.data() : nullptr;
		const auto* norm_D02_ptr = compute_K_ ? norm_D02.data() : nullptr;
		const auto* norm_D03_ptr = compute_K_ ? norm_D03.data() : nullptr;
		const auto* norm_D12_ptr = compute_K_ ? norm_D12.data() : nullptr;
		const auto* norm_D13_ptr = compute_K_ ? norm_D13.data() : nullptr;
		const auto t_preLoopShell_1 = mpqc::now();
		dur_preLoopShell_ += mpqc::duration_in_s(t_preLoopShell_0, t_preLoopShell_1);

		const auto t_fock_tile_0 = mpqc::now();
		ExEnv::out0() << "\nI'm the only task.\n" << std::endl;
    // compute contributions to all Fock matrices
    {
			const auto t_preLoop_0 = mpqc::now();
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
			if (same_c0c1 && same_c2c3 && same_c0c2 ) {
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
			const auto nshells1 = cluster1.size();
			const auto nshells2 = cluster2.size();
			const auto nshells3 = cluster3.size();
			ExEnv::out0() << "\n# of shells = \n" << nshells0 << std::endl;

      // this is the index of the first basis functions for each shell *in this shell cluster*
      auto cf0_offset = 0;
      // this is the index of the first basis functions for each shell *in the basis set*
      auto bf0_offset = rng0.first;
			// this is the index of the first shell *in this shell cluster*
			auto sh0 = 0;
			const auto t_preLoop_1 = mpqc::now();
			dur_preLoop_ += mpqc::duration_in_s(t_preLoop_0, t_preLoop_1);

      // loop over unique shell sets
      // N.B. skip nonunique shell sets that did not get eliminated by unique cluster set iteration
      for (const auto& shell0 : cluster0) {
				const auto t_preLoopInts_s0_0 = mpqc::now();
				const auto nf0 = shell0.size();
        auto cf1_offset = 0;
        auto bf1_offset = rng1.first;
				auto sh1 = 0;
				const auto t_preLoopInts_s0_1 = mpqc::now();
				dur_preLoopInts_s0_ += mpqc::duration_in_s(t_preLoopInts_s0_0, t_preLoopInts_s0_1);

				const auto &pl0 = bra_shellpair_list[sh0];
				for (const auto& shell1 : cluster1) {
					const auto nf1 = shell1.size();
					if (std::find(pl0.begin(), pl0.end(), sh1) == pl0.end()) {
						cf1_offset += nf1;
						bf1_offset += nf1;
						sh1++;
						continue;
					}

					const auto t_preLoopInts_s1_0 = mpqc::now();
					// skip if shell set is nonunique
					if (bf0_offset < bf1_offset)
						break;  // assuming basis functions increase monotonically in the basis

					const auto multiplicity01 = bf0_offset == bf1_offset ? 1.0 : 2.0;

					const auto sh01 = sh0 * nshells1 + sh1;  // index of {sh0, sh1} in norm_D01
					const auto Dnorm01 =
							(norm_D01_ptr != nullptr) ? norm_D01_ptr[sh01] : 0.0;

					auto cf2_offset = 0;
					auto bf2_offset = rng2.first;
					auto sh2 = 0;
					const auto t_preLoopInts_s1_1 = mpqc::now();
					dur_preLoopInts_s1_ += mpqc::duration_in_s(t_preLoopInts_s1_0, t_preLoopInts_s1_1);

					for (const auto& shell2 : cluster2) {
						const auto t_preLoopInts_s2_0 = mpqc::now();
						// skip if shell set is nonunique
						if (bf0_offset < bf2_offset) break;

						const auto sh02 = sh0 * nshells2 + sh2;  // index of {sh0, sh2} in norm_D02
						const auto sh12 = sh1 * nshells2 + sh2;  // index of {sh1, sh2} in norm_D12
						const auto Dnorm02 =
								(norm_D02_ptr != nullptr) ? norm_D02_ptr[sh02] : 0.0;
						const auto Dnorm12 =
								(norm_D12_ptr != nullptr) ? norm_D12_ptr[sh12] : 0.0;
						const auto Dnorm012 = std::max({Dnorm12, Dnorm02, Dnorm01});

						const auto nf2 = shell2.size();
						auto cf3_offset = 0;
						auto bf3_offset = rng3.first;
						auto sh3 = 0;
						const auto t_preLoopInts_s2_1 = mpqc::now();
						dur_preLoopInts_s2_ += mpqc::duration_in_s(t_preLoopInts_s2_0, t_preLoopInts_s2_1);

						const auto &pl2 = ket_shellpair_list[sh2];
						for (const auto& shell3 : cluster3) {

							const auto nf3 = shell3.size();

							if (std::find(pl2.begin(), pl2.end(), sh3) == pl2.end()) {
								cf3_offset += nf3;
								bf3_offset += nf3;
								sh3++;
								continue;
							}
//							const auto t_preLoopInts_s3_0 = mpqc::now();


							// skip if shell set is nonunique
							if (bf2_offset < bf3_offset ||
									(bf0_offset == bf2_offset && bf1_offset < bf3_offset))
								break;

							const auto sh03 = sh0 * nshells3 + sh3;  // index of {sh0, sh3} in norm_D03
							const auto sh13 = sh1 * nshells3 + sh3;  // index of {sh1, sh3} in norm_D13
							const auto sh23 = sh2 * nshells3 + sh3;  // index of {sh2, sh3} in norm_D23

							const auto Dnorm03 =
									(norm_D03_ptr != nullptr) ? norm_D03_ptr[sh03] : 0.0;
							const auto Dnorm13 =
									(norm_D13_ptr != nullptr) ? norm_D13_ptr[sh13] : 0.0;
							const auto Dnorm23 =
									(norm_D23_ptr != nullptr) ? norm_D23_ptr[sh23] : 0.0;

							const auto Dnorm0123 = std::max({Dnorm03, Dnorm13, Dnorm23, Dnorm012});

//							const auto t_preLoopInts_s3_1 = mpqc::fenced_now(this->get_world());
//							dur_preLoopInts_s3_ += mpqc::duration_in_s(t_preLoopInts_s3_0, t_preLoopInts_s3_1);

//							const auto t_skip_0 = mpqc::now();
							if (screen.skip(bf0_offset, bf1_offset, bf2_offset, bf3_offset, Dnorm0123)) {
								cf3_offset += nf3;
								bf3_offset += nf3;
								sh3++;
								continue;
							}
//							const auto t_skip_1 = mpqc::now();
//							dur_skip_ += mpqc::duration_in_s(t_skip_0, t_skip_1);

							num_ints_computed_ += nf0 * nf1 * nf2 * nf3;

							const auto multiplicity23 =
									bf2_offset == bf3_offset ? 1.0 : 2.0;
							const auto multiplicity0213 =
									(bf0_offset == bf2_offset && bf1_offset == bf3_offset)
											? 1.0
											: 2.0;
							const auto multiplicity =
									multiplicity01 * multiplicity23 * multiplicity0213;

							const auto t_ints_0 = mpqc::now();
							// compute shell set
							engine.compute2<libint2::Operator::coulomb,
															libint2::BraKet::xx_xx, 0>(
									shell0, shell1, shell2, shell3);
							const auto* eri_0123 = computed_shell_sets[0];
							const auto t_ints_1 = mpqc::now();
							dur_int_compute_ += mpqc::duration_in_s(t_ints_0, t_ints_1);

							if (eri_0123 != nullptr) {  // if the shell set is not screened out

								const auto t_contr_0 = mpqc::now();
								for (auto f0 = 0, f0123 = 0; f0 != nf0; ++f0) {
									const auto cf0 =
											f0 + cf0_offset;  // basis function index in the tile (i.e. shell cluster)
									for (auto f1 = 0; f1 != nf1; ++f1) {
										const auto cf1 = f1 + cf1_offset;
										const auto cf01 =
												cf0 * rng1_size + cf1;  // index of {cf0,cf1} in D01 or F01
										for (auto f2 = 0; f2 != nf2; ++f2) {
											const auto cf2 = f2 + cf2_offset;
											const auto cf02 =
													cf0 * rng2_size + cf2;  // index of {cf0,cf2} in D02 or F02
											const auto cf12 =
													cf1 * rng2_size + cf2;  // index of {cf1,cf2} in D12 or F12
											for (auto f3 = 0; f3 != nf3; ++f3, ++f0123) {
												const auto cf3 = f3 + cf3_offset;
												const auto cf03 =
														cf0 * rng3_size + cf3;  // index of {cf0,cf3} in D03 or F03
												const auto cf13 =
														cf1 * rng3_size + cf3;  // index of {cf1,cf3} in D13 or F13
												const auto cf23 =
														cf2 * rng3_size + cf3;  // index of {cf2,cf3} in D23 or F23

												const auto value = eri_0123[f0123];

												const auto value_scaled_by_multiplicity =
														value * multiplicity;

												if (compute_J_) {
													F01_ptr[cf01] +=
															(D23_ptr != nullptr)
																	? D23_ptr[cf23] * value_scaled_by_multiplicity
																	: 0.0;
													F23_ptr[cf23] +=
															(D01_ptr != nullptr)
																	? D01_ptr[cf01] * value_scaled_by_multiplicity
																	: 0.0;
												}
												if (compute_K_) {
													F02_ptr[cf02] -=
															(D13_ptr != nullptr)
																	? 0.25 * D13_ptr[cf13] *
																				value_scaled_by_multiplicity
																	: 0.0;
													F13_ptr[cf13] -=
															(D02_ptr != nullptr)
																	? 0.25 * D02_ptr[cf02] *
																				value_scaled_by_multiplicity
																	: 0.0;
													F03_ptr[cf03] -=
															(D12_ptr != nullptr)
																	? 0.25 * D12_ptr[cf12] *
																				value_scaled_by_multiplicity
																	: 0.0;
													F12_ptr[cf12] -=
															(D03_ptr != nullptr)
																	? 0.25 * D03_ptr[cf03] *
																				value_scaled_by_multiplicity
																	: 0.0;
												}
											}
										}
									}
								}
								const auto t_contr_1 = mpqc::now();
								dur_contr_ += mpqc::duration_in_s(t_contr_0, t_contr_1);
							}

							cf3_offset += nf3;
							bf3_offset += nf3;
							sh3++;

//							const auto t_s3_1 = mpqc::now();
//							dur_s3_ += mpqc::duration_in_s(t_s3_0, t_s3_1);
						}

						cf2_offset += nf2;
						bf2_offset += nf2;
						sh2++;
					}

					cf1_offset += nf1;
					bf1_offset += nf1;
					sh1++;
				}

        cf0_offset += nf0;
        bf0_offset += nf0;
				sh0++;
      }
    }
		const auto t_fock_tile_1 = mpqc::now();
		dur_fock_tile_ += mpqc::duration_in_s(t_fock_tile_0, t_fock_tile_1);

		const auto t_accumul_0 = mpqc::now();
    const auto me = this->get_world().rank();
    assert(me == pmap_D_->rank());
    // accumulate the contributions by submitting tasks to the owners of their
    // tiles
    if (compute_J_) {
      const auto proc01 = pmap_D_->owner(tile_idx[0] * ntiles + tile_idx[1]);
      const auto proc23 = pmap_D_->owner(tile_idx[2] * ntiles + tile_idx[3]);
      WorldObject_::task(proc01, &FourCenterFockBuilder_::accumulate_task, F01,
                         tile_idx[0], tile_idx[1],
                         madness::TaskAttributes::hipri());
      WorldObject_::task(proc23, &FourCenterFockBuilder_::accumulate_task, F23,
                         tile_idx[2], tile_idx[3],
                         madness::TaskAttributes::hipri());
    }
    if (compute_K_) {
      const auto proc02 = pmap_D_->owner(tile_idx[0] * ntiles + tile_idx[2]);
      const auto proc03 = pmap_D_->owner(tile_idx[0] * ntiles + tile_idx[3]);
      const auto proc12 = pmap_D_->owner(tile_idx[1] * ntiles + tile_idx[2]);
      const auto proc13 = pmap_D_->owner(tile_idx[1] * ntiles + tile_idx[3]);
      WorldObject_::task(proc02, &FourCenterFockBuilder_::accumulate_task, F02,
                         tile_idx[0], tile_idx[2],
                         madness::TaskAttributes::hipri());
      WorldObject_::task(proc03, &FourCenterFockBuilder_::accumulate_task, F03,
                         tile_idx[0], tile_idx[3],
                         madness::TaskAttributes::hipri());
      WorldObject_::task(proc12, &FourCenterFockBuilder_::accumulate_task, F12,
                         tile_idx[1], tile_idx[2],
                         madness::TaskAttributes::hipri());
      WorldObject_::task(proc13, &FourCenterFockBuilder_::accumulate_task, F13,
                         tile_idx[1], tile_idx[3],
                         madness::TaskAttributes::hipri());
    }
		const auto t_accumul_1 = mpqc::now();
		dur_accumul_ += mpqc::duration_in_s(t_accumul_0, t_accumul_1);

	}

	// TODO compute norms in a parallel fashion
	array_type compute_shellblock_norm(const Basis& bs, const array_type& D) const {

		auto& world = this->get_world();
		// make trange1
		const auto& shells_Vec = bs.cluster_shells();
		auto blocking = std::vector<int64_t> {0};
		for (const auto& shells : shells_Vec) {
			const auto nshell = shells.size();
			auto next = blocking.back() + nshell;
			blocking.emplace_back(next);
		}
		auto trange1 = TA::TiledRange1(blocking.begin(), blocking.end());

		auto eig_D = ::mpqc::array_ops::array_to_eigen(D);
		// compute shell block norms
		const auto shells = bs.flattened_shells();
		const auto nshell = shells.size();
		RowMatrixXd norm_D(nshell, nshell);
		for (auto sh1 = 0, sh1_first = 0; sh1 != nshell; ++sh1) {
			const auto sh1_size = shells[sh1].size();
			for (auto sh2 = 0, sh2_first = 0; sh2 != nshell; ++sh2) {
				const auto sh2_size = shells[sh2].size();

				norm_D(sh1, sh2) = eig_D.block(sh1_first, sh2_first, sh1_size, sh2_size).template lpNorm<Eigen::Infinity>();

				sh2_first += sh2_size;
			}

			sh1_first += sh1_size;
		}

		return array_ops::eigen_to_array<Tile, Policy>(world, norm_D, trange1, trange1);
	}

	shellpair_list_t compute_shellpair_list(const ShellVec &shv1, const ShellVec &_shv2 = std::vector<Shell>({Shell()}), double threshold = 1e-12) const {
		const ShellVec &shv2 = ((_shv2.size() == 1 && _shv2[0] == Shell()) ? shv1 : _shv2);
		const auto nsh1 = shv1.size();
		const auto nsh2 = shv2.size();
		const auto shv1_equiv_shv2 = (&shv1 == &shv2);

		// determine max # of primitives in a shell cluster
		auto max_nprim = [](const ShellVec &shv) {
			size_t n = 0;
			for (auto shell : shv)
				n = std::max(shell.nprim(), n);
			return n;
		};
		const auto max_nprim_1 = max_nprim(shv1);
		const auto max_nprim_2 = max_nprim(shv2);

		// determine max angular momentum of a shell cluster
		auto max_l = [](const ShellVec &shv) {
			int l = 0;
			for (auto shell : shv)
				for (auto c: shell.contr)
					l = std::max(c.l, l);
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

				if (significant)
					result[s1].emplace_back(s2);
			}
		}

		// resort shell list in increasing order
		for (auto s1 = 0l; s1 != nsh1; ++s1) {
			auto &list = result[s1];
			std::sort(list.begin(), list.end());
		}

		return result;
	}

};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
