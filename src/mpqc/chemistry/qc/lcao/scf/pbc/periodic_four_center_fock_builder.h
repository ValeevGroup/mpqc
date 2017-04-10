#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class ReferencePeriodicFourCenterFockBuilder
		: public PeriodicFockBuilder<Tile, Policy> {
 public:
	using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;

	ReferencePeriodicFourCenterFockBuilder(Factory &ao_factory)
			: ao_factory_(ao_factory) {}

	~ReferencePeriodicFourCenterFockBuilder() {}

	array_type operator()(array_type const &D, double) override {
		// feed density matrix to Factory
		ao_factory_.set_density(D);

		array_type J, K, G;
		J = ao_factory_.compute_direct(L"(μ ν| J|κ λ)");
		K = ao_factory_.compute_direct(L"(μ ν| K|κ λ)");

		G("mu, nu") = 2.0 * J("mu, nu") - K("mu, nu");

		return G;
	}

	void register_fock(const array_type &fock,
										 FormulaRegistry<array_type> &registry) override {
		registry.insert(Formula(L"(κ|F|λ)"), fock);
	}

 private:
	Factory &ao_factory_;
};

template <typename Tile, typename Policy>
class PeriodicFourCenterFockBuilder
		: public PeriodicFockBuilder<Tile, Policy>,
			public madness::WorldObject<PeriodicFourCenterFockBuilder<Tile, Policy>> {
 public:
	using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
	using WorldObject_ =
			madness::WorldObject<FourCenterFockBuilder<Tile, Policy>>;
	using FourCenterFockBuilder_ = FourCenterFockBuilder<Tile, Policy>;

	using Basis = ::mpqc::lcao::gaussian::Basis;
	using BasisVector = std::vector<Basis>;
	using Shell = typename ::mpqc::lcao::gaussian::Shell;
	using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
	using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
	using func_offset_list = std::unordered_map<size_t, std::tuple<size_t, size_t>>;

	PeriodicFourCenterFockBuilder(madness::World &world,
																std::shared_ptr<const Basis> basis,
																Vector3d &dcell, Vector3i &R_max,
																Vector3i &RJ_max, Vector3i &RD_max,
																bool compute_J, bool compute_K,
																std::string screen = "schwarz",
																double screen_threshold = 1.0e-10)
			: WorldObject_(world),
				compute_J_(compute_J),
				compute_K_(compute_K),
				dcell_(dcell),
				R_max_(R_max),
				RJ_max_(RJ_max),
				RD_max_(RD_max),
				basis0_(basis),
				screen_(screen),
				screen_threshold_(screen_threshold) {
		using ::mpqc::lcao::detail::direct_ord_idx;
		R_size_ = 1 + direct_ord_idx(R_max_(0), R_max_(1), R_max_(2), R_max_);
		RJ_size_ = 1 + direct_ord_idx(RJ_max_(0), RJ_max_(1), RJ_max_(2), RJ_max_);
		RD_size_ = 1 + direct_ord_idx(RD_max_(0), RD_max_(1), RD_max_(2), RD_max_);

		assert(basis0_ != nullptr);
		// WorldObject mandates this is called from the ctor
		WorldObject_::process_pending();
	}

	~ReferencePeriodicFourCenterFockBuilder() {}

	array_type operator()(array_type const &D, double target_precision) override {
		num_ints_computed_ = 0;

		auto vec_G = std::vector<array_type>(RJ_size_, array_type());
		using ::mpqc::lcao::detail::direct_vector;
		using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
		using ::mpqc::lcao::gaussian::detail::make_screener;

		Vector3d zero_shift_base(0.0, 0.0, 0.0);
		auto basisR = shift_basis_origin(*basis0_, zero_shift_base, R_max_, dcell_);
		// make shell block norm of D
		shblk_norm_D_ = compute_shellblock_norm(*basis0_, *basisR, D);

		for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
			auto &G = vec_G[RJ];
			auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);

			auto basisRJ = shift_basis_origin(*basis0_, vec_RJ);
			auto basisRD = shift_basis_origin(*basis0_, vec_RJ, RD_max_, dcell_);

			G = compute_JK_abcd(D, basisR, basisRJ, basisRD, target_precision);
		}

		if (RJ_size_ > 0) {
			for (auto RJ = 1; RJ < RJ_size_; ++RJ) {
				vec_G[0]("mu, nu") += vec_G[RJ]("mu, nu");
			}
		}

		return vec_G[0];
	}

	void register_fock(const array_type &fock,
										 FormulaRegistry<array_type> &registry) override {
		registry.insert(Formula(L"(κ|F|λ)"), fock);
	}

	array_type compute_JK_abcd(array_type &D, std::vector<Basis> basisR,
														 std::vector<Basis> basisRJ,
														 std::vector<Basis> basisRD,
														 double target_precision) const {
		// prepare input data
		auto &compute_world = this->get_world();
		const auto me = compute_world.rank();
		target_precision_ = target_precision;

		// # of tiles per basis
		auto ntiles0 = basis0_->nclusters();
		auto ntilesR = basisR->nclusters();
		auto ntilesRJ = basisRJ->nclusters();
		auto ntilesRD = basisRD->nclusters();

		trange_D_ = D.trange();
		pmap_D_ = D.pmap();

		const auto ntile_tasks =
				static_cast<unit64_t>(ntiles0 * ntilesR * ntilesRD * ntilesRJ);
		auto pmap = std::make_shared<const TA::detail::BlockedPmap>(compute_world,
																																ntile_tasks);

		// make the engine pool
		auto oper_type = libint2::Operator::coulomb;
		if (compute_J_) {
			j_engines_ = ::mpqc::lcao::gaussian::make_engine_pool(
					oper_type,
					utility::make_array_of_refs(*basis0_, *basisR, *basisRJ, *basisRD),
					libint2::BraKet::xx_xx);
			auto bases = ::mpqc::lcao::gaussian::BasisVector{
					{*basis0_, *basisR, *basisRJ, *basisRD}};
			j_p_screener_ = ::mpqc::lcao::gaussian::detail::make_screener(
					compute_world, j_engines_, bases, screen_, screen_threshold_);
		}
		if (compute_K_) {
			k_engines_ = ::mpqc::lcao::gaussian::make_engine_pool(
					oper_type,
					utility::make_array_of_refs(*basis0_, *basisRJ, *basisR, *basisRD),
					libint2::BraKet::xx_xx);
			auto bases = ::mpqc::lcao::gaussian::BasisVector{
					{*basis0_, *basisRJ, *basisR, *basisRD}};
			k_p_screener_ = ::mpqc::lcao::gaussian::detail::make_screener(
					compute_world, k_engines_, bases, screen_, screen_threshold_);
		}

		auto empty = TA::Future<Tile>(Tile());
		for (auto tile0 = 0ul, tile0123 = 0ul; tile0 != ntiles0; ++tile0) {
			for (auto tile1 = 0ul; tile1 != ntilesR; ++tile1) {
				for (auto tile2 = 0ul; tile2 != ntilesRJ; ++tile2) {
					for (auto tile3 = 0ul; tile3 != ntilesRD; ++tile3, ++tile0123) {
						auto D23 = (!compute_J_ || D.is_zero({tile2, tile3}))
													 ? empty
													 : D.find({tile2, tile3});
						auto D13 = (!compute_K_ || D.is_zero({tile1, tile3}))
													 ? empty
													 : D.find({tile1, tile3});
						if (pmap->is_local(tile0123))
							WorldObject_::task(
									me, &PeriodicFourCenterFockBuilder::compute_task_abcd, D23,
									D13, basisR, basisRJ, basisRD,
									std::array<size_t, 4>{{tile0, tile1, tile2, tile3}});
					}
				}
			}
		}

		compute_world.gop.fence();

		// cleanup
		if (compute_J_) j_engines_.reset();
		if (compute_K_) k_engines_.reset();

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
				local_tile_norms.push_back(
						std::make_pair(std::array<size_t, 2>{{i, j}}, ij_norm));
			}
			shape = decltype(shape)(compute_world, local_tile_norms, trange_D_);
		}

		array_type G(compute_world, trange_D_, shape, pmap_D_);

		// copy results of local reduction tasks into G
		for (const auto &local_tile : local_fock_tiles_) {
			// if this tile was not truncated away
			if (!G.shape().is_zero(local_tile.first))
				G.set(local_tile.first, local_tile.second);
		}
		// set the remaining local tiles to 0 (this should only be needed for dense
		// policy)
		G.fill_local(0.0, true);
		local_fock_tiles_.clear();

		return G;
	}

 private:
	// set by ctor
	const bool compute_J_;
	const bool compute_K_;
	std::shared_ptr<const Basis> basis0_;
	const std::string screen_;
	const double screen_threshold_;
	const Vector3d dcell_;
	const Vector3i R_max_;
	const Vector3i RJ_max_;
	const Vector3i RD_max_;
	const int64_t R_size_;
	const int64_t RJ_size_;
	const int64_t RD_size_;

	// mutated by compute_ functions
	mutable std::shared_ptr<lcao::Screener> j_p_screener_;
	mutable std::shared_ptr<lcao::Screener> k_p_screener_;
	mutable madness::ConcurrentHashMap<std::size_t, Tile> local_fock_tiles_;
	mutable TA::TiledRange trange_D_;
	mutable std::shared_ptr<TA::Pmap> pmap_D_;
	mutable double target_precision_ = 0.0;
	mutable ::mpqc::lcao::gaussian::ShrPool<libint2::Engine> j_engines_;
	mutable ::mpqc::lcao::gaussian::ShrPool<libint2::Engine> k_engines_;
	mutable array_type shblk_norm_D_;
	mutable std::atomic<size_t> num_ints_computed_{0};

	void accumulate_task(Tile fock_matrix_tile, long tile0, long tile1) {
		const auto ntiles = trange_D_.dim(1).tile_extent();
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
					[](TA::detail::numeric_t<Tile> &l,
						 const TA::detail::numeric_t<Tile> r) { l += r; },
					size, acc->second.data(), fock_matrix_tile.data());
		}
		acc.release();  // END OF CRITICAL SECTION
	}

	void compute_task_abcd(Tile D23, Tile D13, std::shared_ptr<Basis> basisR,
												 std::shared_ptr<Basis> basisRJ,
												 std::shared_ptr<Basis> basisRD,
												 std::array<size_t, 4> tile_idx) {
		const auto tile0 = tile_idx[0];
		const auto tile1 = tile_idx[1];
		const auto tile2 = tile_idx[2];
		const auto tile3 = tile_idx[3];

		// 1-d tile ranges
		const auto &tr0 = trange_D_.dim(0);
		const auto &tr1 = trange_D_.dim(1);
		const auto ntiles0 = tr0.tile_extent();
		const auto ntiles1 = tr1.tile_extent();
		const auto &rng0 = tr0.tile(tile0);
		const auto &rng1 = tr1.tile(tile1);  // this is for Coulomb term
		const auto &rng2 = tr0.tile(tile2);
		const auto &rng3 = tr1.tile(tile3);
		const auto rng0_size = rng0.second - rng0.first;
		const auto rng1_size = rng1.second - rng1.first;
		const auto rng2_size = rng2.second - rng2.first;
		const auto rng3_size = rng3.second - rng3.first;

		// 2-d tile ranges describing the Fock contribution blocks produced by this
		auto rng23 = compute_J_ ? TA::Range({rng2, rng3}) : TA::Range();
		auto rng13 = compute_K_ ? TA::Range({rng1, rng3}) : TA::Range();

		// initialize contribution to the Fock matrices
		auto F23 = compute_J_ ? Tile(std::move(rng23), 0.0) : Tile();
		auto F13 = compute_K_ ? Tile(std::move(rng13), 0.0) : Tile();

		// find shell block norm of D
		auto norm_D23 = (!compute_J_ || shblk_norm_D_.is_zero({tile2, tile3}))
												? Tile()
												: shblk_norm_D_.find({tile2, tile3});
		auto norm_D13 = (!compute_K_ || shblk_norm_D_.is_zero({tile1, tile3}))
												? Tile()
												: shblk_norm_D_.find({tile1, tile3});

		// grab ptrs to tile data to make addressing more efficient
		auto *F23_ptr = compute_J_ ? F23.data() : nullptr;
		auto *F13_ptr = compute_K_ ? F13.data() : nullptr;
		const auto *D23_ptr = compute_J_ ? D23.data() : nullptr;
		const auto *D13_ptr = compute_K_ ? D13.data() : nullptr;
		const auto *norm_D23_ptr = compute_J_ ? norm_D23.data() : nullptr;
		const auto *norm_D13_ptr = compute_K_ ? norm_D13.data() : nullptr;

		// compute contributions to all Fock matrices
		{
			const auto engine_precision = target_precision_;

			auto &j_screen = *j_p_screener_;
			auto j_engine = j_engines_->local();
			j_engine.set_precision(engine_precision);
			const auto &j_computed_shell_sets = j_engine.results();

			auto &k_screen = *k_p_screener_;
			auto k_engine = k_engines_->local();
			k_engine.set_precision(engine_precision);
			const auto &k_computed_shell_sets = k_engine.results();

			const auto &basis0 = basis0_;
			// shell clusters for this tile
			if (compute_J_) {
				const auto& cluster0 = basis0->cluster_shells()[tile0];
				const auto& cluster1 = basisR->cluster_shells()[tile1];
				const auto& cluster2 = basisRJ->cluster_shells()[tile2];
				const auto& cluster3 = basisRD->cluster_shells()[tile3];
			} else {
				const auto& cluster0 = basis0->cluster_shells()[tile0];
				const auto& cluster1 = basisRJ->cluster_shells()[tile2];
				const auto& cluster2 = basisR->cluster_shells()[tile1];
				const auto& cluster3 = basisRD->cluster_shells()[tile3];
			}

			// make unique shell pair list
			shellpair_list_t bra_shellpair_list, ket_shellpair_list;
			bra_shellpair_list = compute_shellpair_list(cluster0, cluster1);
			ket_shellpair_list = compute_shellpair_list(cluster2, cluster3);

			// number of shells in each cluster
			const auto nshells0 = cluster0.size();
			const auto nshells1 = cluster1.size();
			const auto nshells2 = cluster2.size();
			const auto nshells3 = cluster3.size();

			// compute offset list of cluster1 and cluster3
			auto offset_list_c1 = compute_func_offset_list(cluster1, rng1.first);
			auto offset_list_c3 = compute_func_offset_list(cluster3, rng3.first);

			// this is the index of the first basis functions for each shell *in this shell cluster*
			auto cf0_offset = 0;
			// this is the index of the first basis functions for each shell *in the basis set*
			auto bf0_offset = rng0.first;

			size_t cf1_offset, bf1_offset, cf3_offset, bf3_offset;

			// loop over unique shell sets
			// N.B. skip nonunique shell sets that did not get eliminated by unique cluster set iteration
			for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
				const auto& shell0 = cluster0[sh0];
				const auto nf0 = shell0.size();

				for (const auto &sh1 : bra_shellpair_list[sh0]) {
					std::tie(cf1_offset, bf1_offset) = offset_list_c1[sh1];
					// skip if shell set is nonunique
//					if (bf0_offset < bf1_offset)
//						break;  // assuming basis functions increase monotonically in the basis

					const auto &shell1 = cluster1[sh1];
					const auto nf1 = shell1.size();

//					const auto multiplicity01 = bf0_offset == bf1_offset ? 1.0 : 2.0;

//					const auto sh01 = sh0 * nshells1 + sh1;  // index of {sh0, sh1} in norm_D01
//					const auto Dnorm01 =
//							(norm_D01_ptr != nullptr) ? norm_D01_ptr[sh01] : 0.0;

					auto cf2_offset = 0;
					auto bf2_offset = rng2.first;

					for (auto sh2 = 0; sh2 != nshells2; ++sh2) {
						// skip if shell set is nonunique
//						if (bf0_offset < bf2_offset) break;

						const auto &shell2 = cluster2[sh2];
						const auto nf2 = shell2.size();

//						const auto sh02 = sh0 * nshells2 + sh2;  // index of {sh0, sh2} in norm_D02
//						const auto sh12 = sh1 * nshells2 + sh2;  // index of {sh1, sh2} in norm_D12
//						const auto Dnorm02 =
//								(norm_D02_ptr != nullptr) ? norm_D02_ptr[sh02] : 0.0;
//						const auto Dnorm12 =
//								(norm_D12_ptr != nullptr) ? norm_D12_ptr[sh12] : 0.0;
//						const auto Dnorm012 = std::max({Dnorm12, Dnorm02, Dnorm01});

						for (const auto &sh3 : ket_shellpair_list[sh2]) {
							std::tie(cf3_offset, bf3_offset) = offset_list_c3[sh3];
							// skip if shell set is nonunique
//							if (bf2_offset < bf3_offset ||
//									(bf0_offset == bf2_offset && bf1_offset < bf3_offset))
//								break;

							const auto &shell3 = cluster3[sh3];
							const auto nf3 = shell3.size();

//							const auto sh03 = sh0 * nshells3 + sh3;  // index of {sh0, sh3} in norm_D03
							const auto sh13 = sh1 * nshells3 + sh3;  // index of {sh1, sh3} in norm_D13
							const auto sh23 = sh2 * nshells3 + sh3;  // index of {sh2, sh3} in norm_D23
//							const auto Dnorm03 =
//									(norm_D03_ptr != nullptr) ? norm_D03_ptr[sh03] : 0.0;
							const auto Dnorm13 =
									(norm_D13_ptr != nullptr) ? norm_D13_ptr[sh13] : 0.0;
							const auto Dnorm23 =
									(norm_D23_ptr != nullptr) ? norm_D23_ptr[sh23] : 0.0;
//							const auto Dnorm0123 = std::max({Dnorm03, Dnorm13, Dnorm23, Dnorm012});
							const auto Dnorm0123 = std::max({Dnorm13, Dnorm23});


							if (compute_J_ && j_screen.skip(bf0_offset, bf1_offset, bf2_offset, bf3_offset, Dnorm0123))
								continue;

							num_ints_computed_ += nf0 * nf1 * nf2 * nf3;

//							const auto multiplicity23 =
//									bf2_offset == bf3_offset ? 1.0 : 2.0;
//							const auto multiplicity0213 =
//									(bf0_offset == bf2_offset && bf1_offset == bf3_offset)
//											? 1.0
//											: 2.0;
//							const auto multiplicity =
//									multiplicity01 * multiplicity23 * multiplicity0213;

							// compute shell set
							j_engine.compute2<libint2::Operator::coulomb,
															libint2::BraKet::xx_xx, 0>(
									shell0, shell1, shell2, shell3);
							const auto* eri_0123 = computed_shell_sets[0];

							if (eri_0123 != nullptr) {  // if the shell set is not screened out

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

//												const auto value_scaled_by_multiplicity =
//														value * multiplicity;

												if (compute_J_) {
													F01_ptr[cf01] +=
															(D23_ptr != nullptr)
																	? D23_ptr[cf23] * value
																	: 0.0;
//													F23_ptr[cf23] +=
//															(D01_ptr != nullptr)
//																	? D01_ptr[cf01] * value_scaled_by_multiplicity
//																	: 0.0;
												}
												if (compute_K_) {
													F02_ptr[cf02] -=
															(D13_ptr != nullptr)
																	? 0.25 * D13_ptr[cf13] *
																				value_scaled_by_multiplicity
																	: 0.0;
//													F13_ptr[cf13] -=
//															(D02_ptr != nullptr)
//																	? 0.25 * D02_ptr[cf02] *
//																				value_scaled_by_multiplicity
//																	: 0.0;
//													F03_ptr[cf03] -=
//															(D12_ptr != nullptr)
//																	? 0.25 * D12_ptr[cf12] *
//																				value_scaled_by_multiplicity
//																	: 0.0;
//													F12_ptr[cf12] -=
//															(D03_ptr != nullptr)
//																	? 0.25 * D03_ptr[cf03] *
//																				value_scaled_by_multiplicity
//																	: 0.0;
												}
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
	}

	// TODO compute norms in a parallel fashion
	/*!
	 * \brief This computes shell-block norm of density matrix \c D
	 * \param bs Basis
	 * \param D density matrix
	 * \return
	 */
	array_type compute_shellblock_norm(const Basis &bs0, const Basis &bs1,
																		 const array_type &D) const {
		auto &world = this->get_world();
		TA::TiledRange1 tr0, tr1;
		// make trange1 for basis0
		{
			const auto &shells_Vec = bs0.cluster_shells();
			auto blocking = std::vector<int64_t>{0};
			for (const auto &shells : shells_Vec) {
				const auto nshell = shells.size();
				auto next = blocking.back() + nshell;
				blocking.emplace_back(next);
			}
			tr0 = TA::TiledRange1(blocking.begin(), blocking.end());
		}
		// make trange1 for basis1
		{
			const auto &shells_Vec = bs1.cluster_shells();
			auto blocking = std::vector<int64_t>{0};
			for (const auto &shells : shells_Vec) {
				const auto nshell = shells.size();
				auto next = blocking.back() + nshell;
				blocking.emplace_back(next);
			}
			tr1 = TA::TiledRange1(blocking.begin(), blocking.end());
		}

		auto eig_D = ::mpqc::array_ops::array_to_eigen(D);
		// compute shell block norms
		const auto shells0 = bs0.flattened_shells();
		const auto nshell0 = shells0.size();
		const auto shells1 = bs1.flattened_shells();
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
			double threshold = 1e-12) const {
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
	func_offset_list compute_func_offset_list(const ShellVec& cluster,
																						const size_t bf_first) const {
		func_offset_list result;

		auto cf_offset = 0;
		auto bf_offset = bf_first;

		const auto nshell = cluster.size();
		for (auto s = 0; s != nshell; ++s) {
			const auto& shell = cluster[s];
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

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
