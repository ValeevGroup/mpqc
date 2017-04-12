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
	using const_data_ptr = typename Tile::allocator_type::const_pointer;

	using WorldObject_ =
			madness::WorldObject<PeriodicFourCenterFockBuilder<Tile, Policy>>;
	using PeriodicFourCenterFockBuilder_ = PeriodicFourCenterFockBuilder<Tile, Policy>;

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
																int64_t R_size, int64_t RJ_size,
																int64_t RD_size,
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
				R_size_(R_size),
				RJ_size_(RJ_size),
				RD_size_(RD_size),
				basis0_(basis),
				screen_(screen),
				screen_threshold_(screen_threshold) {

		assert(basis0_ != nullptr && "No basis is provided");
		assert((compute_J_ || compute_K_) && "No Coulomb && No Exchange");
		// WorldObject mandates this is called from the ctor
		WorldObject_::process_pending();
	}

	~PeriodicFourCenterFockBuilder() {}

	array_type operator()(array_type const &D, double target_precision) override {
		J_num_ints_computed_ = 0;
		K_num_ints_computed_ = 0;

		auto vec_G = std::vector<array_type>(RJ_size_, array_type());
		using ::mpqc::lcao::detail::direct_vector;
		using ::mpqc::lcao::gaussian::detail::make_screener;
		using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

		Vector3d zero_shift_base(0.0, 0.0, 0.0);
		basisR_ = shift_basis_origin(*basis0_, zero_shift_base,
																					 R_max_, dcell_);

		for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
			auto &G = vec_G[RJ];
			auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);

			basisRJ_ = shift_basis_origin(*basis0_, vec_RJ);
			basisRD_ = shift_basis_origin(*basis0_, vec_RJ, RD_max_, dcell_);

			G = compute_JK_abcd(D, target_precision);
		}

		if (RJ_size_ > 0) {
			for (auto RJ = 1; RJ < RJ_size_; ++RJ) {
				vec_G[0]("mu, nu") += vec_G[RJ]("mu, nu");
			}
		}

		ExEnv::out0() << "\n# of integrals in Coulomb = "
									<< J_num_ints_computed_ << std::endl;
		ExEnv::out0() << "# of integrals in Exchange = "
									<< K_num_ints_computed_ << std::endl;

		return vec_G[0];
	}

	void register_fock(const array_type &fock,
										 FormulaRegistry<array_type> &registry) override {
		registry.insert(Formula(L"(κ|F|λ)"), fock);
	}

	array_type compute_JK_abcd(array_type const &D,
														 double target_precision) const {
		// prepare input data
		auto &compute_world = this->get_world();
		const auto me = compute_world.rank();
		target_precision_ = target_precision;

		const auto basis0 = *basis0_;
		const auto basisR = *basisR_;
		const auto basisRJ = *basisRJ_;
		const auto basisRD = *basisRD_;

		// # of tiles per basis
		auto ntiles0 = basis0.nclusters();
		auto ntilesR = basisR.nclusters();
		auto ntilesRJ = basisRJ.nclusters();
		auto ntilesRD = basisRD.nclusters();

		trange_D_ = D.trange();
		pmap_D_ = D.pmap();

		const auto ntile_tasks =
				static_cast<uint64_t>(ntiles0 * ntilesR * ntilesRD * ntilesRJ);
		auto pmap = std::make_shared<const TA::detail::BlockedPmap>(compute_world,
																																ntile_tasks);

		// make the engine pool
		auto oper_type = libint2::Operator::coulomb;
		if (compute_J_) {
			j_engines_ = ::mpqc::lcao::gaussian::make_engine_pool(
					oper_type,
					utility::make_array_of_refs(basis0, basisR, basisRJ, basisRD),
					libint2::BraKet::xx_xx);
			auto bases = ::mpqc::lcao::gaussian::BasisVector{
					{basis0, basisR, basisRJ, basisRD}};
			j_p_screener_ = ::mpqc::lcao::gaussian::detail::make_screener(
					compute_world, j_engines_, bases, screen_, screen_threshold_);
		}
		if (compute_K_) {
			k_engines_ = ::mpqc::lcao::gaussian::make_engine_pool(
					oper_type,
					utility::make_array_of_refs(basis0, basisRJ, basisR, basisRD),
					libint2::BraKet::xx_xx);
			auto bases = ::mpqc::lcao::gaussian::BasisVector{
					{basis0, basisRJ, basisR, basisRD}};
			k_p_screener_ = ::mpqc::lcao::gaussian::detail::make_screener(
					compute_world, k_engines_, bases, screen_, screen_threshold_);
		}

		auto empty = TA::Future<Tile>(Tile());
		for (auto tile0 = 0ul, tile0123 = 0ul; tile0 != ntiles0; ++tile0) {
			for (auto tile1 = 0ul; tile1 != ntilesR; ++tile1) {
				for (auto tile2 = 0ul; tile2 != ntilesRJ; ++tile2) {
					for (auto tile3 = 0ul; tile3 != ntilesRD; ++tile3, ++tile0123) {
						auto D_RJRD = (D.is_zero({tile2, tile3}))
													 ? empty
													 : D.find({tile2, tile3});
						if (pmap->is_local(tile0123))
							WorldObject_::task(
									me, &PeriodicFourCenterFockBuilder::compute_task_abcd, D_RJRD,
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
				const auto i = ij / ntilesR;
				const auto j = ij % ntilesR;
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
	mutable std::shared_ptr<Basis> basisR_;
	mutable std::shared_ptr<Basis> basisRJ_;
	mutable std::shared_ptr<Basis> basisRD_;
	mutable madness::ConcurrentHashMap<std::size_t, Tile> local_fock_tiles_;
	mutable TA::TiledRange trange_D_;
	mutable std::shared_ptr<TA::Pmap> pmap_D_;
	mutable double target_precision_ = 0.0;
	mutable ::mpqc::lcao::gaussian::ShrPool<libint2::Engine> j_engines_;
	mutable ::mpqc::lcao::gaussian::ShrPool<libint2::Engine> k_engines_;
	mutable array_type shblk_norm_D_;
	mutable std::atomic<size_t> J_num_ints_computed_{0};
	mutable std::atomic<size_t> K_num_ints_computed_{0};

	void accumulate_task(Tile fock_matrix_tile, long tile0, long tile1) {
		const auto ntiles1 = trange_D_.dim(1).tile_extent();
		const auto tile01 = tile0 * ntiles1 + tile1;
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

	void compute_task_abcd(Tile D_RJRD, std::array<size_t, 4> tile_idx) {
		const auto tile0 = tile_idx[0];
		const auto tileR = tile_idx[1];
		const auto tileRJ = tile_idx[2];
		const auto tileRD = tile_idx[3];

		// 1-d tile ranges
		const auto &tr0 = trange_D_.dim(0);
		const auto &tr1 = trange_D_.dim(1);
		const auto ntiles0 = tr0.tile_extent();
		const auto ntiles1 = tr1.tile_extent();
		const auto &rng0 = tr0.tile(tile0);
		const auto &rngR = tr1.tile(tileR);
		const auto &rngRJ = tr0.tile(tileRJ);
		const auto &rngRD = tr1.tile(tileRD);
		const auto rng0_size = rng0.second - rng0.first;
		const auto rngR_size = rngR.second - rngR.first;
		const auto rngRJ_size = rngRJ.second - rngRJ.first;
		const auto rngRD_size = rngRD.second - rngRD.first;

		// 2-d tile ranges describing the Fock contribution blocks produced by this
//		auto rng0R = compute_J_ ? TA::Range({rng0, rngR}) : TA::Range();
		auto rng0R = TA::Range({rng0, rngR});

		// initialize contribution to the Fock matrices
//		auto F0R = compute_J_ ? Tile(std::move(rng0R), 0.0) : Tile();
//		auto F01K = compute_K_ ? Tile(std::move(rng0R), 0.0) : Tile();
		auto F0R = Tile(std::move(rng0R), 0.0);

		// grab ptrs to tile data to make addressing more efficient
//		auto *F01_ptr = compute_J_ ? F0R.data() : nullptr;
//		auto *F02_ptr = compute_K_ ? F01K.data() : nullptr;
//		const auto *D23_ptr = compute_J_ ? D_RJRD.data() : nullptr;
		auto *F0R_ptr = F0R.data();
		const auto *D_RJRD_ptr = D_RJRD.data();

		// compute Coulomb and/or Exchange contributions to all Fock matrices
		{
			const auto engine_precision = target_precision_;

			const auto &basis0 = basis0_;
			const auto &basisR = basisR_;
			const auto &basisRJ = basisRJ_;
			const auto &basisRD = basisRD_;

			// shell clusters for this tile
			const auto& cluster0 = basis0->cluster_shells()[tile0];
			const auto& clusterR = basisR->cluster_shells()[tileR];
			const auto& clusterRJ = basisRJ->cluster_shells()[tileRJ];
			const auto& clusterRD = basisRD->cluster_shells()[tileRD];

			// compute shell block norms of density tiles
			auto norm_D_RJRD =
					(D_RJRD_ptr == nullptr)
							? Tile()
							: compute_shellblock_norm_of_tile(D_RJRD_ptr, clusterRJ, clusterRD,
																								rngRJ_size, rngRD_size);

			// grab ptrs to tile data to make addressing more efficient
			const auto* norm_D_RJRD_ptr = norm_D_RJRD.data();

			// number of shells in each cluster
			const auto nshells0 = cluster0.size();
			const auto nshellsR = clusterR.size();
			const auto nshellsRJ = clusterRJ.size();
			const auto nshellsRD = clusterRD.size();

			// compute Coulomb contributions to all Fock matrices
			if (compute_J_) {

				auto &screen = *j_p_screener_;
				auto engine = j_engines_->local();
				engine.set_precision(engine_precision);
				const auto &computed_shell_sets = engine.results();

				// make non-negligible shell pair list
				shellpair_list_t bra_shellpair_list, ket_shellpair_list;
				bra_shellpair_list = compute_shellpair_list(cluster0, clusterR);
				ket_shellpair_list = compute_shellpair_list(clusterRJ, clusterRD);

				// compute offset list of cluster1 and cluster3
				auto offset_list_bra1 = compute_func_offset_list(clusterR, rngR.first);
				auto offset_list_ket1 = compute_func_offset_list(clusterRD, rngRD.first);

				// this is the index of the first basis functions for each shell *in this shell cluster*
				auto cf0_offset = 0;
				// this is the index of the first basis functions for each shell *in the basis set*
				auto bf0_offset = rng0.first;
				size_t cf1_offset, bf1_offset, cf3_offset, bf3_offset;

				// loop over all shell sets
				for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
					const auto& shell0 = cluster0[sh0];
					const auto nf0 = shell0.size();

					for (const auto &sh1 : bra_shellpair_list[sh0]) {
						std::tie(cf1_offset, bf1_offset) = offset_list_bra1[sh1];
						// skip if shell set is nonunique
	//					if (bf0_offset < bf1_offset)
	//						break;  // assuming basis functions increase monotonically in the basis

						const auto &shell1 = clusterR[sh1];
						const auto nf1 = shell1.size();

	//					const auto multiplicity01 = bf0_offset == bf1_offset ? 1.0 : 2.0;

	//					const auto sh01 = sh0 * nshells1 + sh1;  // index of {sh0, sh1} in norm_D01
	//					const auto Dnorm01 =
	//							(norm_D01_ptr != nullptr) ? norm_D01_ptr[sh01] : 0.0;

						auto cf2_offset = 0;
						auto bf2_offset = rngRJ.first;

						for (auto sh2 = 0; sh2 != nshellsRJ; ++sh2) {
							// skip if shell set is nonunique
	//						if (bf0_offset < bf2_offset) break;

							const auto &shell2 = clusterRJ[sh2];
							const auto nf2 = shell2.size();

	//						const auto sh02 = sh0 * nshells2 + sh2;  // index of {sh0, sh2} in norm_D02
	//						const auto sh12 = sh1 * nshells2 + sh2;  // index of {sh1, sh2} in norm_D12
	//						const auto Dnorm02 =
	//								(norm_D02_ptr != nullptr) ? norm_D02_ptr[sh02] : 0.0;
	//						const auto Dnorm12 =
	//								(norm_D12_ptr != nullptr) ? norm_D12_ptr[sh12] : 0.0;
	//						const auto Dnorm012 = std::max({Dnorm12, Dnorm02, Dnorm01});

							for (const auto &sh3 : ket_shellpair_list[sh2]) {
								std::tie(cf3_offset, bf3_offset) = offset_list_ket1[sh3];
								// skip if shell set is nonunique
	//							if (bf2_offset < bf3_offset ||
	//									(bf0_offset == bf2_offset && bf1_offset < bf3_offset))
	//								break;

								const auto &shell3 = clusterRD[sh3];
								const auto nf3 = shell3.size();

	//							const auto sh03 = sh0 * nshells3 + sh3;  // index of {sh0, sh3} in norm_D03
	//							const auto sh13 = sh1 * nshells3 + sh3;  // index of {sh1, sh3} in norm_D13
								const auto sh23 = sh2 * nshellsRD + sh3;  // index of {sh2, sh3} in norm_D23
	//							const auto Dnorm03 =
	//									(norm_D03_ptr != nullptr) ? norm_D03_ptr[sh03] : 0.0;
	//							const auto Dnorm13 =
	//									(norm_D13_ptr != nullptr) ? norm_D13_ptr[sh13] : 0.0;
								const auto Dnorm23 =
										(norm_D_RJRD_ptr != nullptr) ? norm_D_RJRD_ptr[sh23] : 0.0;
	//							const auto Dnorm0123 = std::max({Dnorm03, Dnorm13, Dnorm23, Dnorm012});
	//							const auto Dnorm0123 = std::max({Dnorm13, Dnorm23});


								if (screen.skip(bf0_offset, bf1_offset, bf2_offset, bf3_offset, Dnorm23))
									continue;

								J_num_ints_computed_ += nf0 * nf1 * nf2 * nf3;

	//							const auto multiplicity23 =
	//									bf2_offset == bf3_offset ? 1.0 : 2.0;
	//							const auto multiplicity0213 =
	//									(bf0_offset == bf2_offset && bf1_offset == bf3_offset)
	//											? 1.0
	//											: 2.0;
	//							const auto multiplicity =
	//									multiplicity01 * multiplicity23 * multiplicity0213;

								// compute shell set
								engine.compute2<libint2::Operator::coulomb,
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
													cf0 * rngR_size + cf1;  // index of {cf0,cf1} in F0R
											for (auto f2 = 0; f2 != nf2; ++f2) {
												const auto cf2 = f2 + cf2_offset;
//												const auto cf02 =
//														cf0 * rngRJ_size + cf2;  // index of {cf0,cf2} in D02 or F02
//												const auto cf12 =
//														cf1 * rngRJ_size + cf2;  // index of {cf1,cf2} in D12 or F12
												for (auto f3 = 0; f3 != nf3; ++f3, ++f0123) {
													const auto cf3 = f3 + cf3_offset;
//													const auto cf03 =
//															cf0 * rngRD_size + cf3;  // index of {cf0,cf3} in D03 or F03
//													const auto cf13 =
//															cf1 * rngRD_size + cf3;  // index of {cf1,cf3} in D13 or F13
													const auto cf23 =
															cf2 * rngRD_size + cf3;  // index of {cf2,cf3} in D_RJRD

													const auto value = eri_0123[f0123];

	//												const auto value_scaled_by_multiplicity =
	//														value * multiplicity;

													F0R_ptr[cf01] +=
															(D_RJRD_ptr != nullptr)
																	? 2.0 * D_RJRD_ptr[cf23] * value
																	: 0.0;

//													if (compute_J_) {
//														F0R_ptr[cf01] +=
//																(D_RJRD_ptr != nullptr)
//																		? 2.0 * D_RJRD_ptr[cf23] * value
//																		: 0.0;
//	//													F23_ptr[cf23] +=
//	//															(D01_ptr != nullptr)
//	//																	? D01_ptr[cf01] * value_scaled_by_multiplicity
//	//																	: 0.0;
//													}
//													if (compute_K_) {
//														F02_ptr[cf02] -=
//																(D13_ptr != nullptr)
//																		? 0.25 * D13_ptr[cf13] *
//																					value_scaled_by_multiplicity
//																		: 0.0;
//	//													F13_ptr[cf13] -=
//	//															(D02_ptr != nullptr)
//	//																	? 0.25 * D02_ptr[cf02] *
//	//																				value_scaled_by_multiplicity
//	//																	: 0.0;
//	//													F03_ptr[cf03] -=
//	//															(D12_ptr != nullptr)
//	//																	? 0.25 * D12_ptr[cf12] *
//	//																				value_scaled_by_multiplicity
//	//																	: 0.0;
//	//													F12_ptr[cf12] -=
//	//															(D03_ptr != nullptr)
//	//																	? 0.25 * D03_ptr[cf03] *
//	//																				value_scaled_by_multiplicity
//	//																	: 0.0;
//													}
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

			if (compute_K_) {
				auto &screen = *k_p_screener_;
				auto engine = k_engines_->local();
				engine.set_precision(engine_precision);
				const auto &computed_shell_sets = engine.results();

				// make non-negligible shell pair list
				shellpair_list_t bra_shellpair_list, ket_shellpair_list;
				bra_shellpair_list = compute_shellpair_list(cluster0, clusterRJ);
				ket_shellpair_list = compute_shellpair_list(clusterR, clusterRD);

				// compute offset list of cluster1 and cluster3
				auto offset_list_bra1 = compute_func_offset_list(clusterRJ, rngRJ.first);
				auto offset_list_ket1 = compute_func_offset_list(clusterRD, rngRD.first);

				// this is the index of the first basis functions for each shell *in this shell cluster*
				auto cf0_offset = 0;
				// this is the index of the first basis functions for each shell *in the basis set*
				auto bf0_offset = rng0.first;
				size_t cf1_offset, bf1_offset, cf3_offset, bf3_offset;

				// loop over all shell sets
				for (auto sh0 = 0; sh0 != nshells0; ++sh0) {
					const auto& shell0 = cluster0[sh0];
					const auto nf0 = shell0.size();

					for (const auto &sh1 : bra_shellpair_list[sh0]) {
						std::tie(cf1_offset, bf1_offset) = offset_list_bra1[sh1];

						const auto &shell1 = clusterRJ[sh1];
						const auto nf1 = shell1.size();

						auto cf2_offset = 0;
						auto bf2_offset = rngR.first;

						for (auto sh2 = 0; sh2 != nshellsR; ++sh2) {
							const auto &shell2 = clusterR[sh2];
							const auto nf2 = shell2.size();

							for (const auto &sh3 : ket_shellpair_list[sh2]) {
								std::tie(cf3_offset, bf3_offset) = offset_list_ket1[sh3];

								const auto &shell3 = clusterRD[sh3];
								const auto nf3 = shell3.size();

								const auto sh23 = sh2 * nshellsRD + sh3;
								const auto Dnorm23 =
										(norm_D_RJRD_ptr != nullptr) ? norm_D_RJRD_ptr[sh23] : 0.0;

								if (screen.skip(bf0_offset, bf1_offset, bf2_offset, bf3_offset, Dnorm23))
									continue;

								K_num_ints_computed_ += nf0 * nf1 * nf2 * nf3;

								// compute shell set
								engine.compute2<libint2::Operator::coulomb,
																libint2::BraKet::xx_xx, 0>(
											shell0, shell1, shell2, shell3);
								const auto* eri_0123 = computed_shell_sets[0];

								if (eri_0123 != nullptr) {  // if the shell set if not screened out

									for (auto f0 = 0, f0123 = 0; f0 != nf0; ++f0) {
										const auto cf0 = f0 + cf0_offset;  // basis function index in the tile (i.e. shell cluster)
										for (auto f1 = 0; f1 != nf1; ++f1) {
											const auto cf1 = f1 + cf1_offset;
											for (auto f2 = 0; f2 != nf2; ++f2) {
												const auto cf2 = f2 + cf2_offset;
												const auto cf02 = cf0 * rngR_size + cf2;  // index of {cf0,cf2} in F0R
												for (auto f3 = 0; f3 != nf3; ++f3, ++f0123) {
													const auto cf3 = f3 + cf3_offset;
													const auto cf13 = cf1 * rngRD_size + cf3;  // index of {cf1, cf3} in D_RJRD

													const auto value = eri_0123[f0123];

													F0R_ptr[cf02] -=
															(D_RJRD_ptr != nullptr)
																	? D_RJRD_ptr[cf13] * value
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

		}

		const auto me = this->get_world().rank();
		assert(me == pmap_D_->rank());
		// accumulate the contributions by submitting tasks to the owners of their
		// tiles
		{
			const auto proc = pmap_D_->owner(tile0 * ntiles1 + tileR);
			WorldObject_::task(proc, &PeriodicFourCenterFockBuilder_::accumulate_task,
												 F0R, tile0, tileR,
												 madness::TaskAttributes::hipri());
		}

//		if (compute_J_) {
//			const auto proc01 = pmap_D_->owner(tile_idx[0] * ntiles1 + tile_idx[1]);
//			const auto proc23 = pmap_D_->owner(tile_idx[2] * ntiles3 + tile_idx[3]);
//			WorldObject_::task(proc01, &PeriodicFourCenterFockBuilder_::accumulate_task, F0R,
//												 tile_idx[0], tile_idx[1],
//												 madness::TaskAttributes::hipri());
//			WorldObject_::task(proc23, &FourCenterFockBuilder_::accumulate_task, F23,
//												 tile_idx[2], tile_idx[3],
//												 madness::TaskAttributes::hipri());
//		}
//		if (compute_K_) {
//			const auto proc02 = pmap_D_->owner(tile_idx[0] * ntiles2 + tile_idx[2]);
//			const auto proc03 = pmap_D_->owner(tile_idx[0] * ntiles + tile_idx[3]);
//			const auto proc12 = pmap_D_->owner(tile_idx[1] * ntiles + tile_idx[2]);
//			const auto proc13 = pmap_D_->owner(tile_idx[1] * ntiles + tile_idx[3]);
//			WorldObject_::task(proc02, &PeriodicFourCenterFockBuilder_::accumulate_task, F01K,
//												 tile_idx[0], tile_idx[2],
//												 madness::TaskAttributes::hipri());
//			WorldObject_::task(proc03, &FourCenterFockBuilder_::accumulate_task, F03,
//												 tile_idx[0], tile_idx[3],
//												 madness::TaskAttributes::hipri());
//			WorldObject_::task(proc12, &FourCenterFockBuilder_::accumulate_task, F12,
//												 tile_idx[1], tile_idx[2],
//												 madness::TaskAttributes::hipri());
//			WorldObject_::task(proc13, &FourCenterFockBuilder_::accumulate_task, F13,
//												 tile_idx[1], tile_idx[3],
//												 madness::TaskAttributes::hipri());
//		}


	}

	/*!
	 * \brief This computes shell block norm of a density tile
	 * \param arg_tile_ptr a const pointer to data of the tile
	 * \param cluster0 a shell cluster (a.k.a. std::vector<Shell>)
	 * \param cluster1 a shell cluster (a.k.a. std::vector<Shell>)
	 * \param nf0 number of functions in \c cluster0
	 * \param nf1 number of functions in \c cluster1
	 * \return a tile filled with shell block norms
	 */
	Tile compute_shellblock_norm_of_tile(const_data_ptr arg_tile_ptr,
																			 const ShellVec& cluster0,
																			 const ShellVec& cluster1,
																			 const size_t nf0,
																			 const size_t nf1) const {
		// # of shells in this shell cluster
		const auto nsh0 = cluster0.size();
		const auto nsh1 = cluster1.size();

		auto range01 = TA::Range({nsh0, nsh1});
		auto shblk_norm = Tile(std::move(range01), 0.0);
		auto* shblk_norm_ptr = shblk_norm.data();

		const auto arg_tile_eig =
				Eigen::Map<const RowMatrixXd>(arg_tile_ptr, nf0, nf1);

		for (auto s0 = 0, s0_first = 0, s01 = 0; s0 != nsh0; ++s0) {
			const auto s0_size = cluster0[s0].size();
			for (auto s1 = 0, s1_first = 0; s1 != nsh1; ++s1, ++s01) {
				const auto s1_size = cluster1[s1].size();
				shblk_norm_ptr[s01] =
						arg_tile_eig.block(s0_first, s1_first, s0_size, s1_size)
								.template lpNorm<Eigen::Infinity>();
				s1_first += s1_size;
			}
			s0_first += s0_size;
		}

		return shblk_norm;
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
