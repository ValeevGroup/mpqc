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
						auto D23 =
								(!compute_J_ || D.is_zero({tile2, tile3})) ? empty : D.find({tile2, tile3});
						auto D13 =
								(!compute_K_ || D.is_zero({tile1, tile3})) ? empty : D.find({tile1, tile3});
						if (pmap->is_local(tile0123))
							WorldObject_::task(
										me, &PeriodicFourCenterFockBuilder::compute_task, D23, D13,
										std::array<size_t, 4>{{tile0, tile1, tile2, tile3}});
					}
				}
			}
		}

		compute_world.gop.fence();

		// cleanup
		if (compute_J_)
			j_engines_.reset();
		if (compute_K_)
			k_engines_.reset();

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

	void compute_task(Tile D23, Tile D13, std::array<size_t, 4> tile_idx) {

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
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_FOCK_BUILDER_H_
