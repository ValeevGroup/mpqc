#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_

#include "mpqc/chemistry/qc/lcao/integrals/density_fitting/cadf_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_three_center_contraction_builder.h"

#include "mpqc/math/external/tiledarray/array_info.h"

#include <boost/functional/hash.hpp>
#include <unordered_set>

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicCADFKBuilder
    : public madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using DirectTArray = typename Factory::DirectTArray;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

  using WorldObject_ =
      madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>>;
  using PeriodicCADFKBuilder_ = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  PeriodicCADFKBuilder(madness::World &world, Factory &ao_factory,
                       const double force_shape_threshold = 0.0)
      : WorldObject_(world),
        ao_factory_(ao_factory),
        trans_(madness::cblas::Trans),
        notrans_(madness::cblas::NoTrans) {
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    print_detail_ = ao_factory_.print_detail();
    shell_pair_threshold_ = ao_factory_.shell_pair_threshold();
    density_threshold_ = ao_factory_.density_threshold();
    force_shape_threshold_ = force_shape_threshold;
    ExEnv::out0() << "\nforce shape threshold = " << force_shape_threshold_
                  << std::endl;
    target_precision_ = std::numeric_limits<double>::epsilon();
    mpqc::time_point t0, t1;

    // by-cluster orbital basis and df basis
    obs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    dfbs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));
    assert(obs_->nclusters() == dfbs_->nclusters());
    ntiles_per_uc_ = obs_->nclusters();

    dcell_ = ao_factory_.unitcell().dcell();
    R_max_ = ao_factory_.R_max();
    RJ_max_ = ao_factory_.RJ_max();
    RD_max_ = ao_factory_.RD_max();
    R_size_ = ao_factory_.R_size();
    RJ_size_ = ao_factory_.RJ_size();
    RD_size_ = ao_factory_.RD_size();

    const Vector3i ref_latt_range = {0, 0, 0};
    const auto natoms_per_uc = ao_factory_.unitcell().natoms();
    Vector3d zero_shift_base(0.0, 0.0, 0.0);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    // determine max lattice range for product density |μ ρ_Rj)
    t0 = mpqc::fenced_now(world);
    {
      // build a temp basisRJ using user-given RJ_max_
      auto basisRJ =
          shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);

      // compute significant shell pair list
      sig_shellpair_list_ = parallel_compute_shellpair_list(
          *obs_, *basisRJ, shell_pair_threshold_);
      // make a list of significant Rj's as in overlap between μ and ρ_Rj
      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        const auto RJ_3D = direct_3D_idx(RJ, RJ_max_);
        const auto nshells = obs_->flattened_shells().size();
        const auto shell1_min = nshells * RJ;
        const auto shell1_max = shell1_min + nshells;

        auto is_significant = false;
        for (auto shell0 = 0; shell0 != nshells; ++shell0) {
          for (const auto &shell1 : sig_shellpair_list_[shell0]) {
            if (shell1 >= shell1_min && shell1 < shell1_max) {
              is_significant = true;
              RJ_list_.emplace_back(RJ_3D);
              break;
            }
          }
          if (is_significant) break;
        }
      }

      ExEnv::out0() << "\nUser specified RJ_max = " << RJ_max_.transpose()
                    << std::endl;
      // renew RJ_max_, RJ_size_, and basisRJ_
      auto x = 0;
      auto y = 0;
      auto z = 0;
      for (const auto &RJ_3D : RJ_list_) {
        x = std::max(x, RJ_3D(0));
        y = std::max(y, RJ_3D(1));
        z = std::max(z, RJ_3D(2));
      }
      RJ_max_ = Vector3i({x, y, z});
      RJ_size_ = 1 + direct_ord_idx(RJ_max_, RJ_max_);
      basisRJ_ = shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);
      ExEnv::out0() << "Updated RJ_max = " << RJ_max_.transpose() << std::endl;
    }
    t1 = mpqc::fenced_now(world);
    auto t_update_rjmax = mpqc::duration_in_s(t0, t1);

    // compute C(X, μ_0, ρ_Rj)
    t0 = mpqc::fenced_now(world);
    {
      X_dfbs_ = shift_basis_origin(*dfbs_, zero_shift_base, RJ_max_, dcell_);
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*X_dfbs_);
      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

      C_bra_ = lcao::cadf_fitting_coefficients<Tile, Policy>(
          M, *obs_, *basisRJ_, *X_dfbs_, natoms_per_uc, ref_latt_range, RJ_max_,
          RJ_max_);
    }
    t1 = mpqc::fenced_now(world);
    auto t_C_bra = mpqc::duration_in_s(t0, t1);

    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };

    RY_max_ = max_latt_range(R_max_, RJ_max_ + RD_max_);
    Y_dfbs_ = shift_basis_origin(*dfbs_, zero_shift_base, RY_max_, dcell_);

    // compute M(X, Y)
    t0 = mpqc::fenced_now(world);
    {
      RYmRX_max_ = RJ_max_ + RY_max_;
      auto shifted_Y_dfbs =
          shift_basis_origin(*dfbs_, zero_shift_base, RYmRX_max_, dcell_);
      M_ = compute_eri2(world, *dfbs_, *shifted_Y_dfbs);
    }
    t1 = mpqc::fenced_now(world);
    auto t_M = mpqc::duration_in_s(t0, t1);

    // make direct integral E_ket_
    t0 = mpqc::fenced_now(world);
    {
      auto bs_array = utility::make_array_of_refs(*Y_dfbs_, *obs_, *basisRJ_);
      auto bs_vector =
          lcao::gaussian::BasisVector{{*Y_dfbs_, *obs_, *basisRJ_}};

      auto oper_type = libint2::Operator::coulomb;
      const auto screen_thresh = ao_factory_.screen_threshold();
      auto screen_engine =
          make_engine_pool(oper_type, bs_array, libint2::BraKet::xx_xx);
      auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
          lcao::gaussian::create_schwarz_screener(world, screen_engine,
                                                  bs_vector, screen_thresh));
      auto engine =
          make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);

      E_ket_ = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                       std::move(screener));
    }
    t1 = mpqc::fenced_now(world);
    auto t_direct_eri3 = mpqc::duration_in_s(t0, t1);

    // misc:
    // 1. determine tiled ranges and pmap of exchange term
    // 2. determine translationally invariant basis, tiled ranges and pmap for
    // Q(Y, nu, rho)
    t0 = mpqc::fenced_now(world);
    {
      basisR_ = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);

      // make TiledRange and Pmap of Exchange
      result_trange_ = ::mpqc::lcao::gaussian::detail::create_trange(
          lcao::gaussian::BasisVector{{*obs_, *basisR_}});
      auto tvolume = result_trange_.tiles_range().volume();
      result_pmap_ = Policy::default_pmap(world, tvolume);

      // determine lattice ranges of Q(Y, nu, rho) after translation
      RYmR_max_ = RY_max_ + R_max_;
      RJmR_max_ = RJ_max_ + R_max_;

      // make basis of Q(Y, nu, rho)
      Q_bs_Y_ = shift_basis_origin(*dfbs_, zero_shift_base, RYmR_max_, dcell_);
      Q_bs_nu_ = obs_;
      Q_bs_rho_ = shift_basis_origin(*obs_, zero_shift_base, RJmR_max_, dcell_);
      basisRD_ = shift_basis_origin(*obs_, zero_shift_base, RD_max_, dcell_);

      // make TiledRange and Pmap of Q(Y, nu, rho)
      Q_trange_ = ::mpqc::lcao::gaussian::detail::create_trange(
          lcao::gaussian::BasisVector{{*Q_bs_Y_, *Q_bs_nu_, *Q_bs_rho_}});
      auto Q_tvolume = Q_trange_.tiles_range().volume();
      Q_pmap_ = Policy::default_pmap(world, Q_tvolume);
    }
    t1 = mpqc::fenced_now(world);
    auto t_misc = mpqc::duration_in_s(t0, t1);

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K init time decomposition:\n"
                    << "\tupdate RJ_max:       " << t_update_rjmax << " s\n"
                    << "\tC(X, μ_0, ρ_Rj):     " << t_C_bra << " s\n"
                    << "\tM(X, Y):             " << t_M << " s\n"
                    << "\tdirect ERI3:         " << t_direct_eri3 << " s\n"
                    << "\tmisc:                " << t_misc << " s" << std::endl;
    }
  }

  ~PeriodicCADFKBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    return compute_K(D, target_precision);
  }

 private:
  Factory &ao_factory_;
  bool print_detail_;
  double force_shape_threshold_;
  double target_precision_ = std::numeric_limits<double>::epsilon();
  double shell_pair_threshold_;
  double density_threshold_;
  Vector3d dcell_;
  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  Vector3i RY_max_;
  int64_t R_size_;
  int64_t RJ_size_;
  int64_t RD_size_;
  size_t ntiles_per_uc_;

  array_type C_bra_;
  array_type M_;
  DirectTArray E_ket_;

  shellpair_list_t sig_shellpair_list_;
  std::vector<Vector3i> RJ_list_;

  std::shared_ptr<Basis> obs_;
  std::shared_ptr<Basis> dfbs_;
  std::shared_ptr<Basis> X_dfbs_;
  std::shared_ptr<Basis> Y_dfbs_;
  std::shared_ptr<Basis> basisRJ_;
  std::shared_ptr<Basis> basisR_;
  std::shared_ptr<Basis> basisRD_;

  std::shared_ptr<Basis> Q_bs_Y_;
  std::shared_ptr<Basis> Q_bs_nu_;
  std::shared_ptr<Basis> Q_bs_rho_;
  TA::TiledRange Q_trange_;
  std::shared_ptr<TA::Pmap> Q_pmap_;
  Vector3i RYmR_max_;
  Vector3i RJmR_max_;

  TA::TiledRange F_trange_;
  std::shared_ptr<TA::Pmap> F_pmap_;
  Vector3i RYmRX_max_;

  TA::TiledRange result_trange_;
  std::shared_ptr<TA::Pmap> result_pmap_;

  madness::ConcurrentHashMap<std::size_t, Tile> local_contr_tiles_;
  madness::ConcurrentHashMap<std::size_t, Tile> global_contr_tiles_;
  std::atomic<size_t> num_ints_computed_{0};

  const madness::cblas::CBLAS_TRANSPOSE trans_;
  const madness::cblas::CBLAS_TRANSPOSE notrans_;

 private:
  array_type compute_K(const array_type &D, double target_precision) {
    auto &world = ao_factory_.world();

    auto t0_k = mpqc::fenced_now(world);
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    array_type K;

    time_point t0, t1;
    double t_Qket = 0.0;
    double t_eval_Eket = 0.0;
    double t_CM = 0.0;
    double t_F = 0.0;
    double t_permute = 0.0;
    double t_K = 0.0;
    num_ints_computed_ = 0;

    // try new method
    {
      t0 = mpqc::fenced_now(world);
      // compute translational invariant Qket
      // Q(Y, ν_R, ρ_Rj) = C(Y, ν_R, σ_(Rj+Rd)) D(ρ_0, σ_Rd)
      auto Q_ket = compute_Q_ket(C_bra_, D);
      t1 = mpqc::fenced_now(world);
      t_Qket += mpqc::duration_in_s(t0, t1);

      detail::print_size_info(Q_ket, "Q(Y,ν,ρ)");

      // compute F(Y, μ_0, ρ_Rj) = 2 * E(Y, μ_0, ρ_Rj) - C(X, μ_0, ρ_Rj) M(X, Y)
      t0 = mpqc::fenced_now(world);
      array_type F;
      {
        auto t0_eri3 = mpqc::fenced_now(world);
        auto forced_norms =
            force_F_norms(Q_ket.shape().data(), E_ket_.array().shape().data());
        auto trange = E_ket_.array().trange();
        TA::SparseShape<float> forced_shape(world, forced_norms, trange);

        F("Y, mu, rho") = (E_ket_("Y, mu, rho")).set_shape(forced_shape);
        F.truncate();

        F("Y, mu, rho") = 2.0 * F("Y, mu, rho");
        auto t1_eri3 = mpqc::fenced_now(world);
        t_eval_Eket = mpqc::duration_in_s(t0_eri3, t1_eri3);

        auto t0_CM = mpqc::fenced_now(world);
        F_trange_ = F.trange();
        F_pmap_ = F.pmap();
        F("Y, mu, rho") -=
            compute_contr_CM(C_bra_, M_, forced_norms)("Y, mu, rho");
        F.truncate();
        auto t1_CM = mpqc::fenced_now(world);
        t_CM = mpqc::duration_in_s(t0_CM, t1_CM);
      }
      t1 = mpqc::fenced_now(world);
      t_F = mpqc::duration_in_s(t0, t1);

      detail::print_size_info(F, "F(Y,μ,ρ)");

      // permute basis indices
      t0 = mpqc::fenced_now(world);
      Q_ket("Y, rho, nu") = Q_ket("Y, nu, rho");
      F("Y, rho, mu") = F("Y, mu, rho");
      t1 = mpqc::fenced_now(world);
      t_permute = mpqc::duration_in_s(t0, t1);

      // compute K(μ_0, ν_R) = F(Y, ρ_Rj, μ_0) Q(Y, ρ_Rj, ν_R)
      t0 = mpqc::fenced_now(world);
      K = compute_contr_FQ(F, Q_ket);
      t1 = mpqc::fenced_now(world);
      t_K = mpqc::duration_in_s(t0, t1);
    }

    auto t1_k = mpqc::fenced_now(world);
    auto t_tot = mpqc::duration_in_s(t0_k, t1_k);

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K time decomposition:\n"
                    << "\tQ_ket(Y, ν_R, ρ_Rj) :     " << t_Qket << " s\n"
                    << "\tF = 2 * E_ket - C M:      " << t_F << " s\n"
                    << "\t  Eval E_ket(Y, μ_0, ρ):  " << t_eval_Eket << " s\n"
                    << "\t  Contract C M:           " << t_CM << " s\n"
                    << "\tPermute F and Q_ket:      " << t_permute << " s\n"
                    << "\tK = F Q_ket:              " << t_K << " s\n"
                    << "\nTotal K builder time:     " << t_tot << " s"
                    << std::endl;
    }

    return K;
  }

  /*!
   * \brief This computes forced norms of F(Y, μ_0, ρ_Rj) based on the sparsity
   * of Q(Y, ν_R, ρ_Rj).
   *
   * Note that K(μ_0, ν_R) = F(Y, μ_0, ρ_Rj) Q(Y, ν_R, ρ_Rj). For specific Y and
   * ρ_Rj, no need to compute F(Y, μ_0, ρ_Rj) if all ν_R in Q(Y, ν_R, ρ_Rj)
   * gives zero. Translational symmetry of Q is used.
   *
   * \param in_norms norms of Q(Y, ν_R, ρ_Rj)
   * \param out_norms norms of F(Y, μ_0, ρ_Rj). Only correct ranges of F norms
   * are needed.
   * \return forced norms of F
   */
  TA::Tensor<float> force_F_norms(TA::Tensor<float> const &in_norms,
                                  TA::Tensor<float> const &out_norms) {
    const auto ntiles_Y = Y_dfbs_->nclusters();
    const auto ntiles_nu = basisR_->nclusters();
    const auto ntiles_rho = basisRJ_->nclusters();

    using SigPair = std::pair<size_t, size_t>;
    std::unordered_set<SigPair, boost::hash<SigPair>> Y_rho;
    Y_rho.reserve(ntiles_Y * ntiles_rho);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    for (auto tile_Y = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
      const auto RY_ord = tile_Y / ntiles_per_uc_;
      const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
      const auto tile_Y_in_uc = tile_Y % ntiles_per_uc_;
      for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
        const auto RJ_ord = tile_rho / ntiles_per_uc_;
        const auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
        const auto tile_rho_in_uc = tile_rho % ntiles_per_uc_;

        SigPair Y_rho_pair(tile_Y, tile_rho);
        for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
          const auto R_ord = tile_nu / ntiles_per_uc_;
          const auto R_3D = direct_3D_idx(R_ord, R_max_);

          const auto RYmR_3D = RY_3D - R_3D;
          const auto RJmR_3D = RJ_3D - R_3D;

          const auto RYmR_ord = direct_ord_idx(RYmR_3D, RYmR_max_);
          const auto RJmR_ord = direct_ord_idx(RJmR_3D, RJmR_max_);

          const auto shifted_Y = tile_Y_in_uc + RYmR_ord * ntiles_per_uc_;
          const auto shifted_nu = tile_nu % ntiles_per_uc_;
          const auto shifted_rho = tile_rho_in_uc + RJmR_ord * ntiles_per_uc_;

          const auto val = in_norms(shifted_Y, shifted_nu, shifted_rho);
          if (val > force_shape_threshold_) {
            Y_rho.insert(Y_rho_pair);
            break;
          }
        }
      }
    }

    const auto &out_range = out_norms.range();
    TA::Tensor<float> out(out_range, 0.0);

    // F("X, mu, rho")
    for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
      for (auto const &Y_rho_pair : Y_rho) {
        const auto tile_Y = Y_rho_pair.first;
        const auto tile_rho = Y_rho_pair.second;
        out(tile_Y, tile_mu, tile_rho) = std::numeric_limits<float>::max();
      }
    }

    return out;
  }

  /*!
   * \brief This computes ERI2 (\c bs0 | \c bs1)
   * \param world
   * \param bs0
   * \param bs1
   * \return
   */
  array_type compute_eri2(madness::World &world,
                          const lcao::gaussian::Basis &bs0,
                          const lcao::gaussian::Basis &bs1) {
    const auto bs_ref_array = utility::make_array_of_refs(bs0, bs1);
    const auto bs_vector = lcao::gaussian::BasisVector{{bs0, bs1}};
    auto engine = lcao::gaussian::make_engine_pool(
        libint2::Operator::coulomb, bs_ref_array, libint2::BraKet::xs_xs);
    return lcao::gaussian::sparse_integrals(world, engine, bs_vector);
  }

  /*!
   * \brief This computes translationally invariant part of Q(Y_Ry, ν_R, ρ_Rj),
   * i.e. Q(Y(Ry-R), ν_0, ρ_(Rj-R)).
   *
   * Note that Q(Y_Ry, ν_R, ρ_Rj) = C(Y_Ry, ν_R, σ_(Rj+Rd)) D(ρ_0, σ_Rd). Thus
   * Q(Y(Ry-R), ν_0, ρ_(Rj-R)) = C(Y_(Ry-R), ν_0, σ_(Rj+Rd-R)) * D(ρ_0, σ_Rd).
   * Translational symmetry of C is also used.
   *
   * \param C translationally invariant C, i.e. C(Y_(Ry-R), ν_0, σ_(Rj+Rd-R))
   * \param D density matrix
   * \return
   */
  array_type compute_Q_ket(const array_type &C, const array_type &D) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();

    array_type C_repl, D_repl;
    C_repl("Y, nu, sig") = C("Y, nu, sig");
    D_repl("rho, sig") = D("rho, sig");
    C_repl.make_replicated();
    world.gop.fence();

    D_repl.make_replicated();
    // must wait till all replicating is finished
    world.gop.fence();

    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    // # of tiles per basis
    const auto ntiles_Y = Q_bs_Y_->nclusters();
    const auto ntiles_nu = Q_bs_nu_->nclusters();
    const auto ntiles_rho = Q_bs_rho_->nclusters();
    const auto ntiles_sig = basisRD_->nclusters();

    const auto Dtile_norms = D_repl.shape().data();
    for (auto tile_Y = 0ul, task = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
      const auto RYmR_ord = tile_Y / ntiles_per_uc_;
      const auto RYmR_3D = direct_3D_idx(RYmR_ord, RYmR_max_);
      if (!is_in_lattice_range(RYmR_3D, RJ_max_)) continue;

      const auto RY_ord_in_C = direct_ord_idx(RYmR_3D, RJ_max_);
      const auto Y_in_C =
          tile_Y % ntiles_per_uc_ + RY_ord_in_C * ntiles_per_uc_;

      for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
        for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
          const auto RJmR_ord = tile_rho / ntiles_per_uc_;
          const auto RJmR_3D = direct_3D_idx(RJmR_ord, RJmR_max_);
          const auto rho_in_D = tile_rho % ntiles_per_uc_;

          for (auto tile_sig = 0ul; tile_sig != ntiles_sig;
               ++tile_sig, ++task) {
            if (task % nproc == me) {
              const auto RD_ord = tile_sig / ntiles_per_uc_;
              const auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
              const auto RJmRpRD_3D = RJmR_3D + RD_3D;
              if (!(RYmR_3D == Vector3i({0, 0, 0})) && !(RYmR_3D == RJmRpRD_3D))
                continue;

              if (!is_in_lattice_range(RJmRpRD_3D, RJ_max_)) continue;

              if (std::find(RJ_list_.begin(), RJ_list_.end(), RJmRpRD_3D) ==
                  RJ_list_.end())
                continue;

              const auto RJmRpRD_ord = direct_ord_idx(RJmRpRD_3D, RJ_max_);
              const auto sig_in_C =
                  tile_sig % ntiles_per_uc_ + RJmRpRD_ord * ntiles_per_uc_;

              // get future of tiles C and D
              std::array<size_t, 3> idx_C = {{Y_in_C, tile_nu, sig_in_C}};
              std::array<size_t, 2> idx_D = {{rho_in_D, tile_sig}};
              if (C_repl.is_zero(idx_C) || D_repl.is_zero(idx_D)) continue;
              if (Dtile_norms(rho_in_D, tile_sig) <= density_threshold_)
                continue;

              auto C_tile = C_repl.find(idx_C);
              auto D_tile = D_repl.find(idx_D);
              WorldObject_::task(
                  me, &PeriodicCADFKBuilder_::compute_Q_ket_task, C_tile,
                  D_tile, std::array<size_t, 3>{{tile_Y, tile_nu, tile_rho}});
            }
          }
        }
      }
    }

    world.gop.fence();

    // collect local tiles
    for (const auto &local_tile : local_contr_tiles_) {
      const auto tile_ord = local_tile.first;
      const auto proc = Q_pmap_->owner(tile_ord);
      WorldObject_::task(proc, &PeriodicCADFKBuilder_::accumulate_global_task,
                         local_tile.second, tile_ord);
    }
    local_contr_tiles_.clear();
    world.gop.fence();

    typename Policy::shape_type shape;
    // compute the shape, if sparse
    if (!decltype(shape)::is_dense()) {
      // extract local contribution to the shape of G, construct global shape
      std::vector<std::pair<std::array<size_t, 3>, double>> global_tile_norms;
      const auto i_stride = ntiles_nu * ntiles_rho;
      const auto j_stride = ntiles_rho;
      for (const auto &global_tile : global_contr_tiles_) {
        const auto tile_ord = global_tile.first;
        const auto i = tile_ord / i_stride;
        const auto jk = tile_ord % i_stride;
        const auto j = jk / j_stride;
        const auto k = jk % j_stride;
        const auto norm = global_tile.second.norm();
        global_tile_norms.push_back(
            std::make_pair(std::array<size_t, 3>{{i, j, k}}, norm));
      }
      shape = decltype(shape)(world, global_tile_norms, Q_trange_);
    }

    array_type result(world, Q_trange_, shape, Q_pmap_);
    for (const auto &global_tile : global_contr_tiles_) {
      if (!result.shape().is_zero(global_tile.first))
        result.set(global_tile.first, global_tile.second);
    }
    result.fill_local(0.0, true);
    global_contr_tiles_.clear();

    return result;
  }

  /*!
   * \brief This computes contraction between a C tile and a D tile.
   * \param C_tile tile of C
   * \param D_tile tile of D
   * \param tile_idx result index
   */
  void compute_Q_ket_task(Tile C_tile, Tile D_tile,
                          std::array<size_t, 3> tile_idx) {
    const auto ext_C = C_tile.range().extent();
    const auto ext_D = D_tile.range().extent();
    assert(ext_C[2] == ext_D[1]);

    if (C_tile.norm() * D_tile.norm() >= target_precision_) {
      const auto tile_Y = tile_idx[0];
      const auto tile_nu = tile_idx[1];
      const auto tile_rho = tile_idx[2];

      const auto &rng_Y = Q_trange_.dim(0).tile(tile_Y);
      const auto &rng_nu = Q_trange_.dim(1).tile(tile_nu);
      const auto &rng_rho = Q_trange_.dim(2).tile(tile_rho);

      const auto rng_Y_size = rng_Y.second - rng_Y.first;
      const auto rng_nu_size = rng_nu.second - rng_nu.first;
      const auto rng_rho_size = rng_rho.second - rng_rho.first;
      assert(rng_Y_size == ext_C[0]);
      assert(rng_nu_size == ext_C[1]);
      assert(rng_rho_size == ext_D[0]);

      const auto result_rng = TA::Range({rng_Y, rng_nu, rng_rho});
      Tile result_tile(result_rng, 0.0);

      TA::math::GemmHelper gh(notrans_, trans_, result_tile.range().rank(),
                              C_tile.range().rank(), D_tile.range().rank());

      int m, k, n;
      gh.compute_matrix_sizes(m, n, k, C_tile.range(), D_tile.range());
      const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
      const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

      // Notice that we reversed notrans and trans. This is because Lapack
      // expects col major matrices.
      TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, C_tile.data(),
                     lda, D_tile.data(), ldb, 0.0, result_tile.data(), n);

      const auto ntiles_nu = Q_bs_nu_->nclusters();
      const auto ntiles_rho = Q_bs_rho_->nclusters();
      const auto ord =
          tile_Y * ntiles_nu * ntiles_rho + tile_nu * ntiles_rho + tile_rho;

      PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
    }
  }

  /*!
   * \brief This computes K(μ_0, ν_R) = F(Y, ρ_Rj, μ_0) Q(Y, ρ_Rj, ν_R).
   * Translational symmetry of Q is used.
   *
   * \param F array F
   * \param Q translationally invariant Q, i.e. Q(Y(Ry-R), ρ_(Rj-R), ν_0)
   * \return
   */
  array_type compute_contr_FQ(const array_type &F, const array_type &Q) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();

    // # of tiles per basis
    const auto ntiles_Y = Y_dfbs_->nclusters();
    const auto ntiles_rho = basisRJ_->nclusters();
    const auto ntiles_nu = basisR_->nclusters();
    const auto ntiles_mu = obs_->nclusters();

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    for (auto tile_Y = 0ul, task = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
      const auto RY_ord = tile_Y / ntiles_per_uc_;
      const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
      const auto tile_Y_in_uc = tile_Y % ntiles_per_uc_;

      for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
        const auto RJ_ord = tile_rho / ntiles_per_uc_;
        const auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
        const auto tile_rho_in_uc = tile_rho % ntiles_per_uc_;

        for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
          const auto R_ord = tile_nu / ntiles_per_uc_;
          const auto R_3D = direct_3D_idx(R_ord, R_max_);

          const auto RYmR_3D = RY_3D - R_3D;
          const auto RJmR_3D = RJ_3D - R_3D;

          const auto RYmR_ord = direct_ord_idx(RYmR_3D, RYmR_max_);
          const auto RJmR_ord = direct_ord_idx(RJmR_3D, RJmR_max_);

          const auto shifted_Y = tile_Y_in_uc + RYmR_ord * ntiles_per_uc_;
          const auto shifted_rho = tile_rho_in_uc + RJmR_ord * ntiles_per_uc_;
          const auto shifted_nu = tile_nu % ntiles_per_uc_;

          std::array<size_t, 3> idx_Q = {{shifted_Y, shifted_rho, shifted_nu}};

          for (auto tile_mu = 0ul; tile_mu != ntiles_mu; ++tile_mu, ++task) {
            if (task % nproc == me) {
              std::array<size_t, 3> idx_F = {{tile_Y, tile_rho, tile_mu}};
              if (F.is_zero(idx_F) || Q.is_zero(idx_Q)) continue;

              auto F_tile = F.find(idx_F);
              auto Q_tile = Q.find(idx_Q);
              WorldObject_::task(
                  me, &PeriodicCADFKBuilder_::compute_contr_FQ_task, F_tile,
                  Q_tile, std::array<size_t, 2>{{tile_mu, tile_nu}});
            }
          }
        }
      }
    }

    world.gop.fence();

    // collect local tiles
    for (const auto &local_tile : local_contr_tiles_) {
      const auto tile_ord = local_tile.first;
      const auto proc = result_pmap_->owner(tile_ord);
      WorldObject_::task(proc, &PeriodicCADFKBuilder_::accumulate_global_task,
                         local_tile.second, tile_ord);
    }
    local_contr_tiles_.clear();
    world.gop.fence();

    typename Policy::shape_type shape;
    // compute the shape, if sparse
    if (!decltype(shape)::is_dense()) {
      // extract local contribution to the shape of G, construct global shape
      std::vector<std::pair<std::array<size_t, 2>, double>> global_tile_norms;
      for (const auto &global_tile : global_contr_tiles_) {
        const auto tile_ord = global_tile.first;
        const auto i = tile_ord / ntiles_nu;
        const auto j = tile_ord % ntiles_nu;
        const auto norm = global_tile.second.norm();
        global_tile_norms.push_back(
            std::make_pair(std::array<size_t, 2>{{i, j}}, norm));
      }
      shape = decltype(shape)(world, global_tile_norms, result_trange_);
    }

    array_type result(world, result_trange_, shape, result_pmap_);
    for (const auto &global_tile : global_contr_tiles_) {
      if (!result.shape().is_zero(global_tile.first))
        result.set(global_tile.first, global_tile.second);
    }
    result.fill_local(0.0, true);
    global_contr_tiles_.clear();

    return result;
  }

  /*!
   * \brief This computes contraction between a F tile and a Q tile.
   * \param F_tile tile of F
   * \param Q_tile tile of Q
   * \param tile_idx result index
   */
  void compute_contr_FQ_task(Tile F_tile, Tile Q_tile,
                             std::array<size_t, 2> tile_idx) {
    const auto ext_F = F_tile.range().extent();
    const auto ext_Q = Q_tile.range().extent();
    assert(ext_F[0] == ext_Q[0]);
    assert(ext_F[1] == ext_Q[1]);

    if (F_tile.norm() * Q_tile.norm() >= target_precision_) {
      const auto tile_mu = tile_idx[0];
      const auto tile_nu = tile_idx[1];

      const auto &rng_mu = result_trange_.dim(0).tile(tile_mu);
      const auto &rng_nu = result_trange_.dim(1).tile(tile_nu);

      const auto rng_mu_rng = rng_mu.second - rng_mu.first;
      const auto rng_nu_rng = rng_nu.second - rng_nu.first;

      assert(rng_mu_rng == ext_F[2]);
      assert(rng_nu_rng == ext_Q[2]);

      const auto result_rng = TA::Range({rng_mu, rng_nu});
      Tile result_tile(result_rng, 0.0);

      TA::math::GemmHelper gh(trans_, notrans_, result_tile.range().rank(),
                              F_tile.range().rank(), Q_tile.range().rank());

      int m, k, n;
      gh.compute_matrix_sizes(m, n, k, F_tile.range(), Q_tile.range());
      const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
      const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

      // Notice that we reversed notrans and trans. This is because Lapack
      // expects col major matrices.
      TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, F_tile.data(),
                     lda, Q_tile.data(), ldb, 0.0, result_tile.data(), n);

      const auto ntiles_nu = basisR_->nclusters();
      const auto ord = tile_mu * ntiles_nu + tile_nu;

      PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
    }
  }

  /*!
   * \brief This computes contraction C(X, μ_0, ρ_Rj) M(X_Rx, Y_Ry).
   * Translational symmetry of M is used.
   *
   * \param C array C
   * \param M translationally invariant M, i.e. M(X_0, Y_(Ry-Rx))
   * \param force_norms forced norms of the result
   * \return
   */
  array_type compute_contr_CM(const array_type &C, const array_type &M,
                              const TA::Tensor<float> &force_norms) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();

    array_type C_repl, M_repl;
    C_repl("X, mu, rho") = C("X, mu, rho");
    M_repl("X, Y") = M("X, Y");
    C_repl.make_replicated();
    world.gop.fence();
    M_repl.make_replicated();
    world.gop.fence();  // must wait till all replicating is finished

    // # of tiles per basis
    const auto ntiles_X = X_dfbs_->nclusters();
    const auto ntiles_Y = Y_dfbs_->nclusters();
    const auto ntiles_rho = basisRJ_->nclusters();
    const auto ntiles_mu = obs_->nclusters();

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    for (auto tile_X = 0ul, task = 0ul; tile_X != ntiles_X; ++tile_X) {
      const auto RX_ord = tile_X / ntiles_per_uc_;
      const auto RX_3D = direct_3D_idx(RX_ord, RJ_max_);
      const auto X_in_uc = tile_X % ntiles_per_uc_;

      for (auto tile_Y = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
        const auto RY_ord = tile_Y / ntiles_per_uc_;
        const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);

        const auto RYmRX_3D = RY_3D - RX_3D;
        const auto RYmRX_ord = direct_ord_idx(RYmRX_3D, RYmRX_max_);
        const auto Y_in_M =
            tile_Y % ntiles_per_uc_ + RYmRX_ord * ntiles_per_uc_;

        std::array<size_t, 2> idx_M = {{X_in_uc, Y_in_M}};
        if (M_repl.is_zero(idx_M)) continue;
        auto M_tile = M_repl.find(idx_M);

        for (auto tile_mu = 0ul; tile_mu != ntiles_mu; ++tile_mu) {
          for (auto tile_rho = 0ul; tile_rho != ntiles_rho;
               ++tile_rho, ++task) {
            if (task % nproc == me) {
              const auto RJ_ord = tile_rho / ntiles_per_uc_;
              const auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
              if (std::find(RJ_list_.begin(), RJ_list_.end(), RJ_3D) ==
                  RJ_list_.end())
                continue;

              std::array<size_t, 3> idx_C = {{tile_X, tile_mu, tile_rho}};
              if (C_repl.is_zero(idx_C) ||
                  force_norms(tile_Y, tile_mu, tile_rho) <
                      force_shape_threshold_)
                continue;

              auto C_tile = C_repl.find(idx_C);

              WorldObject_::task(
                  me, &PeriodicCADFKBuilder_::compute_contr_CM_task, C_tile,
                  M_tile, std::array<size_t, 3>{{tile_Y, tile_mu, tile_rho}});
            }
          }
        }
      }
    }

    world.gop.fence();

    // collect local tiles
    for (const auto &local_tile : local_contr_tiles_) {
      const auto tile_ord = local_tile.first;
      const auto proc = F_pmap_->owner(tile_ord);
      WorldObject_::task(proc, &PeriodicCADFKBuilder_::accumulate_global_task,
                         local_tile.second, tile_ord);
    }
    local_contr_tiles_.clear();
    world.gop.fence();

    typename Policy::shape_type shape;
    // compute the shape, if sparse
    if (!decltype(shape)::is_dense()) {
      // extract local contribution to the shape of G, construct global shape
      std::vector<std::pair<std::array<size_t, 3>, double>> global_tile_norms;
      const auto i_stride = ntiles_mu * ntiles_rho;
      const auto j_stride = ntiles_rho;
      for (const auto &global_tile : global_contr_tiles_) {
        const auto tile_ord = global_tile.first;
        const auto i = tile_ord / i_stride;
        const auto jk = tile_ord % i_stride;
        const auto j = jk / j_stride;
        const auto k = jk % j_stride;
        const auto norm = global_tile.second.norm();
        global_tile_norms.push_back(
            std::make_pair(std::array<size_t, 3>{{i, j, k}}, norm));
      }
      shape = decltype(shape)(world, global_tile_norms, F_trange_);
    }

    array_type result(world, F_trange_, shape, F_pmap_);
    for (const auto &global_tile : global_contr_tiles_) {
      if (!result.shape().is_zero(global_tile.first))
        result.set(global_tile.first, global_tile.second);
    }
    result.fill_local(0.0, true);
    global_contr_tiles_.clear();

    return result;
  }

  /*!
   * \brief This computes contraction between a C tile and a M tile.
   * \param C_tile tile of C
   * \param M_tile tile of M
   * \param tile_idx result index
   */
  void compute_contr_CM_task(Tile C_tile, Tile M_tile,
                             std::array<size_t, 3> tile_idx) {
    const auto ext_C = C_tile.range().extent();
    const auto ext_M = M_tile.range().extent();
    assert(ext_M[0] == ext_C[0]);

    if (C_tile.norm() * M_tile.norm() >= target_precision_) {
      const auto tile_Y = tile_idx[0];
      const auto tile_mu = tile_idx[1];
      const auto tile_rho = tile_idx[2];

      const auto &rng_Y = F_trange_.dim(0).tile(tile_Y);
      const auto &rng_mu = F_trange_.dim(1).tile(tile_mu);
      const auto &rng_rho = F_trange_.dim(2).tile(tile_rho);

      const auto rng_Y_size = rng_Y.second - rng_Y.first;
      const auto rng_mu_size = rng_mu.second - rng_mu.first;
      const auto rng_rho_size = rng_rho.second - rng_rho.first;

      assert(rng_Y_size == ext_M[1]);
      assert(rng_mu_size == ext_C[1]);
      assert(rng_rho_size == ext_C[2]);

      const auto result_rng = TA::Range({rng_Y, rng_mu, rng_rho});
      Tile result_tile(result_rng, 0.0);

      TA::math::GemmHelper gh(trans_, notrans_, result_tile.range().rank(),
                              M_tile.range().rank(), C_tile.range().rank());

      int m, k, n;
      gh.compute_matrix_sizes(m, n, k, M_tile.range(), C_tile.range());
      const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
      const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

      // Notice that we reversed notrans and trans. This is because Lapack
      // expects col major matrices.
      TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, M_tile.data(),
                     lda, C_tile.data(), ldb, 0.0, result_tile.data(), n);

      const auto ntiles_mu = obs_->nclusters();
      const auto ntiles_rho = basisRJ_->nclusters();
      const auto ord =
          tile_Y * ntiles_mu * ntiles_rho + tile_mu * ntiles_rho + tile_rho;

      PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
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

  void dump_shape_02_1(array_type const &Q, std::string const &name,
                       int rj = 0) {
    std::cout << "Dumping shape for 3-d tensor (02-1) to file " << name << "_"
              << rj << ".csv ...\n";
    auto norms = Q.shape().data();
    auto range = norms.range();
    auto ext = range.extent_data();

    Eigen::MatrixXd M(ext[0] * ext[2], ext[1]);
    M.setZero();

    for (auto X = 0; X < ext[0]; ++X) {
      for (auto mu = 0; mu < ext[1]; ++mu) {
        for (auto sig = 0; sig < ext[2]; ++sig) {
          M(X * ext[2] + sig, mu) = norms(X, mu, sig);
        }
      }
    }

    std::ofstream outfile(name + "_" + std::to_string(rj) + ".csv");
    auto ncols = M.cols();
    auto nrows = M.rows();
    for (auto i = 0; i < nrows; ++i) {
      for (auto j = 0; j < ncols; ++j) {
        outfile << M(i, j);
        if (j != ncols - 1) {
          outfile << ", ";
        } else {
          outfile << "\n";
        }
      }
    }
  }

  void dump_shape_01_2(array_type const &C, std::string const &name,
                       int rj = 0) {
    std::cout << "Dumping shape for 3-d tensor (01-2) to file " << name << "_"
              << rj << ".csv ...\n";
    auto norms = C.shape().data();
    auto range = norms.range();
    auto ext = range.extent_data();

    Eigen::MatrixXd M(ext[0] * ext[1], ext[2]);
    M.setZero();

    for (auto X = 0; X < ext[0]; ++X) {
      for (auto mu = 0; mu < ext[1]; ++mu) {
        for (auto rho = 0; rho < ext[2]; ++rho) {
          M(X * ext[1] + mu, rho) = norms(X, mu, rho);
        }
      }
    }

    std::ofstream outfile(name + "_" + std::to_string(rj) + ".csv");
    auto ncols = M.cols();
    auto nrows = M.rows();
    for (auto i = 0; i < nrows; ++i) {
      for (auto j = 0; j < ncols; ++j) {
        outfile << M(i, j);
        if (j != ncols - 1) {
          outfile << ", ";
        } else {
          outfile << "\n";
        }
      }
    }
  }

  void dump_shape_0_12(array_type const &C, std::string const &name,
                       int rj = 0) {
    std::cout << "Dumping shape for 3-d tensor (0-12) to file " << name << "_"
              << rj << ".csv ...\n";
    auto norms = C.shape().data();
    auto range = norms.range();
    auto ext = range.extent_data();

    Eigen::MatrixXd M(ext[0], ext[1] * ext[2]);
    M.setZero();

    for (auto X = 0; X < ext[0]; ++X) {
      for (auto mu = 0; mu < ext[1]; ++mu) {
        for (auto rho = 0; rho < ext[2]; ++rho) {
          M(X, mu * ext[2] + rho) = norms(X, mu, rho);
        }
      }
    }

    std::ofstream outfile(name + "_" + std::to_string(rj) + ".csv");
    auto ncols = M.cols();
    auto nrows = M.rows();
    for (auto i = 0; i < nrows; ++i) {
      for (auto j = 0; j < ncols; ++j) {
        outfile << M(i, j);
        if (j != ncols - 1) {
          outfile << ", ";
        } else {
          outfile << "\n";
        }
      }
    }
  }

  void dump_shape_d2(array_type const &M, std::string const &name, int rj = 0) {
    std::cout << "Dumping matrix to file " << name << "_" << rj << ".csv ...\n";
    auto norms = M.shape().data();
    auto range = norms.range();
    auto ext = range.extent_data();

    Eigen::MatrixXd M_eig(ext[0], ext[1]);
    M_eig.setZero();

    for (auto r = 0; r < ext[0]; ++r) {
      for (auto c = 0; c < ext[1]; ++c) {
        M_eig(r, c) = norms(r, c);
      }
    }

    std::ofstream outfile(name + "_" + std::to_string(rj) + ".csv");
    auto ncols = M_eig.cols();
    auto nrows = M_eig.rows();
    for (auto i = 0; i < nrows; ++i) {
      for (auto j = 0; j < ncols; ++j) {
        outfile << M_eig(i, j);
        if (j != ncols - 1) {
          outfile << ", ";
        } else {
          outfile << "\n";
        }
      }
    }
  }

  void dump_matrix(array_type const &M, std::string const &name, int rj = 0) {
    std::cout << "Dumping matrix to file " << name << "_" << rj << ".csv ...\n";

    RowMatrixXd M_eig = array_ops::array_to_eigen(M);

    std::ofstream outfile(name + "_" + std::to_string(rj) + ".csv");
    auto ncols = M_eig.cols();
    auto nrows = M_eig.rows();
    for (auto i = 0; i < nrows; ++i) {
      for (auto j = 0; j < ncols; ++j) {
        outfile << M_eig(i, j);
        if (j != ncols - 1) {
          outfile << ", ";
        } else {
          outfile << "\n";
        }
      }
    }
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

    auto &world = ao_factory_.world();
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
   * \brief This determines if a unit cell is included by the give lattice range
   * \param in_idx 3D index of a unit cell
   * \param range lattice range
   * \param center center of \range
   * \return
   */
  bool is_in_lattice_range(Vector3i const &in_idx, Vector3i const &range,
                           Vector3i const &center = {0, 0, 0}) {
    if (in_idx(0) <= center(0) + range(0) &&
        in_idx(0) >= center(0) - range(0) &&
        in_idx(1) <= center(1) + range(1) &&
        in_idx(1) >= center(1) - range(1) &&
        in_idx(2) <= center(2) + range(2) && in_idx(2) >= center(2) - range(2))
      return true;
    else
      return false;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
