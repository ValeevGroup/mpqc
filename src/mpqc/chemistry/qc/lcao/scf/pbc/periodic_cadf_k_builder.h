#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_

#include "mpqc/chemistry/qc/lcao/basis/util.h"
#include "mpqc/chemistry/qc/lcao/integrals/density_fitting/cadf_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

#include "mpqc/math/external/tiledarray/array_info.h"

#include <boost/functional/hash.hpp>
#include <unordered_set>

namespace mpqc {
namespace lcao {
namespace scf {

/// PeriodicCADFKBuilder computes the exchange term in periodic HF using
/// the concentric atomic density fitting (CADF) approximation.
///
/// K(μ_0, ρ_Rρ) = (μ_0 ν_Rν | ρ_Rρ σ_(Rρ+Rσ)) D(ν_Rν, σ_(Rρ+Rσ))
/// The left/right product density in the 2-body 4-center ERI is approximated
/// by CADF, i.e. |μ_0 ν_Rν) = Sum_X C(X, μ_0, ν_Rν) |X) where X is on the
/// center of either μ0 or νR_ν. Dunlap's robust formula is used.
template <typename Tile, typename Policy, typename Factory>
class PeriodicCADFKBuilder
    : public madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using DirectTArray = typename Factory::DirectTArray;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using shellpair_list_t = std::vector<std::vector<size_t>>;

  using WorldObject_ =
      madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>>;
  using PeriodicCADFKBuilder_ = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  /*!
   * \brief This constructs PeriodicCADFKBuilder using a PeriodicAOFactory
   * object.
   * \param world MADNESS world object
   * \param ao_factory PeriodicAOFactory object
   * \param force_shape_threshold threshold used to construct the shape of
   * F(Υ, μ, ν) using the shape of Q(Y, ρ, ν).
   */
  PeriodicCADFKBuilder(Factory &ao_factory)
      : WorldObject_(ao_factory.world()), ao_factory_(ao_factory) {
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    print_detail_ = ao_factory_.print_detail();
    screen_threshold_ = ao_factory_.screen_threshold();
    density_threshold_ = ao_factory_.density_threshold();
    force_hermiticity_ = ao_factory_.force_hermiticity();

    force_shape_threshold_ = ao_factory_.keyval()->template value<double>(
        "force_shape_threshold", std::numeric_limits<double>::epsilon());
    cadf_sparse_threshold_ = ao_factory_.keyval()->template value<double>(
        "cadf_sparse_threshold", std::numeric_limits<double>::epsilon());
    cadf_contr_threshold_ = ao_factory_.keyval()->template value<double>(
        "cadf_contraction_threshold", std::numeric_limits<double>::epsilon());
    ExEnv::out0() << "\nforce shape threshold = " << force_shape_threshold_
                  << "\nCADF sparse threshold = " << cadf_sparse_threshold_
                  << "\nCADF contraction threshold = " << cadf_contr_threshold_
                  << std::endl;

    // by-cluster orbital basis and df basis
    obs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    dfbs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));
    assert(obs_->nclusters() == dfbs_->nclusters());

    ntiles_per_uc_ = obs_->nclusters();
    natoms_per_uc_ = ao_factory_.unitcell().natoms();
    dcell_ = ao_factory_.unitcell().dcell();
    R_max_ = ao_factory_.R_max();
    RJ_max_ = ao_factory_.RJ_max();
    RD_max_ = ao_factory_.RD_max();
    R_size_ = ao_factory_.R_size();

    init();
  }

  /*!
   * \brief This constructs PeriodicCADFKBuilder using specific basis sets,
   * lattice params, and thresholds.
   * \param world MADNESS world object
   * \param ao_factory PeriodicAOFactory object
   * \param obs orbital basis set
   * \param dfbs density fitting basis set
   * \param dcell unit cell parameters
   * \param R_max range of expansion of Bloch Gaussians in AO Gaussians
   * \param RJ_max range of Coulomb operation
   * \param RD_max range of density representation
   * \param shell_pair_threshold threshold for screeing non-negligible shell
   * pairs
   * \param screen_threshold threshold for schwarz screening.
   * \param density_threshold threshold for screening density blocks
   * \param target_precision controls libint engine precision and the sparsity
   * of C.D, C.M and F.Q products
   * \param print_detail print more details if true
   * \param force_shape_threshold threshold used to construct the shape of
   * F(Υ, μ, ν) using the shape of Q(Y, ρ, ν).
   */
  PeriodicCADFKBuilder(
      madness::World &world, Factory &ao_factory, std::shared_ptr<Basis> obs,
      std::shared_ptr<Basis> dfbs, const Vector3d &dcell, const Vector3i &R_max,
      const Vector3i &RJ_max, const Vector3i &RD_max,
      const double screen_threshold = 1.0e-20,
      const double density_threshold = Policy::shape_type::threshold(),
      const double cadf_contr_threshold = std::numeric_limits<double>::epsilon(),
      const bool print_detail = false, const double force_shape_threshold = 0.0)
      : WorldObject_(world),
        ao_factory_(ao_factory),
        obs_(obs),
        dfbs_(dfbs),
        dcell_(dcell),
        R_max_(R_max),
        RJ_max_(RJ_max),
        RD_max_(RD_max),
        print_detail_(print_detail),
        screen_threshold_(screen_threshold),
        density_threshold_(density_threshold),
        force_shape_threshold_(force_shape_threshold),
        cadf_contr_threshold_(cadf_contr_threshold) {
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    ExEnv::out0() << "\nforce shape threshold = " << force_shape_threshold_
                  << std::endl;

    MPQC_ASSERT(obs_->nclusters() == dfbs_->nclusters());
    ntiles_per_uc_ = obs_->nclusters();
    natoms_per_uc_ = ao_factory_.unitcell().natoms();

    using ::mpqc::detail::direct_ord_idx;
    R_size_ = 1 + direct_ord_idx(R_max_, R_max_);

    // by-cluster orbital basis and df basis
    assert(obs_->nclusters() == dfbs_->nclusters());

    init();
  }

  ~PeriodicCADFKBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    return compute_K(D, target_precision);
  }

  Vector3i K_lattice_range() { return Rrho_max_; }

 private:
  // the following private member variables are set up by ctor
  Factory &ao_factory_;
  std::shared_ptr<Basis> obs_;
  std::shared_ptr<Basis> dfbs_;
  Vector3d dcell_;
  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  size_t ntiles_per_uc_;
  size_t natoms_per_uc_;
  bool print_detail_;
  double screen_threshold_;
  double density_threshold_;

  double cadf_sparse_threshold_;
  double force_shape_threshold_;
  double cadf_contr_threshold_ = std::numeric_limits<double>::epsilon();
  bool force_hermiticity_;

  int64_t R_size_;
  std::shared_ptr<Basis> basisR_;
  Vector3i truncated_RD_max_;
  Vector3i Rrho_max_;
  Vector3i RX_max_;
  Vector3i RY_max_;
  Vector3i R1m2_max_;
  Vector3i RYmRX_max_;
  int64_t Rrho_size_;
  int64_t RY_size_;
  int64_t R1m2_size_;
  int64_t ref_uc_ord_;

  array_type C_;
  array_type M_;
  DirectTArray eri3_;

  TA::TiledRange Q_trange_;
  TA::TiledRange F_trange_;
  TA::TiledRange result_trange_;
  std::shared_ptr<TA::Pmap> Q_pmap_;
  std::shared_ptr<TA::Pmap> F_pmap_;
  std::shared_ptr<TA::Pmap> result_pmap_;

  madness::ConcurrentHashMap<std::size_t, Tile> local_contr_tiles_;
  madness::ConcurrentHashMap<std::size_t, Tile> global_contr_tiles_;

 private:
  void init() {
    auto &world = this->get_world();

    mpqc::time_point t0, t1;

    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

    // ordinal # of the reference unit cell in R_max_ lattice range
    ref_uc_ord_ = (R_size_ - 1) / 2;
    // make compound basis for ν in product density |μ_0 ν_Rν)
    basisR_ = shift_basis_origin(*obs_, Vector3d::Zero(), R_max_, dcell_);

    // compute C(X_Rx, μ_0, ν_Rν)
    t0 = mpqc::fenced_now(world);
    {
      // X is on the center of either μ_0 or ν_Rν. Thus X_Rx should have the
      // same lattice range as ν_Rν.
      RX_max_ = R_max_;
      auto X_dfbs =
          shift_basis_origin(*dfbs_, Vector3d::Zero(), RX_max_, dcell_);
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*X_dfbs);
      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);
      const Vector3i ref_lattice_range = {0, 0, 0};

      C_ = lcao::cadf_fitting_coefficients<Tile, Policy>(
          M, *obs_, *basisR_, *X_dfbs, natoms_per_uc_, ref_lattice_range,
          R_max_, RX_max_);
    }
    t1 = mpqc::fenced_now(world);
    auto t_C = mpqc::duration_in_s(t0, t1);

    mpqc::detail::print_size_info(C_, "C(X,μ,ν)");

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K init time decomposition:\n"
                    << "\tC(X_Rx, μ_0, ν_Rν):  " << t_C << " s\n"
                    << std::endl;
    }

    truncated_RD_max_ = RD_max_;  // make them equal for initialization
    // do not forget to update RD-dependent variables
    update_RD_dependent_variables(RD_max_);
  }

  array_type compute_K(const array_type &D, double target_precision) {
    auto &world = this->get_world();

    auto t0_k = mpqc::fenced_now(world);

    array_type K;

    time_point t0, t1;
    double t_Q = 0.0;
    double t_eval_eri3 = 0.0;
    double t_CM = 0.0;
    double t_F = 0.0;
    double t_permute = 0.0;
    double t_K = 0.0;
    double t_symm = 0.0;

    // Update lattice range of density representation
    ExEnv::out0() << "\nTruncating lattice range of density representation\n";
    Vector3i old_RD_max = truncated_RD_max_;

    using ::mpqc::pbc::detail::truncate_lattice_range;
    truncated_RD_max_ = truncate_lattice_range(D, RD_max_, density_threshold_);
    // Update RD-dependent variables if RD_max is changed
    if (truncated_RD_max_ != old_RD_max) {
      ExEnv::out0() << "\nLattice range of density representation is changed. "
                       "Update RD-dependent variables!"
                    << std::endl;
      update_RD_dependent_variables(truncated_RD_max_);
    } else {
      ExEnv::out0() << "\nLattice range of density representation is not "
                       "changed. No need to update RD-dependent variables!"
                    << std::endl;
    }

    // compute exchange
    {
      t0 = mpqc::fenced_now(world);
      // compute translational invariant Q
      // Q(Y, ν_R, ρ_Rj) = C(Y, ν_R, σ_(Rj+Rd)) D(ρ_0, σ_Rd)
      auto Q = compute_Q(C_, D);
      t1 = mpqc::fenced_now(world);
      t_Q += mpqc::duration_in_s(t0, t1);

      mpqc::detail::print_size_info(Q, "Q(Y,ν,ρ)");

      // compute F(Y, μ_0, ρ_Rj) = 2 * E(Y, μ_0, ρ_Rj) - C(X, μ_0, ρ_Rj) M(X, Y)
      t0 = mpqc::fenced_now(world);
      array_type F;
      {
        auto t0_eri3 = mpqc::fenced_now(world);
        auto forced_norms =
            force_F_norms(Q.shape().data(), eri3_.array().shape().data());
        auto trange = eri3_.array().trange();
        TA::SparseShape<float> forced_shape(world, forced_norms, trange);

        F("Y, mu, nu") = (eri3_("Y, mu, nu")).set_shape(forced_shape);
        F.truncate();

        F("Y, mu, nu") = 2.0 * F("Y, mu, nu");
        auto t1_eri3 = mpqc::fenced_now(world);
        t_eval_eri3 = mpqc::duration_in_s(t0_eri3, t1_eri3);

        auto t0_CM = mpqc::fenced_now(world);
        F("Y, mu, nu") -= compute_contr_CM(C_, M_, forced_norms)("Y, mu, nu");
        F.truncate();

        auto t1_CM = mpqc::fenced_now(world);
        t_CM = mpqc::duration_in_s(t0_CM, t1_CM);
      }
      t1 = mpqc::fenced_now(world);
      t_F = mpqc::duration_in_s(t0, t1);

      mpqc::detail::print_size_info(F, "F(Y,μ,ρ)");

      // permute basis indices
      t0 = mpqc::fenced_now(world);
      Q("Y, nu, rho") = Q("Y, rho, nu");
      F("Y, nu, mu") = F("Y, mu, nu");
      t1 = mpqc::fenced_now(world);
      t_permute = mpqc::duration_in_s(t0, t1);

      // compute K(μ_0, ν_R) = F(Y, ρ_Rj, μ_0) Q(Y, ρ_Rj, ν_R)
      t0 = mpqc::fenced_now(world);
      array_type K_unsymm;
      K_unsymm = compute_contr_FQ(F, Q);
      t1 = mpqc::fenced_now(world);
      t_K = mpqc::duration_in_s(t0, t1);

      t0 = mpqc::fenced_now(world);
      // symmetrize K if required
      if (force_hermiticity_) {
        // force hermiticity of exchange
        K = pbc::detail::symmetrize_matrix(K_unsymm);
      } else {
        // leave exchange as non-hermitian (due to finite lattice range)
        K = K_unsymm;
      }
      t1 = mpqc::fenced_now(world);
      t_symm = mpqc::duration_in_s(t0, t1);
    }

    auto t1_k = mpqc::fenced_now(world);
    auto t_tot = mpqc::duration_in_s(t0_k, t1_k);

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K time decomposition:\n"
                    << "\tQ(Y, ν_R, ρ_Rj) :         " << t_Q << " s\n"
                    << "\tF = 2 * ERI3 - C M:       " << t_F << " s\n"
                    << "\t  Eval E(Y, μ_0, ρ):      " << t_eval_eri3 << " s\n"
                    << "\t  Contract C M:           " << t_CM << " s\n"
                    << "\tPermute F and Q:          " << t_permute << " s\n"
                    << "\tK = F Q:                  " << t_K << " s\n"
                    << "\tSymmetrize K:             " << t_symm << " s\n"
                    << "\nTotal K builder time:     " << t_tot << " s"
                    << std::endl;
    }

    return K;
  }

  /*!
   * \brief This computes forced norms of F(Y_Ry, μ_0, ν_Rν) based on the
   * sparsity of Q(Y_Ry, ρ_Rρ, ν_Rν).
   *
   * Note that K(μ_0, ρ_Rρ) = F(Y_Ry, μ_0, ν_Rν) Q(Y_Ry, ρ_Rρ, ν_Rν). For
   * specific Y_Ry and ν_Rν, no need to compute F(Y_Ry, μ_0, ν_Rν) if all ρ_Rρ
   * in Q(Y_Ry, ρ_Rρ, ν_Rν) give zero. Translational symmetry of Q is used.
   *
   * \param in_norms norms of translational invariant Q,
   * i.e. Q(Y_(Ry-Rρ), ρ_0, ν_(Rν-Rρ))
   * \param out_norms norms of F(Y_Ry, μ_0, ν_Rν). Only correct ranges of F
   * norms are needed.
   * \return forced norms of F
   */
  TA::Tensor<float> force_F_norms(TA::Tensor<float> const &in_norms,
                                  TA::Tensor<float> const &out_norms) {
    const auto ext_F = out_norms.range().extent();
    const auto ntiles_Y = ext_F[0];
    const auto ntiles_nu = ext_F[2];

    using SigPair = std::pair<size_t, size_t>;
    std::unordered_set<SigPair, boost::hash<SigPair>> Y_nu;
    Y_nu.reserve(ntiles_Y * ntiles_nu);

    using ::mpqc::detail::direct_3D_idx;
    using ::mpqc::detail::direct_ord_idx;
    using ::mpqc::detail::is_in_lattice_range;

    // determine Y-ν pairs that are necessary to compute in F(Y_Ry, μ_0, ν_Rν)
    for (auto RY_ord = int64_t(0); RY_ord != RY_size_; ++RY_ord) {
      const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
      for (auto Y = 0ul; Y != ntiles_per_uc_; ++Y) {
        const size_t Y_in_F = Y + RY_ord * ntiles_per_uc_;
        for (auto R1_ord = int64_t(0); R1_ord != R_size_; ++R1_ord) {
          const auto R1_3D = direct_3D_idx(R1_ord, R_max_);
          for (auto nu = 0ul; nu != ntiles_per_uc_; ++nu) {
            const size_t nu_in_F = nu + R1_ord * ntiles_per_uc_;

            SigPair Y_nu_pair(Y_in_F, nu_in_F);
            auto Y_nu_exist = false;
            for (auto R2_ord = int64_t(0); R2_ord != Rrho_size_; ++R2_ord) {
              const auto R2_3D = direct_3D_idx(R2_ord, Rrho_max_);
              const auto RYm2_3D = RY_3D - R2_3D;
              if (!is_in_lattice_range(RYm2_3D, R_max_)) {
                continue;
              }

              const auto R1m2_3D = R1_3D - R2_3D;
              if (!is_in_lattice_range(R1m2_3D, R1m2_max_)) {
                continue;
              }

              const auto R1m2_ord = direct_ord_idx(R1m2_3D, R1m2_max_);
              const auto RYm2_ord = direct_ord_idx(RYm2_3D, R_max_);
              const size_t Y_in_Q = Y + RYm2_ord * ntiles_per_uc_;
              const size_t nu_in_Q = nu + R1m2_ord * ntiles_per_uc_;
              for (auto rho = 0ul; rho != ntiles_per_uc_; ++rho) {
                const auto val = in_norms(Y_in_Q, rho, nu_in_Q);
                if (val > force_shape_threshold_) {
                  Y_nu.insert(Y_nu_pair);
                  Y_nu_exist = true;
                  break;
                }
              }
              if (Y_nu_exist) {
                break;
              }
            }
          }
        }
      }
    }

    const auto &out_range = out_norms.range();
    TA::Tensor<float> out(out_range, 0.0);

    // force F(Y_Ry, μ_0, ν_Rν) norms
    for (auto mu = 0ul; mu != ntiles_per_uc_; ++mu) {
      for (auto const &Y_nu_pair : Y_nu) {
        const auto Y = Y_nu_pair.first;
        const auto nu = Y_nu_pair.second;
        out(Y, mu, nu) = std::numeric_limits<float>::max();
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
   * \brief This computes translationally invariant part of Q(Y_Ry, ρ_Rρ, ν_Rν),
   * i.e. Q(Y_(Ry-Rρ), ρ_0, ν_(Rν-Rρ)).
   *
   * Note Q(Y_Ry, ρ_Rρ, ν_Rν) = C(Y_Ry, ρ_Rρ, σ_(Rρ+Rσ)) D(ν_Rν, σ_(Rρ+Rσ)),
   * thus Q(Y_(Ry-Rρ), ρ_0, ν_(Rν-Rρ)) = C(Y_(Ry-Rρ), ρ_0, σ_Rσ) D(ν_0,
   * σ_(Rρ+Rσ-Rν)). Translational symmetry of C is also used. Y is now on the
   * center of either ρ_0 or σ_Rσ.
   *
   * \param C translationally invariant C, i.e. C(Y_(Ry-Rρ), ρ_0, σ_Rσ)
   * \param D density matrix
   * \return
   */
  array_type compute_Q(const array_type &C, const array_type &D) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();

    array_type C_repl, D_repl;
    C_repl("Y, rho, sig") = C("Y, rho, sig");
    D_repl("nu, sig") = D("nu, sig");
    C_repl.make_replicated();
    D_repl.make_replicated();
    world.gop.fence();  // must wait till all replicating is finished

    using ::mpqc::detail::direct_3D_idx;
    using ::mpqc::detail::direct_ord_idx;
    using ::mpqc::detail::direct_vector;
    using ::mpqc::detail::is_in_lattice_range;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

    const auto &Dnorms = D_repl.shape().data();
    const auto &Cnorms = C_repl.shape().data();
    auto task_id = 0ul;
    for (auto R1m2_ord = int64_t(0); R1m2_ord != R1m2_size_; ++R1m2_ord) {
      const auto R1m2_3D = direct_3D_idx(R1m2_ord, R1m2_max_);
      for (auto R3_ord = int64_t(0); R3_ord != R_size_; ++R3_ord) {
        const auto R3_3D = direct_3D_idx(R3_ord, R_max_);
        const auto RD_3D = R3_3D - R1m2_3D;  // unit cell index for D: Rρ+Rσ-Rν
        if (!is_in_lattice_range(RD_3D, truncated_RD_max_)) {
          continue;
        }

        const auto RD_ord = direct_ord_idx(RD_3D, RD_max_);
        for (auto rho = 0ul; rho != ntiles_per_uc_; ++rho) {
          const auto Y_iij = rho + ref_uc_ord_ * ntiles_per_uc_;
          for (auto sigma = 0ul; sigma != ntiles_per_uc_; ++sigma) {
            const auto sigma_in_C = sigma + R3_ord * ntiles_per_uc_;
            const auto sigma_in_D = sigma + RD_ord * ntiles_per_uc_;
            const auto Y_jij = sigma_in_C;
            for (auto nu = 0ul; nu != ntiles_per_uc_; ++nu) {
              const std::array<size_t, 2> idx_D{
                  {size_t(nu), size_t(sigma_in_D)}};
              if (Dnorms(idx_D) < density_threshold_) {
                continue;
              }

              const auto nu_in_Q = nu + R1m2_ord * ntiles_per_uc_;
              const std::array<size_t, 3> idx_C_jij{
                  {size_t(Y_jij), size_t(rho), size_t(sigma_in_C)}};
              const std::array<size_t, 3> idx_Q_iij{
                  {size_t(Y_iij), size_t(rho), size_t(nu_in_Q)}};
              const std::array<size_t, 3> idx_Q_jij{
                  {size_t(Y_jij), size_t(rho), size_t(nu_in_Q)}};

              if (rho == sigma && R3_ord == ref_uc_ord_) {
                if (C_repl.is_zero(idx_C_jij) ||
                    (Cnorms(idx_C_jij) * Dnorms(idx_D) <
                     cadf_contr_threshold_)) {
                  continue;
                }

                if (task_id % nproc == me) {
                  auto C_jij = C_repl.find(idx_C_jij);
                  auto D_tile = D_repl.find(idx_D);
                  WorldObject_::task(me,
                                     &PeriodicCADFKBuilder_::compute_Q_task_ii,
                                     C_jij, D_tile, idx_Q_jij);
                }
                task_id++;

              } else {
                const std::array<size_t, 3> idx_C_iij{
                    {size_t(Y_iij), size_t(rho), size_t(sigma_in_C)}};
                bool compute_iij = (Cnorms(idx_C_iij) * Dnorms(idx_D) >=
                    cadf_contr_threshold_);
                bool compute_jij = (Cnorms(idx_C_jij) * Dnorms(idx_D) >=
                    cadf_contr_threshold_);
                if ((C_repl.is_zero(idx_C_iij) && C_repl.is_zero(idx_C_jij)) ||
                    (!compute_iij && !compute_jij)) {
                  continue;
                }

                if (task_id % nproc == me) {
                  auto C_iij = C_repl.find(idx_C_iij);
                  auto C_jij = C_repl.find(idx_C_jij);
                  auto D_tile = D_repl.find(idx_D);
                  WorldObject_::task(me,
                                     &PeriodicCADFKBuilder_::compute_Q_task_ij,
                                     C_iij, C_jij, D_tile, idx_Q_iij, idx_Q_jij,
                                     compute_iij, compute_jij);
                }
                task_id++;
              }
            }
          }
        }
      }
    }

    world.gop.fence();

    const auto ntiles_rho = Q_trange_.dim(1).tile_extent();
    const auto ntiles_nu = Q_trange_.dim(2).tile_extent();
    if (world.size() > 1) {
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
        const auto i_stride = ntiles_rho * ntiles_nu;
        const auto j_stride = ntiles_nu;
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
        if (shape[global_tile.first] >= cadf_sparse_threshold_)
          result.set(global_tile.first, global_tile.second);
      }
      result.fill_local(0.0, true);
      global_contr_tiles_.clear();
      world.gop.fence();

      result.truncate();
      world.gop.fence();

      return result;

    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 3>, double>> local_tile_norms;
        const auto i_stride = ntiles_rho * ntiles_nu;
        const auto j_stride = ntiles_nu;
        for (const auto &local_tile : local_contr_tiles_) {
          const auto tile_ord = local_tile.first;
          const auto i = tile_ord / i_stride;
          const auto jk = tile_ord % i_stride;
          const auto j = jk / j_stride;
          const auto k = jk % j_stride;
          const auto norm = local_tile.second.norm();
          local_tile_norms.push_back(
              std::make_pair(std::array<size_t, 3>{{i, j, k}}, norm));
        }
        shape = decltype(shape)(world, local_tile_norms, Q_trange_);
      }

      array_type result(world, Q_trange_, shape, Q_pmap_);
      for (const auto &local_tile : local_contr_tiles_) {
        if (shape[local_tile.first] >= cadf_sparse_threshold_)
          result.set(local_tile.first, local_tile.second);
      }
      result.fill_local(0.0, true);
      local_contr_tiles_.clear();
      world.gop.fence();

      result.truncate();
      world.gop.fence();

      return result;
    }
  }

  /*!
   * \brief This computes contraction between a C tile and a D tile.
   * \param C_tile tile of C
   * \param D_tile tile of D
   * \param tile_idx result index
   */
  void compute_Q_task_ii(Tile C_tile, Tile D_tile,
                         std::array<size_t, 3> tile_idx) {
    const auto ext_C = C_tile.range().extent();
    const auto ext_D = D_tile.range().extent();
    assert(ext_C[2] == ext_D[1]);

    const auto tile_Y = tile_idx[0];
    const auto tile_rho = tile_idx[1];
    const auto tile_nu = tile_idx[2];

    const auto &rng_Y = Q_trange_.dim(0).tile(tile_Y);
    const auto &rng_rho = Q_trange_.dim(1).tile(tile_rho);
    const auto &rng_nu = Q_trange_.dim(2).tile(tile_nu);

    const auto rng_Y_size = rng_Y.second - rng_Y.first;
    const auto rng_rho_size = rng_rho.second - rng_rho.first;
    const auto rng_nu_size = rng_nu.second - rng_nu.first;
    assert(rng_Y_size == ext_C[0]);
    assert(rng_rho_size == ext_C[1]);
    assert(rng_nu_size == ext_D[0]);

    const auto result_rng = TA::Range({rng_Y, rng_rho, rng_nu});
    Tile result_tile(result_rng, 0.0);

    TA::math::GemmHelper gh(madness::cblas::NoTrans, madness::cblas::Trans,
                            result_tile.range().rank(), C_tile.range().rank(),
                            D_tile.range().rank());

    int m, k, n;
    gh.compute_matrix_sizes(m, n, k, C_tile.range(), D_tile.range());
    const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
    const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

    // Notice that we reversed notrans and trans. This is because Lapack
    // expects col major matrices.
    TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, C_tile.data(),
                   lda, D_tile.data(), ldb, 0.0, result_tile.data(), n);

    const auto ntiles_rho = Q_trange_.dim(1).tile_extent();
    const auto ntiles_nu = Q_trange_.dim(2).tile_extent();
    const auto ord =
        tile_Y * ntiles_rho * ntiles_nu + tile_rho * ntiles_nu + tile_nu;

    PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
  }

  /*!
   * \brief This computes contractions between two C tiles (C_iij and C_jij) and
   * a D tile, i.e. Q_iij = C_iij*D and Q_jij = C_jij*D.
   * \param C_iij tile of C_iij
   * \param C_jij tile of C_jij
   * \param D_tile tile of D
   * \param idx_Q_iij result index of Q_iij
   * \param idx_Q_jij result index of Q_jij
   */
  void compute_Q_task_ij(Tile C_iij, Tile C_jij, Tile D_tile,
                         std::array<size_t, 3> idx_Q_iij,
                         std::array<size_t, 3> idx_Q_jij, bool compute_iij,
                         bool compute_jij) {
    const auto ext_C_iij = C_iij.range().extent();
    const auto ext_C_jij = C_jij.range().extent();
    const auto ext_D = D_tile.range().extent();
    assert(idx_Q_iij[0] != idx_Q_jij[0] && idx_Q_iij[1] == idx_Q_jij[1] &&
           idx_Q_iij[2] == idx_Q_jij[2]);
    assert(ext_C_iij[2] == ext_D[1] && ext_C_jij[2] == ext_D[1]);

    auto create_Q_tile = [&](Tile &C, Tile &D, std::array<size_t, 3> &idx_Q) {
      const auto tile_Y = idx_Q[0];
      const auto tile_rho = idx_Q[1];
      const auto tile_nu = idx_Q[2];

      const auto &rng_Y = Q_trange_.dim(0).tile(tile_Y);
      const auto &rng_rho = Q_trange_.dim(1).tile(tile_rho);
      const auto &rng_nu = Q_trange_.dim(2).tile(tile_nu);

      const auto rng_Y_size = rng_Y.second - rng_Y.first;
      const auto rng_rho_size = rng_rho.second - rng_rho.first;
      const auto rng_nu_size = rng_nu.second - rng_nu.first;
      const auto ext_C = C.range().extent();
      const auto ext_D = D.range().extent();
      assert(rng_Y_size == ext_C[0]);
      assert(rng_rho_size == ext_C[1]);
      assert(rng_nu_size == ext_D[0]);

      const auto result_rng = TA::Range({rng_Y, rng_rho, rng_nu});
      Tile result_tile(result_rng, 0.0);

      TA::math::GemmHelper gh(madness::cblas::NoTrans, madness::cblas::Trans,
                              result_tile.range().rank(), C.range().rank(),
                              D.range().rank());
      int m, k, n;
      gh.compute_matrix_sizes(m, n, k, C.range(), D.range());
      const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
      const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

      // Notice that we reversed notrans and trans. This is because Lapack
      // expects col major matrices.
      TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, C.data(), lda,
                     D.data(), ldb, 0.0, result_tile.data(), n);

      const auto ntiles_rho = Q_trange_.dim(1).tile_extent();
      const auto ntiles_nu = Q_trange_.dim(2).tile_extent();
      const auto ord =
          tile_Y * ntiles_rho * ntiles_nu + tile_rho * ntiles_nu + tile_nu;

      PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
    };

    if (compute_iij) {
      create_Q_tile(C_iij, D_tile, idx_Q_iij);
    }

    if (compute_jij) {
      create_Q_tile(C_jij, D_tile, idx_Q_jij);
    }
  }

  /*!
   * \brief This computes K(μ_0, ρ_Rρ) = F(Y_Ry, ν_Rν, μ_0) Q(Y_Ry, ν_Rν, ρ_Rρ).
   * Translational symmetry of Q is used.
   *
   * \param F array F
   * \param Q translationally invariant Q, i.e. Q(Y(Ry-Rρ), ν_(Rν-Rρ), ρ_0)
   * \return
   */
  array_type compute_contr_FQ(const array_type &F, const array_type &Q) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();

    using ::mpqc::detail::direct_3D_idx;
    using ::mpqc::detail::direct_ord_idx;
    using ::mpqc::detail::is_in_lattice_range;

    const auto &Fnorms = F.shape().data();
    const auto &Qnorms = Q.shape().data();
    auto task_id = 0ul;
    for (auto RY_ord = int64_t(0); RY_ord != RY_size_; ++RY_ord) {
      const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
      for (auto R2_ord = int64_t(0); R2_ord != Rrho_size_; ++R2_ord) {
        const auto R2_3D = direct_3D_idx(R2_ord, Rrho_max_);
        const auto RYm2_3D = RY_3D - R2_3D;
        if (!is_in_lattice_range(RYm2_3D, R_max_)) {
          continue;
        }

        const auto RYm2_ord = direct_ord_idx(RYm2_3D, R_max_);
        for (auto R1_ord = int64_t(0); R1_ord != R_size_; ++R1_ord) {
          const auto R1_3D = direct_3D_idx(R1_ord, R_max_);
          const auto R1m2_3D = R1_3D - R2_3D;
          if (!is_in_lattice_range(R1m2_3D, R1m2_max_)) {
            continue;
          }

          const auto R1m2_ord = direct_ord_idx(R1m2_3D, R1m2_max_);
          for (auto Y = 0ul; Y != ntiles_per_uc_; ++Y) {
            const auto Y_in_Q = Y + RYm2_ord * ntiles_per_uc_;
            const auto Y_in_F = Y + RY_ord * ntiles_per_uc_;
            for (auto nu = 0ul; nu != ntiles_per_uc_; ++nu) {
              const auto nu_in_Q = nu + R1m2_ord * ntiles_per_uc_;
              const auto nu_in_F = nu + R1_ord * ntiles_per_uc_;
              for (auto rho = 0ul; rho != ntiles_per_uc_; ++rho) {
                const std::array<size_t, 3> idx_Q{
                    {size_t(Y_in_Q), size_t(nu_in_Q), size_t(rho)}};
                if (Q.is_zero(idx_Q)) {
                  continue;
                }
                const auto rho_in_K = rho + R2_ord * ntiles_per_uc_;
                for (auto mu = 0ul; mu != ntiles_per_uc_; ++mu) {
                  const std::array<size_t, 3> idx_F{
                      {size_t(Y_in_F), size_t(nu_in_F), size_t(mu)}};
                  if (F.is_zero(idx_F) || (Fnorms(idx_F) * Qnorms(idx_Q) <
                      cadf_contr_threshold_)) {
                    continue;
                  }
                  if (task_id % nproc == me) {
                    const std::array<size_t, 2> idx_K{
                        {size_t(mu), size_t(rho_in_K)}};
                    auto F_tile = F.find(idx_F);
                    auto Q_tile = Q.find(idx_Q);
                    WorldObject_::task(
                        me, &PeriodicCADFKBuilder_::compute_contr_FQ_task,
                        F_tile, Q_tile, idx_K);
                  }
                  task_id++;
                }
              }
            }
          }
        }
      }
    }

    world.gop.fence();

    const auto ntiles_rho = result_trange_.dim(1).tile_extent();
    if (world.size() > 1) {
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
          const auto i = tile_ord / ntiles_rho;
          const auto j = tile_ord % ntiles_rho;
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

    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 2>, double>> local_tile_norms;
        for (const auto &local_tile : local_contr_tiles_) {
          const auto tile_ord = local_tile.first;
          const auto i = tile_ord / ntiles_rho;
          const auto j = tile_ord % ntiles_rho;
          const auto norm = local_tile.second.norm();
          local_tile_norms.push_back(
              std::make_pair(std::array<size_t, 2>{{i, j}}, norm));
        }
        shape = decltype(shape)(world, local_tile_norms, result_trange_);
      }

      array_type result(world, result_trange_, shape, result_pmap_);
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

    const auto tile_mu = tile_idx[0];
    const auto tile_rho = tile_idx[1];

    const auto &rng_mu = result_trange_.dim(0).tile(tile_mu);
    const auto &rng_rho = result_trange_.dim(1).tile(tile_rho);

    const auto rng_mu_rng = rng_mu.second - rng_mu.first;
    const auto rng_rho_rng = rng_rho.second - rng_rho.first;

    assert(rng_mu_rng == ext_F[2]);
    assert(rng_rho_rng == ext_Q[2]);

    const auto result_rng = TA::Range({rng_mu, rng_rho});
    Tile result_tile(result_rng, 0.0);

    TA::math::GemmHelper gh(madness::cblas::Trans, madness::cblas::NoTrans,
                            result_tile.range().rank(), F_tile.range().rank(),
                            Q_tile.range().rank());

    int m, k, n;
    gh.compute_matrix_sizes(m, n, k, F_tile.range(), Q_tile.range());
    const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
    const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

    // Notice that we reversed notrans and trans. This is because Lapack
    // expects col major matrices.
    TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, F_tile.data(),
                   lda, Q_tile.data(), ldb, 0.0, result_tile.data(), n);

    const auto ntiles_rho = result_trange_.dim(1).tile_extent();
    const auto ord = tile_mu * ntiles_rho + tile_rho;

    PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
  }

  /*!
   * \brief This computes contraction C(X_Rx, μ_0, ν_Rν) M(X_Rx, Y_Ry).
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
    C_repl("X, mu, nu") = C("X, mu, nu");
    M_repl("X, Y") = M("X, Y");
    C_repl.make_replicated();
    world.gop.fence();
    M_repl.make_replicated();
    world.gop.fence();  // must wait till all replicating is finished

    using ::mpqc::detail::direct_3D_idx;
    using ::mpqc::detail::direct_ord_idx;
    using ::mpqc::detail::is_in_lattice_range;

    const auto &Cnorms = C_repl.shape().data();
    const auto &Mnorms = M_repl.shape().data();
    auto task_id = 0ul;
    for (auto R1_ord = int64_t(0); R1_ord != R_size_; ++R1_ord) {
      const auto R1_3D = direct_3D_idx(R1_ord, R_max_);
      for (auto mu = 0ul; mu != ntiles_per_uc_; ++mu) {
        const auto X_in_C_iij = mu + ref_uc_ord_ * ntiles_per_uc_;
        for (auto nu = 0ul; nu != ntiles_per_uc_; ++nu) {
          const auto nu_R1 = nu + R1_ord * ntiles_per_uc_;
          const auto X_in_C_jij = nu_R1;
          const std::array<size_t, 3> idx_C_iij{
              {size_t(X_in_C_iij), size_t(mu), size_t(nu_R1)}};
          const std::array<size_t, 3> idx_C_jij{
              {size_t(X_in_C_jij), size_t(mu), size_t(nu_R1)}};

          for (auto RY_ord = int64_t(0); RY_ord != RY_size_; ++RY_ord) {
            const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
            const auto RY_ord_in_M = direct_ord_idx(RY_3D, RYmRX_max_);
            const auto RYm1_3D = RY_3D - R1_3D;
            if (!is_in_lattice_range(RYm1_3D, RYmRX_max_)) {
              continue;
            }

            const auto RYm1_ord = direct_ord_idx(RYm1_3D, RYmRX_max_);
            for (auto Y = 0ul; Y != ntiles_per_uc_; ++Y) {
              const auto Y_in_F = Y + RY_ord * ntiles_per_uc_;
              const std::array<size_t, 3> idx_F{
                  {size_t(Y_in_F), size_t(mu), size_t(nu_R1)}};
              if (force_norms(idx_F) < force_shape_threshold_) {
                continue;
              }

              const auto Y_in_M_iij = Y + RY_ord_in_M * ntiles_per_uc_;
              const auto Y_in_M_jij = Y + RYm1_ord * ntiles_per_uc_;
              const std::array<size_t, 2> idx_M_iij{
                  {size_t(mu), size_t(Y_in_M_iij)}};
              const std::array<size_t, 2> idx_M_jij{
                  {size_t(nu), size_t(Y_in_M_jij)}};

              if (mu == nu && R1_ord == ref_uc_ord_) {
                if (C_repl.is_zero(idx_C_iij) || M_repl.is_zero(idx_M_iij) ||
                    (Cnorms(idx_C_iij) * Mnorms(idx_M_iij) <
                        cadf_contr_threshold_)) {
                  continue;
                }

                if (task_id % nproc == me) {
                  auto C_iij = C_repl.find(idx_C_iij);
                  auto M_iij = M_repl.find(idx_M_iij);
                  WorldObject_::task(
                      me, &PeriodicCADFKBuilder_::compute_contr_CM_task_ii,
                      C_iij, M_iij, idx_F);
                }
                task_id++;
              } else {
                bool compute_iij = (Cnorms(idx_C_iij) * Mnorms(idx_M_iij) >=
                    cadf_contr_threshold_);
                bool compute_jij = (Cnorms(idx_C_jij) * Mnorms(idx_M_jij) >=
                    cadf_contr_threshold_);
                if (((C_repl.is_zero(idx_C_iij) || M_repl.is_zero(idx_M_iij)) &&
                     (C_repl.is_zero(idx_C_jij) ||
                      M_repl.is_zero(idx_M_jij))) ||
                    (!compute_iij && !compute_jij)) {
                  continue;
                }

                if (task_id % nproc == me) {
                  auto C_iij = C_repl.find(idx_C_iij);
                  auto C_jij = C_repl.find(idx_C_jij);
                  auto M_iij = M_repl.find(idx_M_iij);
                  auto M_jij = M_repl.find(idx_M_jij);
                  WorldObject_::task(
                      me, &PeriodicCADFKBuilder_::compute_contr_CM_task_ij,
                      C_iij, C_jij, M_iij, M_jij, idx_F, compute_iij,
                      compute_jij);
                }
                task_id++;
              }
            }
          }
        }
      }
    }

    world.gop.fence();

    const auto ntiles_mu = F_trange_.dim(1).tile_extent();
    const auto ntiles_nu = F_trange_.dim(2).tile_extent();
    if (world.size() > 1) {
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
        const auto i_stride = ntiles_mu * ntiles_nu;
        const auto j_stride = ntiles_nu;
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
        if (shape[global_tile.first] >= cadf_sparse_threshold_)
          result.set(global_tile.first, global_tile.second);
      }
      result.fill_local(0.0, true);
      global_contr_tiles_.clear();
      world.gop.fence();

      result.truncate();
      world.gop.fence();

      return result;

    } else {
      typename Policy::shape_type shape;
      // compute the shape, if sparse
      if (!decltype(shape)::is_dense()) {
        // extract local contribution to the shape of G, construct global shape
        std::vector<std::pair<std::array<size_t, 3>, double>> local_tile_norms;
        const auto i_stride = ntiles_mu * ntiles_nu;
        const auto j_stride = ntiles_nu;
        for (const auto &local_tile : local_contr_tiles_) {
          const auto tile_ord = local_tile.first;
          const auto i = tile_ord / i_stride;
          const auto jk = tile_ord % i_stride;
          const auto j = jk / j_stride;
          const auto k = jk % j_stride;
          const auto norm = local_tile.second.norm();
          local_tile_norms.push_back(
              std::make_pair(std::array<size_t, 3>{{i, j, k}}, norm));
        }
        shape = decltype(shape)(world, local_tile_norms, F_trange_);
      }

      array_type result(world, F_trange_, shape, F_pmap_);
      for (const auto &local_tile : local_contr_tiles_) {
        if (shape[local_tile.first] >= cadf_sparse_threshold_)
          result.set(local_tile.first, local_tile.second);
      }
      result.fill_local(0.0, true);
      local_contr_tiles_.clear();
      world.gop.fence();

      result.truncate();
      world.gop.fence();

      return result;
    }
  }

  /*!
   * \brief This computes contraction between a C tile and a M tile.
   * \param C_tile tile of C
   * \param M_tile tile of M
   * \param tile_idx result index
   */
  void compute_contr_CM_task_ii(Tile C_tile, Tile M_tile,
                                std::array<size_t, 3> tile_idx) {
    const auto ext_C = C_tile.range().extent();
    const auto ext_M = M_tile.range().extent();
    assert(ext_M[0] == ext_C[0]);

    const auto tile_Y = tile_idx[0];
    const auto tile_mu = tile_idx[1];
    const auto tile_nu = tile_idx[2];

    const auto &rng_Y = F_trange_.dim(0).tile(tile_Y);
    const auto &rng_mu = F_trange_.dim(1).tile(tile_mu);
    const auto &rng_nu = F_trange_.dim(2).tile(tile_nu);

    const auto rng_Y_size = rng_Y.second - rng_Y.first;
    const auto rng_mu_size = rng_mu.second - rng_mu.first;
    const auto rng_nu_size = rng_nu.second - rng_nu.first;

    assert(rng_Y_size == ext_M[1]);
    assert(rng_mu_size == ext_C[1]);
    assert(rng_nu_size == ext_C[2]);

    const auto result_rng = TA::Range({rng_Y, rng_mu, rng_nu});
    Tile result_tile(result_rng, 0.0);

    TA::math::GemmHelper gh(madness::cblas::Trans, madness::cblas::NoTrans,
                            result_tile.range().rank(), M_tile.range().rank(),
                            C_tile.range().rank());

    int m, k, n;
    gh.compute_matrix_sizes(m, n, k, M_tile.range(), C_tile.range());
    const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
    const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

    // Notice that we reversed notrans and trans. This is because Lapack
    // expects col major matrices.
    TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, M_tile.data(),
                   lda, C_tile.data(), ldb, 0.0, result_tile.data(), n);

    const auto ntiles_mu = F_trange_.dim(1).tile_extent();
    const auto ntiles_nu = F_trange_.dim(2).tile_extent();
    const auto ord =
        tile_Y * ntiles_mu * ntiles_nu + tile_mu * ntiles_nu + tile_nu;

    PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
  }

  /*!
   * \brief This computes contractions between two C tiles (C_iij and C_jij) and
   * two M tiles (M_iij and M_jij), i.e. F = C_iij*M_iij + C_jij*M_jij.
   * \param C_iij tile of C_iij
   * \param C_jij tile of C_jij
   * \param M_iij tile of M_iij
   * \param M_jij tile of M_jij
   * \param idx_F result index of F
   */
  void compute_contr_CM_task_ij(Tile C_iij, Tile C_jij, Tile M_iij, Tile M_jij,
                                std::array<size_t, 3> idx_F, bool compute_iij,
                                bool compute_jij) {
    const auto ext_C_iij = C_iij.range().extent();
    const auto ext_C_jij = C_jij.range().extent();
    const auto ext_M_iij = M_iij.range().extent();
    const auto ext_M_jij = M_jij.range().extent();
    assert(ext_C_iij[1] == ext_C_jij[1] && ext_C_iij[2] == ext_C_jij[2]);
    assert(ext_M_iij[1] == ext_M_jij[1]);
    assert(ext_M_iij[0] == ext_C_iij[0] && ext_M_jij[0] == ext_C_jij[0]);

    const auto tile_Y = idx_F[0];
    const auto tile_mu = idx_F[1];
    const auto tile_nu = idx_F[2];

    const auto &rng_Y = F_trange_.dim(0).tile(tile_Y);
    const auto &rng_mu = F_trange_.dim(1).tile(tile_mu);
    const auto &rng_nu = F_trange_.dim(2).tile(tile_nu);

    const auto rng_Y_size = rng_Y.second - rng_Y.first;
    const auto rng_mu_size = rng_mu.second - rng_mu.first;
    const auto rng_nu_size = rng_nu.second - rng_nu.first;

    assert(rng_Y_size == ext_M_iij[1]);
    assert(rng_mu_size == ext_C_iij[1]);
    assert(rng_nu_size == ext_C_iij[2]);

    const auto result_rng = TA::Range({rng_Y, rng_mu, rng_nu});
    Tile result_tile(result_rng, 0.0);

    auto add_to_result_tile = [&](Tile &C, Tile &M) {
      // set up gemm helper
      TA::math::GemmHelper gh(madness::cblas::Trans, madness::cblas::NoTrans,
                              result_tile.range().rank(), M.range().rank(),
                              C.range().rank());

      int m, k, n;
      gh.compute_matrix_sizes(m, n, k, M.range(), C.range());
      const auto lda = (gh.left_op() == madness::cblas::NoTrans ? k : m);
      const auto ldb = (gh.right_op() == madness::cblas::NoTrans ? n : k);

      // Notice that we reversed notrans and trans. This is because Lapack
      // expects col major matrices.
      TA::math::gemm(gh.left_op(), gh.right_op(), m, n, k, 1.0, M.data(), lda,
                     C.data(), ldb, 1.0, result_tile.data(), n);
    };

    if (compute_iij) {
      add_to_result_tile(C_iij, M_iij);
    }

    if (compute_jij) {
      add_to_result_tile(C_jij, M_jij);
    }

    const auto ntiles_mu = F_trange_.dim(1).tile_extent();
    const auto ntiles_nu = F_trange_.dim(2).tile_extent();
    const auto ord =
        tile_Y * ntiles_mu * ntiles_nu + tile_mu * ntiles_nu + tile_nu;

    PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
  }

  void accumulate_global_task(Tile arg_tile, size_t tile_ord) {
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

  void accumulate_local_task(Tile arg_tile, size_t tile_ord) {
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
        for (auto nu = 0; nu < ext[2]; ++nu) {
          M(X * ext[1] + mu, nu) = norms(X, mu, nu);
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
        for (auto nu = 0; nu < ext[2]; ++nu) {
          M(X, mu * ext[2] + nu) = norms(X, mu, nu);
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

    RowMatrixXd M_eig = math::array_to_eigen(M);

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
   * \brief This updates RD-dependent variables
   * \param RD_max
   */
  void update_RD_dependent_variables(const Vector3i &RD_max) {
    if (print_detail_) {
      ExEnv::out0() << "\nUpdating RD-dependent variables:" << std::endl;
    }
    using ::mpqc::detail::direct_ord_idx;
    using ::mpqc::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::make_engine_pool;
    auto &world = this->get_world();
    time_point t0, t1;

    // compute M(X_Rx, Y_Ry) = M(X_0, Y_(Ry-Rx))
    t0 = mpqc::fenced_now(world);
    {
      // compute lattice range for Rρ in |ρ_Rρ σ_(Rρ+Rσ))
      Rrho_max_ = RD_max + 2 * R_max_;
      Rrho_size_ = 1 + direct_ord_idx(Rrho_max_, Rrho_max_);
      // Y is on the center of either ρ_Rρ or σ_(Rρ+Rσ). Thus Y_Ry should have
      // the same lattice range as σ_(Rρ+Rσ).
      RY_max_ = Rrho_max_;
      RY_size_ = 1 + direct_ord_idx(RY_max_, RY_max_);
      RYmRX_max_ = Rrho_max_;
      auto shifted_Y_dfbs =
          shift_basis_origin(*dfbs_, Vector3d::Zero(), RYmRX_max_, dcell_);
      M_ = compute_eri2(world, *dfbs_, *shifted_Y_dfbs);
    }
    t1 = mpqc::fenced_now(world);
    auto t_M = mpqc::duration_in_s(t0, t1);

    // make direct integral eri3_ = (μ_0 ν_Rν | Y_Ry)
    t0 = mpqc::fenced_now(world);
    {
      auto Y_dfbs =
          shift_basis_origin(*dfbs_, Vector3d::Zero(), RY_max_, dcell_);
      auto bs_array = utility::make_array_of_refs(*Y_dfbs, *obs_, *basisR_);
      auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs, *obs_, *basisR_}};

      auto oper_type = libint2::Operator::coulomb;
      auto screen_engine =
          make_engine_pool(oper_type, bs_array, libint2::BraKet::xx_xx);
      auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
          lcao::gaussian::create_schwarz_screener(
              world, screen_engine, bs_vector, screen_threshold_));
      auto engine =
          make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);
      eri3_ = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                      std::move(screener));

      // make TiledRange and Pamp of F(Y_Ry, μ_0, ν_Rν) (same as eri3)
      F_trange_ = ::mpqc::lcao::gaussian::detail::create_trange(bs_vector);
      auto F_tvolume = F_trange_.tiles_range().volume();
      F_pmap_ = Policy::default_pmap(world, F_tvolume);
    }
    t1 = mpqc::fenced_now(world);
    auto t_direct_eri3 = mpqc::duration_in_s(t0, t1);

    // misc:
    // 1. determine tiled ranges and pmap of the exchange term
    // 2. determine translationally invariant basis, tiled ranges and pmap for
    // Q(Y, rho, nu)
    t0 = mpqc::fenced_now(world);
    {
      // make TiledRange and Pmap of exchange K(μ_0, ρ_Rρ)
      auto basisRrho =
          shift_basis_origin(*obs_, Vector3d::Zero(), Rrho_max_, dcell_);
      result_trange_ = ::mpqc::lcao::gaussian::detail::create_trange(
          lcao::gaussian::BasisVector{{*obs_, *basisRrho}});
      auto tvolume = result_trange_.tiles_range().volume();
      result_pmap_ = Policy::default_pmap(world, tvolume);

      // make basis of Q(Y_(Ry-Rρ), ρ_0, ν(Rν-Rρ))
      R1m2_max_ = RD_max + R_max_;
      R1m2_size_ = 1 + direct_ord_idx(R1m2_max_, R1m2_max_);
      auto Q_bs_Y =
          shift_basis_origin(*dfbs_, Vector3d::Zero(), R_max_, dcell_);
      auto Q_bs_nu =
          shift_basis_origin(*obs_, Vector3d::Zero(), R1m2_max_, dcell_);

      // make TiledRange and Pmap of Q(Y_(Ry-Rρ), ρ_0, ν(Rν-Rρ))
      Q_trange_ = ::mpqc::lcao::gaussian::detail::create_trange(
          lcao::gaussian::BasisVector{{*Q_bs_Y, *obs_, *Q_bs_nu}});
      auto Q_tvolume = Q_trange_.tiles_range().volume();
      Q_pmap_ = Policy::default_pmap(world, Q_tvolume);
    }
    t1 = mpqc::fenced_now(world);
    auto t_misc = mpqc::duration_in_s(t0, t1);

    if (print_detail_) {
      ExEnv::out0() << "\tM(X_0, Y_(Ry-Rx)):        " << t_M << " s\n"
                    << "\tdirect ERI3:              " << t_direct_eri3 << " s\n"
                    << "\tmisc:                     " << t_misc << " s"
                    << std::endl;
    }
  }
};

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
