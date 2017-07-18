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
  using const_data_ptr = typename Tile::allocator_type::const_pointer;

  using DirectTArray = typename Factory::DirectTArray;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;
  using Qmatrix = ::mpqc::lcao::gaussian::Qmatrix;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
  using shellpair_val_list_t =
      std::unordered_map<size_t, std::vector<std::pair<size_t, double>>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;

  using WorldObject_ =
      madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>>;
  using PeriodicCADFKBuilder_ = PeriodicCADFKBuilder<Tile, Policy, Factory>;
  using FutureVectorTile = TA::Future<std::vector<TA::Future<Tile>>>;

  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;

  template <int rank>
  using norm_type = std::vector<std::pair<std::array<int, rank>, float>>;

  PeriodicCADFKBuilder(madness::World &world, Factory &ao_factory,
                       const double force_shape_threshold = 0.0)
      : WorldObject_(world), ao_factory_(ao_factory) {
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    print_detail_ = ao_factory_.print_detail();
    shell_pair_threshold_ = ao_factory_.shell_pair_threshold();
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
    RJ_max_ = R_max_;  // replace RJ range with R range
    RD_max_ = ao_factory_.RD_max();
    R_size_ = ao_factory_.R_size();
    RJ_size_ = R_size_;  // replace RJ range with R range
    RD_size_ = ao_factory_.RD_size();

    const Vector3i ref_latt_range = {0, 0, 0};
    const auto natoms_per_uc = ao_factory_.unitcell().natoms();
    Vector3d zero_shift_base(0.0, 0.0, 0.0);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    X_dfbs_ = shift_basis_origin(*dfbs_, zero_shift_base, RJ_max_, dcell_);

    // compute C(X, μ_0, ρ_Rj)
    t0 = mpqc::fenced_now(world);
    {
      basisRJ_ = shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);

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
    M_ = compute_eri2(world, *X_dfbs_, *Y_dfbs_);
    t1 = mpqc::fenced_now(world);
    double t_M = mpqc::duration_in_s(t0, t1);

    auto oper_type = libint2::Operator::coulomb;
    const auto screen = ao_factory_.screen();
    const auto screen_thresh = ao_factory_.screen_threshold();
    auto screen_engine = make_engine_pool(
        oper_type, utility::make_array_of_refs(*dfbs_, *obs_, *obs_),
        libint2::BraKet::xx_xx);

    // make screener and engines for eri3
    t0 = mpqc::fenced_now(world);
    {
      engines_ = make_engine_pool(
          oper_type, utility::make_array_of_refs(*dfbs_, *obs_, *obs_),
          libint2::BraKet::xs_xx);

      eri3_X_dfbs_ =
          shift_basis_origin(*dfbs_, zero_shift_base, RJ_max_ + R_max_, dcell_);
      eri3_bs0_ = obs_;
      eri3_bs1_ = shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);
      auto basis_vector =
          lcao::gaussian::BasisVector{{*eri3_X_dfbs_, *eri3_bs0_, *eri3_bs1_}};
      p_screener_ = lcao::gaussian::detail::make_screener(
          world, screen_engine, basis_vector, screen, screen_thresh);

      eri3_bs0_shell_offset_map_ = compute_shell_offset(*eri3_bs0_);
      eri3_bs1_shell_offset_map_ = compute_shell_offset(*eri3_bs1_);
      Q_X_shell_offset_map_ = compute_shell_offset(*X_dfbs_);

      basisR_ = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);

      eri3_X_trange1_ = eri3_X_dfbs_->create_trange1();
      eri3_bs0_trange1_ = eri3_bs0_->create_trange1();
      eri3_bs1_trange1_ = eri3_bs1_->create_trange1();
      X_dfbs_trange1_ = X_dfbs_->create_trange1();

      // compute significant shell pair list
      sig_shellpair_list_ = parallel_compute_shellpair_list(*obs_, *basisRJ_);
      nu_shellpair_list_ = parallel_compute_shellpair_list(*basisRJ_, *obs_);
      // make a list of significant Rj's as in overlap between μ and ρ_Rj
      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        const auto nshells = obs_->flattened_shells().size();
        const auto shell1_min = nshells * RJ;
        const auto shell1_max = shell1_min + nshells;

        auto is_significant = false;
        for (auto shell0 = 0; shell0 != nshells; ++shell0) {
          for (const auto &shell1 : sig_shellpair_list_[shell0]) {
            if (shell1 >= shell1_min && shell1 < shell1_max) {
              is_significant = true;
              RJ_list_.emplace_back(RJ);
              break;
            }
          }
          if (is_significant) break;
        }
      }

      // make TiledRange of Exchange
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
      Q_trange_ = ::mpqc::lcao::gaussian::detail::create_trange(lcao::gaussian::BasisVector{{*Q_bs_Y_, *Q_bs_nu_, *Q_bs_rho_}});
      auto Q_tvolume = Q_trange_.tiles_range().volume();
      Q_pmap_ = Policy::default_pmap(world, Q_tvolume);

      // make direct integral E_ket_
      {
        auto bs_array =
            utility::make_array_of_refs(*Y_dfbs_, *obs_, *basisRJ_);
        auto bs_vector =
            lcao::gaussian::BasisVector{{*Y_dfbs_, *obs_, *basisRJ_}};

        auto screen_engine = make_engine_pool(
            oper_type, bs_array, libint2::BraKet::xx_xx);
        auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(lcao::gaussian::create_schwarz_screener(world, screen_engine, bs_vector, screen_thresh));

        auto engine =
            make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);

        E_ket_ = lcao::gaussian::direct_sparse_integrals(
            world, engine, bs_vector, std::move(screener));
      }

    }
    t1 = mpqc::fenced_now(world);
    auto t_misc = mpqc::duration_in_s(t0, t1);

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K init time decomposition:\n"
                    << "\tC(X, μ_0, ρ_Rj):     " << t_C_bra << " s\n"
                    << "\tM(X, Y):             " << t_M << " s\n"
                    << "\tengine, screen ... : " << t_misc << " s" << std::endl;
    }
  }

  ~PeriodicCADFKBuilder() { engines_.reset(); }

  array_type operator()(array_type const &D, double target_precision) {
    return compute_K(D, target_precision);
  }

 private:
  Factory &ao_factory_;
  bool print_detail_;
  double force_shape_threshold_;
  double target_precision_ = std::numeric_limits<double>::epsilon();
  double shell_pair_threshold_;
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

  std::unique_ptr<PTC_Builder> three_center_builder_;
  shellpair_list_t sig_shellpair_list_;
  shellpair_list_t nu_shellpair_list_;
  shellpair_val_list_t sig_X_shellpair_val_list_;
  std::vector<int64_t> RJ_list_;
  std::shared_ptr<lcao::Screener> p_screener_;
  Engine engines_;

  std::shared_ptr<Basis> obs_;
  std::shared_ptr<Basis> dfbs_;
  std::shared_ptr<Basis> X_dfbs_;
  std::shared_ptr<Basis> Y_dfbs_;
  std::shared_ptr<Basis> basisRJ_;
  std::shared_ptr<Basis> basisR_;
  std::shared_ptr<Basis> basisRD_;
  std::shared_ptr<Basis> eri3_X_dfbs_;
  std::shared_ptr<Basis> eri3_bs0_;
  std::shared_ptr<Basis> eri3_bs1_;

  std::shared_ptr<Basis> Q_bs_Y_;
  std::shared_ptr<Basis> Q_bs_nu_;
  std::shared_ptr<Basis> Q_bs_rho_;
  TA::TiledRange Q_trange_;
  std::shared_ptr<TA::Pmap> Q_pmap_;
  Vector3i RYmR_max_;
  Vector3i RJmR_max_;

  TA::TiledRange1 basisRD_trange1_;
  TA::TiledRange1 X_dfbs_trange1_;
  TA::TiledRange1 eri3_X_trange1_;
  TA::TiledRange1 eri3_bs0_trange1_;
  TA::TiledRange1 eri3_bs1_trange1_;
  TA::TiledRange result_trange_;
  std::shared_ptr<TA::Pmap> result_pmap_;
  std::unordered_map<size_t, size_t> eri3_bs0_shell_offset_map_;
  std::unordered_map<size_t, size_t> eri3_bs1_shell_offset_map_;
  std::unordered_map<size_t, size_t> Q_sig_shell_offset_map_;
  std::unordered_map<size_t, size_t> Q_X_shell_offset_map_;

  madness::ConcurrentHashMap<std::size_t, Tile> local_contr_tiles_;
  madness::ConcurrentHashMap<std::size_t, Tile> global_contr_tiles_;
  std::atomic<size_t> num_ints_computed_{0};

 private:
  array_type compute_K(const array_type &D, double target_precision) {
    auto &world = ao_factory_.world();
    const auto me = world.rank();

    auto t0_k = mpqc::fenced_now(world);
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    array_type K;

    time_point t0, t1;
    double t_Qket = 0.0;
    double t_eval_Eket = 0.0;
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

      // compute F(Y, μ_0, ρ_Rj) = E(Y, μ_0, ρ_Rj) - C(X, μ_0, ρ_Rj) M(X, Y)
      t0 = mpqc::fenced_now(world);
      array_type F;
      {
        auto t0 = mpqc::fenced_now(world);
        auto F_norms = force_F_norms_new(Q_ket.shape().data(),
                                         E_ket_.array().shape().data());
        auto trange = E_ket_.array().trange();
        TA::SparseShape<float> forced_shape(world, F_norms, trange);

        F("Y, mu, rho") = (E_ket_("Y, mu, rho")).set_shape(forced_shape);
        F.truncate();

        F("Y, mu, rho") = 2.0 * F("Y, mu, rho");
        auto t1 = mpqc::fenced_now(world);
        t_eval_Eket = mpqc::duration_in_s(t0, t1);


        t0 = mpqc::fenced_now(world);
        F("Y, mu, rho") -=
            (M_("X, Y") * C_bra_("X, mu, rho")).set_shape(forced_shape);
        F.truncate();
      }
      t1 = mpqc::fenced_now(world);
      t_F = mpqc::duration_in_s(t0, t1);


      // permute basis indices
      t0 = mpqc::fenced_now(world);
      Q_ket("Y, rho, nu") = Q_ket("Y, nu, rho");
      F("Y, rho, mu") = F("Y, mu, rho");
      t1 = mpqc::fenced_now(world);
      t_permute = mpqc::duration_in_s(t0, t1);

      t0 = mpqc::fenced_now(world);
      K = compute_contr_FQ_v3(F, Q_ket);
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
                    << "\tPermute F and Q_ket:      " << t_permute << " s\n"
                    << "\tK = F Q_ket:              " << t_K << " s\n"
                    << "\nTotal K builder time:     " << t_tot << " s"
                    << std::endl;
    }

    return K;
  }

  TA::Tensor<float> force_F_norms(TA::Tensor<float> const &in,
                                  TA::Tensor<float> const &out_norms,
                                  bool translate) {
    const auto ntiles_Y = Y_dfbs_->nclusters();
    const auto ntiles_rho = ntiles_per_uc_;
    const auto ntiles_nu = basisR_->nclusters();

    using SigPair = std::pair<int64_t, int64_t>;
    std::unordered_set<SigPair, boost::hash<SigPair>> Y_rho;
    Y_rho.reserve(ntiles_Y * ntiles_rho);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    if (translate) {
      for (auto tile_Y = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
        const auto RY_ord = tile_Y / ntiles_per_uc_;
        const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
        const auto tile_Y_in_uc = tile_Y % ntiles_per_uc_;
        for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
          SigPair Y_rho_pair(tile_Y, tile_rho);
          for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
            const auto R_ord = tile_nu / ntiles_per_uc_;
            const auto R_3D = direct_3D_idx(R_ord, R_max_);

            const auto mR_3D = Vector3i(0, 0, 0) - R_3D;
            const auto RYmR_3D = RY_3D - R_3D;

            if (is_in_lattice_range(mR_3D, RD_max_) &&
                is_in_lattice_range(RYmR_3D, RJ_max_)) {
              const auto RYmR_ord = direct_ord_idx(RYmR_3D, RJ_max_);
              const auto mR_ord = direct_ord_idx(mR_3D, RD_max_);

              const auto shifted_Y = tile_Y_in_uc + RYmR_ord * ntiles_per_uc_;
              const auto shifted_rho = tile_rho + mR_ord * ntiles_per_uc_;
              const auto shifted_nu = tile_nu % ntiles_per_uc_;

              const auto val = in(shifted_Y, shifted_rho, shifted_nu);
              if (val > force_shape_threshold_) {
                Y_rho.insert(Y_rho_pair);
                break;
              }
            }
          }
        }
      }
    } else {
      for (auto tile_Y = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
        for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
          for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
            const auto val = in(tile_Y, tile_nu, tile_rho);
            if (val > force_shape_threshold_) {
              SigPair Y_rho_pair(tile_Y, tile_rho);
              Y_rho.insert(Y_rho_pair);
              break;
            }
          }
        }
      }
    }

    const auto &out_range = out_norms.range();
    TA::Tensor<float> out(out_range, 0.0);

    // F("X, rho, mu")
    for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
      for (auto const &Y_rho_pair : Y_rho) {
        const auto tile_Y = Y_rho_pair.first;
        const auto tile_rho = Y_rho_pair.second;
        out(tile_Y, tile_rho, tile_mu) = std::numeric_limits<float>::max();
      }
    }

    return out;
  }

  TA::Tensor<float> force_F_norms_new(TA::Tensor<float> const &in,
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

          const auto val = in(shifted_Y, shifted_nu, shifted_rho);
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

  array_type compute_eri2(madness::World &world,
                          const lcao::gaussian::Basis &bs0,
                          const lcao::gaussian::Basis &bs1) {
    const auto bs_ref_array = utility::make_array_of_refs(bs0, bs1);
    const auto bs_vector = lcao::gaussian::BasisVector{{bs0, bs1}};
    auto engine = lcao::gaussian::make_engine_pool(
        libint2::Operator::coulomb, bs_ref_array, libint2::BraKet::xs_xs);
    return lcao::gaussian::sparse_integrals(world, engine, bs_vector);
  }

  array_type compute_Q_ket(const array_type &C, const array_type &D,
                           const int64_t RJ_ord) {
    auto &world = C.world();

    array_type C_repl, D_repl;
    C_repl("Y, nu, sig") = C("Y, nu, sig");
    D_repl("rho, sig") = D("rho, sig");
    C_repl.make_replicated();
    D_repl.make_replicated();
    // must wait till all replicating is finished
    world.gop.fence();

    const auto &C_norms = C_repl.shape().data();

    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    const auto shifted_Y_latt_range = RJ_max_;

    const auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
    const auto vec_RJ = direct_vector(RJ_ord, RJ_max_, dcell_);
    const auto bs1 = shift_basis_origin(*obs_, vec_RJ);

    const auto trange = lcao::gaussian::detail::create_trange(
        lcao::gaussian::BasisVector{{*Y_dfbs_, *basisR_, *bs1}});
    const auto tvolume = trange.tiles_range().volume();
    const auto pmap = Policy::default_pmap(world, tvolume);

    // force Q norms
    const size_t ntiles_Y = Y_dfbs_->nclusters();
    const size_t ntiles_nu = basisR_->nclusters();
    const size_t ntiles_rho = bs1->nclusters();
    const size_t ntiles_sig = ntiles_per_uc_ * RD_size_;

    TA::Range range;
    range = TA::Range(std::array<size_t, 3>{{ntiles_Y, ntiles_nu, ntiles_rho}});
    TA::Tensor<float> norms(range, 0.0);

    using SigPair = std::pair<size_t, size_t>;
    using TilePairList = std::unordered_map<size_t, std::vector<SigPair>>;
    TilePairList significant_Y_rho;

    for (auto nu = 0ul; nu != ntiles_nu; ++nu) {
      auto R_ord = nu / ntiles_per_uc_;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_per_uc_;

      significant_Y_rho.insert(std::make_pair(nu, std::vector<SigPair>()));

      for (auto Y = 0ul; Y != ntiles_Y; ++Y) {
        auto RY_ord = Y / ntiles_per_uc_;
        auto RY_3D = direct_3D_idx(RY_ord, RY_max_);

        auto RYmR_3D = RY_3D - R_3D;
        if (!is_in_lattice_range(RYmR_3D, shifted_Y_latt_range)) continue;

        auto RYmR_ord = direct_ord_idx(RYmR_3D, shifted_Y_latt_range);
        auto Y_in_C = Y % ntiles_per_uc_ + RYmR_ord * ntiles_per_uc_;

        for (auto rho = 0ul; rho != ntiles_rho; ++rho) {
          for (auto sig = 0ul; sig != ntiles_sig; ++sig) {
            auto RD_ord = sig / ntiles_per_uc_;
            auto RJpRD_3D = RJ_3D + direct_3D_idx(RD_ord, RD_max_);

            if (!(RY_3D == R_3D) && !(RY_3D == RJpRD_3D)) continue;

            auto RJpRDmR_3D = RJpRD_3D - R_3D;
            if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

            auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
            if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) ==
                RJ_list_.end())
              continue;

            auto sig_in_C = sig % ntiles_per_uc_ + RJpRDmR_ord * ntiles_per_uc_;
            auto val = C_norms(Y_in_C, nu_in_C, sig_in_C);
            if (val >= force_shape_threshold_) {
              significant_Y_rho[nu].emplace_back(std::make_pair(Y, rho));
              norms(Y, nu, rho) = std::numeric_limits<float>::max();
              break;
            }
          }
        }
      }
    }

    TA::SparseShape<float> shape(world, norms, trange);

    array_type Q(world, trange, shape, pmap);

    auto create_task_Q_ket_tile = [&](Vector3i &RY_3D, Vector3i &R_3D,
                                      std::array<size_t, 3> external_tile_idx,
                                      std::array<size_t, 3> external_offset,
                                      std::array<size_t, 3> external_extent) {
      const auto Y_in_C = external_tile_idx[0];
      const auto nu_in_C = external_tile_idx[1];
      const auto rho = external_tile_idx[2];

      const auto offset_Y = external_offset[0];
      const auto offset_nu = external_offset[1];
      const auto offset_rho = external_offset[2];

      const auto ext_Y = external_extent[0];
      const auto ext_nu = external_extent[1];
      const auto ext_rho = external_extent[2];

      RowMatrixXd out_eig(ext_Y * ext_nu, ext_rho);
      out_eig.setZero();
      for (auto sig = 0ul; sig != ntiles_sig; ++sig) {
        auto RD_ord = sig / ntiles_per_uc_;
        auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
        auto RJpRD_3D = RJ_3D + RD_3D;

        if (!(RY_3D == R_3D) && !(RY_3D == RJpRD_3D)) continue;

        auto RJpRDmR_3D = RJpRD_3D - R_3D;
        if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

        auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
        if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) ==
            RJ_list_.end())
          continue;

        auto sig_in_C = sig % ntiles_per_uc_ + RJpRDmR_ord * ntiles_per_uc_;

        std::array<size_t, 3> idx_C = {{Y_in_C, nu_in_C, sig_in_C}};
        std::array<size_t, 2> idx_D = {{rho, sig}};
        if (D_repl.is_zero(idx_D) || C_repl.is_zero(idx_C)) continue;

        Tile C = C_repl.find(idx_C);
        Tile D = D_repl.find(idx_D);
        const auto C_ext = C.range().extent();
        const auto D_ext = D.range().extent();
        assert(C_ext[2] == D_ext[1]);
        RowMatrixXd C_eig = TA::eigen_map(C, C_ext[0] * C_ext[1], C_ext[2]);
        RowMatrixXd D_eig = TA::eigen_map(D, D_ext[0], D_ext[1]);
        RowMatrixXd D_trans = D_eig.transpose();
        out_eig += C_eig * D_trans;
      }

      TA::Range out_range;
      std::array<size_t, 3> lb = {{offset_Y, offset_nu, offset_rho}};
      std::array<size_t, 3> ub = {
          {offset_Y + ext_Y, offset_nu + ext_nu, offset_rho + ext_rho}};
      out_range = TA::Range(lb, ub);

      TA::TensorD out_tile(out_range, 0.0);
      TA::eigen_map(out_tile, ext_Y * ext_nu, ext_rho) = out_eig;

      return out_tile;
    };

    for (auto nu = 0ul; nu != ntiles_nu; ++nu) {
      auto R_ord = nu / ntiles_per_uc_;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_per_uc_;

      const auto &rng_nu = trange.dim(1).tile(nu);
      const auto offset_nu = rng_nu.first;
      const auto ext_nu = rng_nu.second - rng_nu.first;

      for (const auto &significant_pair : significant_Y_rho[nu]) {
        auto Y = significant_pair.first;
        auto rho = significant_pair.second;

        std::array<size_t, 3> tile_idx = {{Y, nu, rho}};
        const size_t tile_ord = trange.tiles_range().ordinal(tile_idx);

        if (Q.is_local(tile_idx) && !Q.is_zero(tile_idx)) {
          auto RY_ord = Y / ntiles_per_uc_;
          auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
          auto Y_in_C = Y % ntiles_per_uc_ +
                        direct_ord_idx(RY_3D - R_3D, shifted_Y_latt_range) *
                            ntiles_per_uc_;

          const auto &rng_Y = trange.dim(0).tile(Y);
          const auto offset_Y = rng_Y.first;
          const auto ext_Y = rng_Y.second - rng_Y.first;

          const auto &rng_rho = trange.dim(2).tile(rho);
          const auto offset_rho = rng_rho.first;
          const auto ext_rho = rng_rho.second - rng_rho.first;

          auto fut_tile = world.taskq.add(
              create_task_Q_ket_tile, RY_3D, R_3D,
              std::array<size_t, 3>{{Y_in_C, nu_in_C, rho}},
              std::array<size_t, 3>{{offset_Y, offset_nu, offset_rho}},
              std::array<size_t, 3>{{ext_Y, ext_nu, ext_rho}});

          Q.set(tile_ord, fut_tile);
        }
      }
    }

    world.gop.fence();
    Q.truncate();

    return Q;
  }

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

    for (auto tile_Y = 0ul, task = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
      const auto RYmR_ord = tile_Y / ntiles_per_uc_;
      const auto RYmR_3D = direct_3D_idx(RYmR_ord, RYmR_max_);
      if (!is_in_lattice_range(RYmR_3D, RJ_max_)) continue;
      const auto RY_ord_in_C = direct_ord_idx(RYmR_3D, RJ_max_);
      const auto Y_in_C = tile_Y % ntiles_per_uc_ + RY_ord_in_C * ntiles_per_uc_;

      for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {

        for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
          const auto RJmR_ord = tile_rho / ntiles_per_uc_;
          const auto RJmR_3D = direct_3D_idx(RJmR_ord, RJmR_max_);
          const auto rho_in_D = tile_rho % ntiles_per_uc_;

          for (auto tile_sig = 0ul; tile_sig != ntiles_sig; ++tile_sig, ++task) {

            if (task % nproc == me) {
              const auto RD_ord = tile_sig / ntiles_per_uc_;
              const auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
              const auto RJmRpRD_3D = RJmR_3D + RD_3D;
              if (!(RYmR_3D == Vector3i({0, 0, 0})) && !(RYmR_3D == RJmRpRD_3D)) continue;

              if (!is_in_lattice_range(RJmRpRD_3D, RJ_max_)) continue;
              const auto RJmRpRD_ord = direct_ord_idx(RJmRpRD_3D, RJ_max_);

              if (std::find(RJ_list_.begin(), RJ_list_.end(), RJmRpRD_ord) == RJ_list_.end()) continue;

              const auto sig_in_C = tile_sig % ntiles_per_uc_ + RJmRpRD_ord * ntiles_per_uc_;

              // get future of tiles C and D
              std::array<size_t, 3> idx_C = {{Y_in_C, tile_nu, sig_in_C}};
              std::array<size_t, 2> idx_D = {{rho_in_D, tile_sig}};
              if (C_repl.is_zero(idx_C) || D_repl.is_zero(idx_D)) continue;

              auto C_tile = C_repl.find(idx_C);
              auto D_tile = D_repl.find(idx_D);
              WorldObject_::task(me, &PeriodicCADFKBuilder_::compute_Q_ket_task, C_tile, D_tile, std::array<size_t, 3>{{tile_Y, tile_nu, tile_rho}});
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

  void compute_Q_ket_task(Tile C_tile, Tile D_tile, std::array<size_t, 3> tile_idx) {
    const auto ext_C = C_tile.range().extent();
    const auto ext_D = D_tile.range().extent();
    assert(ext_C[2] == ext_D[1]);

    RowMatrixXd C_eig = TA::eigen_map(C_tile, ext_C[0] * ext_C[1], ext_C[2]);
    RowMatrixXd D_eig = TA::eigen_map(D_tile, ext_D[0], ext_D[1]);
    RowMatrixXd result_eig = C_eig * D_eig.transpose();

    const auto tile_Y = tile_idx[0];
    const auto tile_nu = tile_idx[1];
    const auto tile_rho = tile_idx[2];

    const auto &tr1_Y = Q_trange_.dim(0);
    const auto &tr1_nu = Q_trange_.dim(1);
    const auto &tr1_rho = Q_trange_.dim(2);

    const auto &rng_Y = tr1_Y.tile(tile_Y);
    const auto &rng_nu = tr1_nu.tile(tile_nu);
    const auto &rng_rho = tr1_rho.tile(tile_rho);

    const auto rng_Y_size = rng_Y.second - rng_Y.first;
    assert(rng_Y_size == ext_C[0]);

    const auto result_rng = TA::Range({rng_Y, rng_nu, rng_rho});
    Tile result_tile(result_rng, 0.0);
    TA::eigen_map(result_tile, ext_C[0] * ext_C[1], ext_D[0]) = result_eig;

    const auto ntiles_nu = Q_bs_nu_->nclusters();
    const auto ntiles_rho = Q_bs_rho_->nclusters();
    const auto ord = tile_Y * ntiles_nu * ntiles_rho + tile_nu * ntiles_rho + tile_rho;

    PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
  }

  array_type compute_contr_EQ(array_type const &Q, int64_t const RJ,
                              double target_precision) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();
    target_precision_ = target_precision;

    // # of tiles per basis
    const auto ntiles_nu = basisR_->nclusters();
    const auto ntiles_X = X_dfbs_->nclusters();
    const auto ntiles_sig = basisRD_->nclusters();

    auto empty = TA::Future<Tile>(Tile());
    for (auto tile_X = 0ul, task = 0ul; tile_X != ntiles_X; ++tile_X) {
      for (auto tile_sig = 0ul; tile_sig != ntiles_sig; ++tile_sig, ++task) {
        if (task % nproc == me) {
          auto vec_fut_tile =
              madness::future_vector_factory<Tile>(ntiles_per_uc_);
          for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
            std::array<size_t, 3> idx_Q = {{tile_X, tile_sig, tile_mu}};
            vec_fut_tile[tile_mu] = Q.is_zero(idx_Q) ? empty : Q.find(idx_Q);
          }

          WorldObject_::task(me, &PeriodicCADFKBuilder_::compute_contr_EQ_task,
                             vec_fut_tile, RJ,
                             std::array<size_t, 2>{{tile_X, tile_sig}});
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

  void compute_contr_EQ_task(std::vector<madness::Future<Tile>> Q_mu_vec,
                             int64_t RJ, std::array<size_t, 2> tile_idx) {
    const auto tile_X = tile_idx[0];
    const auto tile_sig = tile_idx[1];

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;
    const auto RJ_3D = direct_3D_idx(RJ, RJ_max_);

    const auto RX_ord = tile_X / ntiles_per_uc_;
    const auto RX_3D = direct_3D_idx(RX_ord, RJ_max_);
    const auto tile_X_in_uc = tile_X % ntiles_per_uc_;

    const auto RD_ord = tile_sig / ntiles_per_uc_;
    const auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
    const auto tile_sig_in_uc = tile_sig % ntiles_per_uc_;

    const auto ntiles_nu = basisR_->nclusters();

    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);

    // get reference to basis sets
    const auto &eri3_basis_X = eri3_X_dfbs_;
    const auto &eri3_basis_nu = eri3_bs0_;
    const auto &eri3_basis_sig = eri3_bs1_;
    const auto &basis_sig = basisRD_;
    const auto &basis_X = X_dfbs_;

    // initialize result tiles with ranges&zeros and grab pointers to Q tiles
    std::vector<Tile> result_tiles(ntiles_per_uc_ * ntiles_nu, Tile());
    std::vector<typename Tile::numeric_type *> Q_ptrs(ntiles_per_uc_);
    auto num_nonzero_Q = 0;
    for (auto tile_mu = 0ul, ord = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
      const auto &rng_mu = tr0.tile(tile_mu);
      Q_ptrs[tile_mu] = Q_mu_vec[tile_mu].get().data();
      if (Q_ptrs[tile_mu] != nullptr) num_nonzero_Q++;

      for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu, ++ord) {
        const auto &rng_nu = tr1.tile(tile_nu);
        // 2-d tile ranges describing the contribution blocks
        // produced by this
        auto result_rng = TA::Range({rng_mu, rng_nu});
        // initialize contribution to the result matrices
        auto &result_tile = result_tiles[ord];
        result_tile = Tile(std::move(result_rng), 0.0);
      }
    }

    if (num_nonzero_Q > 0) {
      auto engine = engines_->local();

      // compute max value of mu for each (X, sigma) pair for Q
      const auto &cluster_X = basis_X->cluster_shells()[tile_X];
      const auto &cluster_sig = basis_sig->cluster_shells()[tile_sig];
      const auto nshells_X = cluster_X.size();
      shellpair_val_list_t Xsig_shpair_val_list;
      for (auto sh_X = 0; sh_X != nshells_X; ++sh_X) {
        Xsig_shpair_val_list.insert(
            std::make_pair(sh_X, std::vector<std::pair<size_t, double>>()));
      }
      for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
        if (Q_ptrs[tile_mu] != nullptr) {
          add_shell_pair_with_val(Xsig_shpair_val_list, Q_mu_vec[tile_mu].get(),
                                  cluster_X, cluster_sig,
                                  shell_pair_threshold_);
        }
      }

      // test
      //      for (auto sh_X = 0; sh_X != nshells_X; ++sh_X) {
      //        std::cout << "sh_X = " << sh_X << "\n";
      //        for (auto sh_sig_with_val : Xsig_shpair_val_list[sh_X]) {
      //          auto sh_sig = sh_sig_with_val.first;
      //          auto val = sh_sig_with_val.second;
      //          std::cout << "  sh_sig = " << sh_sig << ", max = " << val <<
      //          "\n";
      //        }
      //      }

      for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
        const auto R_ord = tile_nu / ntiles_per_uc_;
        const auto R_3D = direct_3D_idx(R_ord, R_max_);

        const auto RJpRDmR_3D = RJ_3D + RD_3D - R_3D;
        if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

        auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
        if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) ==
            RJ_list_.end())
          continue;

        const auto RXmR_3D = RX_3D - R_3D;
        const auto RXmR_ord = direct_ord_idx(RXmR_3D, RJ_max_ + R_max_);

        const auto eri3_tile_X = tile_X_in_uc + RXmR_ord * ntiles_per_uc_;
        const auto eri3_tile_nu = tile_nu % ntiles_per_uc_;
        const auto eri3_tile_sig =
            tile_sig_in_uc + RJpRDmR_ord * ntiles_per_uc_;

        //        ExEnv::out0() << "\ntile_X = " << tile_X << ", tile_sig = " <<
        //        tile_sig << std::endl; ExEnv::out0() << "tile_nu = " <<
        //        tile_nu << std::endl; ExEnv::out0() << "eri3_tile_X = " <<
        //        eri3_tile_X
        //                      << ", eri3_tile_sig = " << eri3_tile_sig
        //                      << ", eri3_tile_nu = " << eri3_tile_nu
        //                      << std::endl;

        // shell clusters for this tile
        const auto &eri3_cluster_X =
            eri3_basis_X->cluster_shells()[eri3_tile_X];
        const auto &eri3_cluster_nu =
            eri3_basis_nu->cluster_shells()[eri3_tile_nu];
        const auto &eri3_cluster_sig =
            eri3_basis_sig->cluster_shells()[eri3_tile_sig];

        // # of shells in each cluster
        const auto eri3_nshells_X = eri3_cluster_X.size();
        const auto eri3_nshells_nu = eri3_cluster_nu.size();
        const auto eri3_nshells_sig = eri3_cluster_sig.size();

        // 1-d tile ranges
        const auto &rng_sig = basisRD_trange1_.tile(tile_sig);
        const auto &eri3_rng_X = eri3_X_trange1_.tile(eri3_tile_X);
        const auto &eri3_rng_nu = eri3_bs0_trange1_.tile(eri3_tile_nu);
        const auto &eri3_rng_sig = eri3_bs1_trange1_.tile(eri3_tile_sig);

        // range sizes
        const auto &rng_nu = tr1.tile(tile_nu);
        const auto rng_size_nu = rng_nu.second - rng_nu.first;
        const auto rng_size_sig = rng_sig.second - rng_sig.first;
        // test
        //        const auto rng_size_X = eri3_rng_X.second - eri3_rng_X.first;
        //        RowMatrixXd eri3_eig(rng_size_X * rng_size_sig, rng_size_nu);
        //        eri3_eig.setZero();

        // compute eri3 * Q contribution to all Fock matrices
        {
          // index of first shell in this cluster
          const auto eri3_sh_offset_nu =
              eri3_bs0_shell_offset_map_[eri3_tile_nu];
          const auto eri3_sh_offset_sig =
              eri3_bs1_shell_offset_map_[eri3_tile_sig];

          // index of last shell in this cluster
          const auto eri3_sh_max_nu = eri3_sh_offset_nu + eri3_nshells_nu;
          const auto eri3_sh_max_sig = eri3_sh_offset_sig + eri3_nshells_sig;

          // determine if this task is worth computing by checking whether there
          // are significant shell pairs in (tile_nu, tile_sig)
          auto is_significant = false;
          {
            auto eri3_sh_nu = eri3_sh_offset_nu;
            for (; eri3_sh_nu != eri3_sh_max_nu; ++eri3_sh_nu) {
              for (const auto &eri3_sh_sig : sig_shellpair_list_[eri3_sh_nu]) {
                if (eri3_sh_sig >= eri3_sh_offset_sig &&
                    eri3_sh_sig < eri3_sh_max_sig) {
                  is_significant = true;
                  break;
                }
              }
              if (is_significant) break;
            }
          }

          if (is_significant) {
            auto &screen = *(p_screener_);
            const auto engine_precision = target_precision_;
            engine.set_precision(engine_precision);
            const auto &computed_shell_sets = engine.results();

            // compute offset list of cluster_sig
            auto eri3_offset_list_sig =
                compute_func_offset_list(eri3_cluster_sig, eri3_rng_sig.first);
            auto eri3_offset_list_X =
                compute_func_offset_list(eri3_cluster_X, eri3_rng_X.first);
            auto eri3_offset_list_nu =
                compute_func_offset_list(eri3_cluster_nu, eri3_rng_nu.first);

            auto cf_offset_X = 0;
            auto eri3_bf_offset_X = eri3_rng_X.first;
            size_t cf_offset_sig, eri3_bf_offset_sig;
            size_t cf_offset_nu, eri3_bf_offset_nu;

            for (auto eri3_sh_X = 0; eri3_sh_X != eri3_nshells_X; ++eri3_sh_X) {
              const auto &eri3_shell_X = eri3_cluster_X[eri3_sh_X];
              const auto eri3_nf_X = eri3_shell_X.size();

              for (const auto &sig_with_maxval :
                   Xsig_shpair_val_list[eri3_sh_X]) {
                const auto eri3_sh_sig = sig_with_maxval.first;
                const auto max_Q = sig_with_maxval.second;

                const auto &eri3_shell_sig = eri3_cluster_sig[eri3_sh_sig];
                const auto eri3_nf_sig = eri3_shell_sig.size();
                std::tie(cf_offset_sig, eri3_bf_offset_sig) =
                    eri3_offset_list_sig[eri3_sh_sig];

                const auto sh_sig_in_basis = eri3_sh_sig + eri3_sh_offset_sig;
                for (const auto &sh_nu_in_basis :
                     nu_shellpair_list_[sh_sig_in_basis]) {
                  if (sh_nu_in_basis < eri3_sh_offset_nu ||
                      sh_nu_in_basis >= eri3_sh_max_nu)
                    continue;

                  const auto eri3_sh_nu = sh_nu_in_basis - eri3_sh_offset_nu;
                  std::tie(cf_offset_nu, eri3_bf_offset_nu) =
                      eri3_offset_list_nu[eri3_sh_nu];

                  const auto &eri3_shell_nu = eri3_cluster_nu[eri3_sh_nu];
                  const auto eri3_nf_nu = eri3_shell_nu.size();

                  if (screen.skip(eri3_bf_offset_X, eri3_bf_offset_nu,
                                  eri3_bf_offset_sig, max_Q))
                    continue;

                  //                  std::cout << "eri3_sh_X = " << eri3_sh_X
                  //                            << ", eri3_sh_sig = " <<
                  //                            eri3_sh_sig
                  //                            << ", eri3_sh_nu = " <<
                  //                            eri3_sh_nu
                  //                            << std::endl;

                  num_ints_computed_ += eri3_nf_X * eri3_nf_nu * eri3_nf_sig;

                  // compute shell set
                  engine.compute(eri3_shell_X, eri3_shell_sig, eri3_shell_nu);
                  const auto &eri3 = computed_shell_sets[0];

                  if (eri3 != nullptr) {
                    for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_;
                         ++tile_mu) {
                      const auto Q_ptr = Q_ptrs[tile_mu];
                      if (Q_ptr == nullptr) continue;

                      const auto &rng_mu = tr0.tile(tile_mu);
                      const auto rng_size_mu = rng_mu.second - rng_mu.first;
                      const auto QX_stride = rng_size_sig * rng_size_mu;

                      auto &result_tile =
                          result_tiles[tile_mu * ntiles_nu + tile_nu];
                      // grab ptrs to tile data to make addressing more
                      // efficient
                      auto *result_ptr = result_tile.data();

                      for (auto f_X = 0, eri3_ord = 0; f_X != eri3_nf_X;
                           ++f_X) {
                        const auto cf_X = f_X + cf_offset_X;
                        const auto Qsig_offset = cf_X * QX_stride;
                        for (auto f_sig = 0; f_sig != eri3_nf_sig; ++f_sig) {
                          const auto cf_sig = f_sig + cf_offset_sig;
                          const auto Qmu_offset =
                              Qsig_offset + cf_sig * rng_size_mu;
                          for (auto f_nu = 0; f_nu != eri3_nf_nu;
                               ++f_nu, ++eri3_ord) {
                            const auto cf_nu = f_nu + cf_offset_nu;
                            const auto eri3_value = eri3[eri3_ord];
                            //                              eri3_eig(cf_X *
                            //                              rng_size_sig +
                            //                              cf_sig, cf_nu) =
                            //                              eri3_value;
                            for (auto cf_mu = 0; cf_mu != rng_size_mu;
                                 ++cf_mu) {
                              const auto cf_mu_nu = cf_mu * rng_size_nu + cf_nu;
                              const auto cf_X_sig_mu = Qmu_offset + cf_mu;
                              result_ptr[cf_mu_nu] +=
                                  eri3_value * Q_ptr[cf_X_sig_mu];
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }

              cf_offset_X += eri3_nf_X;
              eri3_bf_offset_X += eri3_nf_X;
            }

            //            std::cout << "\neri3 matrix = \n" << eri3_eig << "\n";
          }
        }
      }

      // accumulate the local contributions
      for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
        const auto Q_ptr = Q_ptrs[tile_mu];
        if (Q_ptr == nullptr) continue;

        const auto tile_nu_offset = tile_mu * ntiles_nu;
        for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
          const auto tile_ord = tile_nu_offset + tile_nu;
          auto &result_tile = result_tiles[tile_ord];
          if (result_tile.data() != nullptr) {
            PeriodicCADFKBuilder_::accumulate_local_task(result_tile, tile_ord);
          }
        }
      }
    }
  }

  array_type compute_contr_FQ(const array_type &F, const array_type &Q) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();

    // # of tiles per basis
    const auto ntiles_Y = Y_dfbs_->nclusters();
    const auto ntiles_nu = basisR_->nclusters();

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    auto empty = TA::Future<Tile>(Tile());
    for (auto tile_Y = 0ul, task = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
      const auto RY_ord = tile_Y / ntiles_per_uc_;
      const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
      const auto tile_Y_in_uc = tile_Y % ntiles_per_uc_;

      for (auto tile_rho = 0ul; tile_rho != ntiles_per_uc_; ++tile_rho) {
        if (task % nproc == me) {
          auto F_mu_vec = madness::future_vector_factory<Tile>(ntiles_per_uc_);
          for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
            std::array<size_t, 3> idx_F = {{tile_Y, tile_rho, tile_mu}};
            F_mu_vec[tile_mu] = F.is_zero(idx_F) ? empty : F.find(idx_F);
          }

          auto Q_nu_vec = madness::future_vector_factory<Tile>(ntiles_nu);
          for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
            const auto R_ord = tile_nu / ntiles_per_uc_;
            const auto R_3D = direct_3D_idx(R_ord, R_max_);

            const auto mR_3D = Vector3i{0, 0, 0} - R_3D;
            const auto RYmR_3D = RY_3D - R_3D;

            if (!is_in_lattice_range(mR_3D, RD_max_) ||
                !is_in_lattice_range(RYmR_3D, RJ_max_)) {
              Q_nu_vec[tile_nu] = empty;
            } else {
              const auto RYmR_ord = direct_ord_idx(RYmR_3D, RJ_max_);
              const auto mR_ord = direct_ord_idx(mR_3D, RD_max_);

              const auto shifted_Y = tile_Y_in_uc + RYmR_ord * ntiles_per_uc_;
              const auto shifted_rho = tile_rho + mR_ord * ntiles_per_uc_;
              const auto shifted_nu = tile_nu % ntiles_per_uc_;
              std::array<size_t, 3> idx_Q = {
                  {shifted_Y, shifted_rho, shifted_nu}};
              Q_nu_vec[tile_nu] = Q.is_zero(idx_Q) ? empty : Q.find(idx_Q);
            }
          }

          WorldObject_::task(me, &PeriodicCADFKBuilder_::compute_contr_FQ_task,
                             F_mu_vec, Q_nu_vec);
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

  void compute_contr_FQ_task(std::vector<madness::Future<Tile>> F_mu_vec,
                             std::vector<madness::Future<Tile>> Q_nu_vec) {
    const auto ntiles_nu = basisR_->nclusters();
    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);

    std::vector<RowMatrixXd> F_eigs(ntiles_per_uc_);
    std::vector<typename Tile::numeric_type *> F_ptrs(ntiles_per_uc_);
    for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
      auto &F_tile = F_mu_vec[tile_mu].get();
      F_ptrs[tile_mu] = F_tile.data();
      if (F_ptrs[tile_mu] != nullptr) {
        const auto ext = F_tile.range().extent_data();
        F_eigs[tile_mu] = TA::eigen_map(F_tile, ext[0] * ext[1], ext[2]);
      }
    }

    for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
      auto &Q_tile = Q_nu_vec[tile_nu].get();
      if (Q_tile.data() != nullptr) {
        const auto ext_Q = Q_tile.range().extent_data();
        RowMatrixXd Q_eig =
            TA::eigen_map(Q_tile, ext_Q[0] * ext_Q[1], ext_Q[2]);
        const auto &rng_nu = tr1.tile(tile_nu);

        for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
          if (F_ptrs[tile_mu] != nullptr) {
            const auto ext_F = F_mu_vec[tile_mu].get().range().extent_data();
            const auto &F_eig = F_eigs[tile_mu];

            RowMatrixXd result_eig = F_eig.transpose() * Q_eig;

            const auto &rng_mu = tr0.tile(tile_mu);
            const auto result_rng = TA::Range({rng_mu, rng_nu});

            Tile result_tile(result_rng, 0.0);
            TA::eigen_map(result_tile, ext_F[2], ext_Q[2]) = result_eig;

            const auto ord = tile_mu * ntiles_nu + tile_nu;
            PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
          }
        }
      }
    }
  }

  void compute_contr_FQ_task_v2(Tile F_tile,
                                Tile Q_tile, std::array<size_t, 2> tile_idx) {
    const auto ext_F = F_tile.range().extent();
    const auto ext_Q = Q_tile.range().extent();
    assert(ext_F[0] == ext_Q[0]);
    assert(ext_F[1] == ext_Q[1]);

    RowMatrixXd F_eig = TA::eigen_map(F_tile, ext_F[0] * ext_F[1], ext_F[2]);
    RowMatrixXd Q_eig = TA::eigen_map(Q_tile, ext_Q[0] * ext_Q[1], ext_Q[2]);
    RowMatrixXd result_eig = F_eig.transpose() * Q_eig;

    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);

    const auto tile_mu = tile_idx[0];
    const auto tile_nu = tile_idx[1];

    const auto &rng_mu = tr0.tile(tile_mu);
    const auto &rng_nu = tr1.tile(tile_nu);
    const auto result_rng = TA::Range({rng_mu, rng_nu});

    Tile result_tile(result_rng, 0.0);
    TA::eigen_map(result_tile, ext_F[2], ext_Q[2]) = result_eig;

    const auto ntiles_nu = basisR_->nclusters();
    const auto ord = tile_mu * ntiles_nu + tile_nu;

    PeriodicCADFKBuilder_::accumulate_local_task(result_tile, ord);
   }

  array_type compute_contr_FQ_v2(const array_type &F, const array_type &Q) {
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

    auto empty = TA::Future<Tile>(Tile());
    for (auto tile_Y = 0ul, task = 0ul; tile_Y != ntiles_Y; ++tile_Y) {
      const auto RY_ord = tile_Y / ntiles_per_uc_;
      const auto RY_3D = direct_3D_idx(RY_ord, RY_max_);
      const auto tile_Y_in_uc = tile_Y % ntiles_per_uc_;

      for (auto tile_rho = 0ul; tile_rho != ntiles_rho; ++tile_rho) {
        const auto RJ_ord = tile_rho / ntiles_per_uc_;
        const auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
        const auto tile_rho_in_uc = tile_rho % ntiles_per_uc_;

        if (task % nproc == me) {
          auto F_mu_vec = madness::future_vector_factory<Tile>(ntiles_mu);
          for (auto tile_mu = 0ul; tile_mu != ntiles_mu; ++tile_mu) {
            std::array<size_t, 3> idx_F = {{tile_Y, tile_rho, tile_mu}};
            F_mu_vec[tile_mu] = F.is_zero(idx_F) ? empty : F.find(idx_F);
          }

          auto Q_nu_vec = madness::future_vector_factory<Tile>(ntiles_nu);
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
            std::array<size_t, 3> idx_Q = {
                {shifted_Y, shifted_rho, shifted_nu}};
            Q_nu_vec[tile_nu] = Q.is_zero(idx_Q) ? empty : Q.find(idx_Q);
          }

          WorldObject_::task(me, &PeriodicCADFKBuilder_::compute_contr_FQ_task,
                             F_mu_vec, Q_nu_vec);
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

  array_type compute_contr_FQ_v3(const array_type &F, const array_type &Q) {
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

    auto empty = TA::Future<Tile>(Tile());
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
              WorldObject_::task(me, &PeriodicCADFKBuilder_::compute_contr_FQ_task_v2, F_tile, Q_tile, std::array<size_t, 2>{{tile_mu, tile_nu}});
            }

          }
        }


//        if (task % nproc == me) {
//          auto F_mu_vec = madness::future_vector_factory<Tile>(ntiles_mu);
//          for (auto tile_mu = 0ul; tile_mu != ntiles_mu; ++tile_mu) {
//            std::array<size_t, 3> idx_F = {{tile_Y, tile_rho, tile_mu}};
//            F_mu_vec[tile_mu] = F.is_zero(idx_F) ? empty : F.find(idx_F);
//          }

//          auto Q_nu_vec = madness::future_vector_factory<Tile>(ntiles_nu);
//          for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu) {
//            const auto R_ord = tile_nu / ntiles_per_uc_;
//            const auto R_3D = direct_3D_idx(R_ord, R_max_);

//            const auto RYmR_3D = RY_3D - R_3D;
//            const auto RJmR_3D = RJ_3D - R_3D;

//            const auto RYmR_ord = direct_ord_idx(RYmR_3D, RYmR_max_);
//            const auto RJmR_ord = direct_ord_idx(RJmR_3D, RJmR_max_);

//            const auto shifted_Y = tile_Y_in_uc + RYmR_ord * ntiles_per_uc_;
//            const auto shifted_rho = tile_rho_in_uc + RJmR_ord * ntiles_per_uc_;
//            const auto shifted_nu = tile_nu % ntiles_per_uc_;
//            std::array<size_t, 3> idx_Q = {
//                {shifted_Y, shifted_rho, shifted_nu}};
//            Q_nu_vec[tile_nu] = Q.is_zero(idx_Q) ? empty : Q.find(idx_Q);
//          }

//          WorldObject_::task(me, &PeriodicCADFKBuilder_::compute_contr_FQ_task,
//                             F_mu_vec, Q_nu_vec);
//        }
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
   * \brief This computes significant (X, sigma) shell pair list
   * of a Q(X, mu, sigma) tile, and grab norm of all nu's
   * for a specific (X, sigma) pair.
   *
   * \param arg_tile Q tile
   * \param cluster_X a shell cluster (a.k.a. std::vector<Shell>)
   * \param cluster_sig a shell cluster (a.k.a. std::vector<Shell>)
   * \param threshold
   * \return a list of (X, sigma) pairs with the norm
   */
  shellpair_val_list_t compute_shell_pair_with_val(
      const Tile &arg_tile, const ShellVec &cluster_X,
      const ShellVec &cluster_sig, const double threshold = 1e-12) {
    // # of shells in this shell cluster
    const auto nsh_X = cluster_X.size();
    const auto nsh_sig = cluster_sig.size();
    const auto ext = arg_tile.range().extent_data();
    RowMatrixXd arg_eig = TA::eigen_map(arg_tile, ext[0] * ext[1], ext[2]);

    shellpair_val_list_t shpair_maxval_list;
    auto f_sig = 0;
    for (auto sh_sig = 0; sh_sig != nsh_sig; ++sh_sig) {
      const auto nf_sh_sig = cluster_sig[sh_sig].size();
      shpair_maxval_list.insert(
          std::make_pair(sh_sig, std::vector<std::pair<size_t, double>>()));

      auto f_X = 0;
      for (auto sh_X = 0; sh_X != nsh_X; ++sh_X) {
        const auto nf_sh_X = cluster_X[sh_X].size();
        // TODO try if maxCoeff is better in screening/accuracy
        const auto val =
            arg_eig.block(f_X * ext[1], f_sig, nf_sh_X * ext[1], nf_sh_sig)
                .template lpNorm<Eigen::Infinity>();

        if (val >= threshold)
          shpair_maxval_list[sh_sig].emplace_back(std::make_pair(sh_X, val));

        f_X += nf_sh_X;
      }

      f_sig += nf_sh_sig;
    }

    // resort shell list in increasing order
    // note that the customized comparator only compares the first element
    // because comparing the second element is meaningless
    for (auto sh_sig = 0; sh_sig != nsh_sig; ++sh_sig) {
      auto &list = shpair_maxval_list[sh_sig];
      std::sort(list.begin(), list.end(),
                [](auto &l, auto &r) { return l.first < r.first; });
    }

    return shpair_maxval_list;
  }

  void add_shell_pair_with_val(shellpair_val_list_t &shellpair_val_list,
                               const Tile &arg_tile, const ShellVec &cluster_X,
                               const ShellVec &cluster_sig,
                               const double threshold = 1e-12) {
    const auto ext = arg_tile.range().extent_data();
    const auto sig_size = ext[1];
    const auto mu_size = ext[2];
    const auto Xsig_size = ext[0] * sig_size;

    std::vector<double> maxes(Xsig_size, 0.0);
    double const *ptr = arg_tile.data();
    for (auto Xsig = 0u; Xsig < Xsig_size; ++Xsig) {
      for (auto mu = 0u; mu < mu_size; ++mu, ++ptr) {
        const auto val = maxes[Xsig];
        maxes[Xsig] = std::max(val, std::abs(*ptr));
      }
    }

    // # of shells in this shell cluster
    const auto nsh_X = cluster_X.size();
    const auto nsh_sig = cluster_sig.size();

    auto f_X_lb = 0;
    const auto iter0 = maxes.begin();
    for (auto sh_X = 0; sh_X != nsh_X; ++sh_X) {
      const auto nf_X = cluster_X[sh_X].size();
      const auto f_X_ub = f_X_lb + nf_X;
      auto &list = shellpair_val_list[sh_X];

      auto f_sig_lb = 0;
      for (auto sh_sig = 0; sh_sig != nsh_sig; ++sh_sig) {
        const auto nf_sig = cluster_sig[sh_sig].size();

        auto iter_old =
            std::find_if(list.begin(), list.end(),
                         [sh_sig](const std::pair<size_t, double> &pair) {
                           return pair.first == sh_sig;
                         });
        const auto iter_sig = iter0 + f_sig_lb;
        auto maxval = 0.0;
        for (auto fs_X = f_X_lb; fs_X != f_X_ub; ++fs_X) {
          const auto lower_iter = iter_sig + fs_X * sig_size;
          const auto max_tmp_it =
              std::max_element(lower_iter, lower_iter + nf_sig);
          maxval = std::max(maxval, *max_tmp_it);
        }

        if (maxval >= threshold) {
          if (iter_old == list.end()) {
            list.emplace_back(std::make_pair(sh_sig, maxval));
          } else {
            iter_old->second = std::max(maxval, iter_old->second);
          }
        }

        f_sig_lb += nf_sig;
      }

      f_X_lb += nf_X;
    }

    // resort shell list in increasing order
    // note that the customized comparator only compares the first element
    // because comparing the second element is meaningless
    for (auto sh_X = 0; sh_X != nsh_X; ++sh_X) {
      auto &list = shellpair_val_list[sh_X];
      std::sort(list.begin(), list.end(),
                [](const auto &l, const auto &r) { return l.first < r.first; });
    }
  }

  /*!
   * \brief This computes significant (X, sigma) shell pair list
   * of a Q(X, mu, sigma) tile, and grab norm of all nu's
   * for a specific (X, sigma) pair.
   *
   * \param arg_tile Q tile
   * \param cluster_X a shell cluster (a.k.a. std::vector<Shell>)
   * \param cluster_sig a shell cluster (a.k.a. std::vector<Shell>)
   * \param threshold
   * \return a list of (X, sigma) pairs with the norm
   */
  shellpair_val_list_t compute_shell_pair_with_val(
      const array_type &arg_array, const double threshold = 1e-12) {
    auto &world = arg_array.world();
    std::mutex mx;
    shellpair_val_list_t shpair_maxval_list;

    const auto ntiles_X = X_dfbs_->nclusters();
    const auto ntiles_sig = basisRD_->nclusters();

    auto max_abs_coeff = [&mx](Tile &arg_tile, double &max_val,
                               std::array<size_t, 2> func_idx,
                               std::array<size_t, 2> nfs) {
      const auto ext = arg_tile.range().extent_data();
      RowMatrixXd arg_eig = TA::eigen_map(arg_tile, ext[0] * ext[1], ext[2]);
      const auto f_sig = func_idx[0];
      const auto f_X = func_idx[1];
      const auto nf_sh_sig = nfs[0];
      const auto nf_sh_X = nfs[1];
      const auto val =
          arg_eig.block(f_X * ext[1], f_sig, nf_sh_X * ext[1], nf_sh_sig)
              .template lpNorm<Eigen::Infinity>();
      mx.lock();
      max_val = std::max(max_val, val);
      mx.unlock();
    };

    auto sh_offset_sig = 0ul;
    for (auto tile_sig = 0ul; tile_sig != ntiles_sig; ++tile_sig) {
      const auto &cluster_sig = basisRD_->cluster_shells()[tile_sig];
      const auto nsh_sig = cluster_sig.size();

      auto f_sig = 0ul;
      for (auto sh_sig = 0ul; sh_sig != nsh_sig; ++sh_sig) {
        const auto sh_sig_in_basis = sh_sig + sh_offset_sig;
        const auto nf_sh_sig = cluster_sig[sh_sig].size();
        shpair_maxval_list.insert(std::make_pair(
            sh_sig_in_basis, std::vector<std::pair<size_t, double>>()));

        auto sh_offset_X = 0ul;
        for (auto tile_X = 0ul; tile_X != ntiles_X; ++tile_X) {
          const auto &cluster_X = X_dfbs_->cluster_shells()[tile_X];
          const auto nsh_X = cluster_X.size();

          auto f_X = 0ul;
          for (auto sh_X = 0ul; sh_X != nsh_X; ++sh_X) {
            const auto sh_X_in_basis = sh_X + sh_offset_X;
            const auto nf_sh_X = cluster_X[sh_X].size();

            auto max_val = 0.0;
            for (auto tile_mu = 0ul; tile_mu != ntiles_per_uc_; ++tile_mu) {
              std::array<size_t, 3> idx_Q = {{tile_X, tile_mu, tile_sig}};
              if (!arg_array.is_zero(idx_Q)) {
                // make it parallel
                Tile arg_tile = arg_array.find(idx_Q);
                max_abs_coeff(arg_tile, max_val,
                              std::array<size_t, 2>{{f_sig, f_X}},
                              std::array<size_t, 2>{{nf_sh_sig, nf_sh_X}});
              }
            }
            if (max_val >= threshold)
              shpair_maxval_list[sh_sig_in_basis].emplace_back(
                  std::make_pair(sh_X_in_basis, max_val));

            f_X += nf_sh_X;
          }

          sh_offset_X += nsh_X;
        }

        f_sig += nf_sh_sig;
      }

      sh_offset_sig += nsh_sig;
    }
    world.gop.fence();

    // resort shell list in increasing order
    // note that the customized comparator only compares the first element
    // because comparing the second element is meaningless

    // # of shells in this basis
    const auto nshells_tot_sig = basisRD_->flattened_shells().size();
    for (auto sh_sig = 0; sh_sig != nshells_tot_sig; ++sh_sig) {
      auto &list = shpair_maxval_list[sh_sig];
      std::sort(list.begin(), list.end(),
                [](auto &l, auto &r) { return l.first < r.first; });
    }

    return shpair_maxval_list;
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
