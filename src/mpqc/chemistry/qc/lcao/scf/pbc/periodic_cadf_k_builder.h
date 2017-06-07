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
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;
  using Qmatrix = ::mpqc::lcao::gaussian::Qmatrix;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using ShellVec = typename ::mpqc::lcao::gaussian::ShellVec;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
  using func_offset_list =
      std::unordered_map<size_t, std::tuple<size_t, size_t>>;

  using WorldObject_ =
      madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>>;
  using PeriodicCADFKBuilder_ = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;

  template <int rank>
  using norm_type = std::vector<std::pair<std::array<int, rank>, float>>;

  PeriodicCADFKBuilder(madness::World &world, Factory &ao_factory)
      : WorldObject_(world), ao_factory_(ao_factory) {
    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    print_detail_ = ao_factory_.print_detail();
    force_shape_threshold_ = 1.0e-12;
    target_precision_ = std::numeric_limits<double>::epsilon();
    mpqc::time_point t0, t1;

    // by-cluster orbital basis and df basis
    obs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    dfbs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));
    assert(obs_->nclusters() == dfbs_->nclusters());

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

    auto Y_latt_range = max_latt_range(R_max_, RJ_max_ + RD_max_);
    Y_dfbs_ = shift_basis_origin(*dfbs_, zero_shift_base, Y_latt_range, dcell_);

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

      basisR_ = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);
      basisRD_ = shift_basis_origin(*obs_, zero_shift_base, RD_max_, dcell_);

      eri3_X_trange1_ = eri3_X_dfbs_->create_trange1();
      eri3_bs0_trange1_ = eri3_bs0_->create_trange1();
      eri3_bs1_trange1_ = eri3_bs1_->create_trange1();
      basisRD_trange1_ = basisRD_->create_trange1();
      X_dfbs_trange1_ = X_dfbs_->create_trange1();

      // compute significant shell pair list
      sig_shellpair_list_ = parallel_compute_shellpair_list(*obs_, *basisRJ_);
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
    }
    t1 = mpqc::fenced_now(world);
    auto t_misc = mpqc::duration_in_s(t0, t1);

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K init time decomposition:\n"
                    << "\tC_bra:               " << t_C_bra << " s\n"
                    << "\tM:                   " << t_M << " s\n"
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
  Vector3d dcell_;
  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  int64_t R_size_;
  int64_t RJ_size_;
  int64_t RD_size_;

  array_type C_bra_;
  array_type M_;

  std::unique_ptr<PTC_Builder> three_center_builder_;
  shellpair_list_t sig_shellpair_list_;
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

  TA::TiledRange1 basisRD_trange1_;
  TA::TiledRange1 X_dfbs_trange1_;
  TA::TiledRange1 eri3_X_trange1_;
  TA::TiledRange1 eri3_bs0_trange1_;
  TA::TiledRange1 eri3_bs1_trange1_;
  TA::TiledRange result_trange_;
  std::shared_ptr<TA::Pmap> result_pmap_;
  std::unordered_map<size_t, size_t> eri3_bs0_shell_offset_map_;
  std::unordered_map<size_t, size_t> eri3_bs1_shell_offset_map_;
  madness::ConcurrentHashMap<std::size_t, Tile> local_contr_tiles_;
  madness::ConcurrentHashMap<std::size_t, Tile> global_contr_tiles_;

 private:
  array_type compute_K(const array_type &D, double target_precision) {
    auto &world = ao_factory_.world();

    auto t0_k = mpqc::fenced_now(world);
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    array_type K;

    auto oper_type = libint2::Operator::coulomb;
    const auto screen_thresh = ao_factory_.screen_threshold();
    auto screen_norm_op = ::mpqc::lcao::gaussian::detail::l2Norm;
    auto screen_engine = make_engine_pool(
        oper_type, utility::make_array_of_refs(*dfbs_, *obs_, *obs_),
        libint2::BraKet::xx_xx);

    const auto &trange_C = C_bra_.trange();
    const auto &X_range = trange_C.dim(0).tiles_range();
    const auto &mu_range = trange_C.dim(1).tiles_range();
    const size_t ntiles_per_uc = obs_->nclusters();

    std::shared_ptr<Qmatrix> screen_Qbra, screen_Qket;
    screen_Qbra = std::make_shared<Qmatrix>(
        Qmatrix(world, screen_engine, *Y_dfbs_, screen_norm_op));

    time_point t0, t1;
    double t_Cblock = 0.0;
    double t_Qbra = 0.0;
    double t_Qket = 0.0;
    double t_Eket = 0.0;
    double t_F = 0.0;
    double t_Kpart1 = 0.0;
    double t_Kpart2 = 0.0;

    for (const auto RJ : RJ_list_) {
      t0 = mpqc::fenced_now(world);
      // grab C(X, μ_0, ρ): the Rj block of C(X, μ_0, ρ_Rj)
      std::vector<size_t> C_low{X_range.first, mu_range.first,
                                RJ * ntiles_per_uc};
      std::vector<size_t> C_up{X_range.second, mu_range.second,
                               (RJ + 1) * ntiles_per_uc};
      array_type C;
      C("X, mu, rho") = C_bra_("X, mu, rho").block(C_low, C_up);
      t1 = mpqc::fenced_now(world);
      t_Cblock += mpqc::duration_in_s(t0, t1);

      // if C(X, μ_0, ρ) is zero, skip the rest of this loop
      const double norm_C = C("X, mu, rho").norm();
      if (norm_C == 0.0) continue;

      // compute Q(X, μ_0, σ_Rd) = C(X, μ_0, ρ) * D(ρ, σ_Rd)
      t0 = mpqc::fenced_now(world);
      array_type Q_bra;
      Q_bra("X, mu, sig") = C("X, mu, rho") * D("rho, sig");
      t1 = mpqc::fenced_now(world);
      t_Qbra += mpqc::duration_in_s(t0, t1);

      // compute Q(Y, ν_R, ρ) = C(X, ν_R, σ(Rj+Rd)) * D(ρ, σ_Rd)
      //
      // note: no need to compute C(X, ν_R, σ(Rj+Rd)). Just find the right
      // tiles from C(X, μ_0, ρ_Rj) due to translational symmetry
      t0 = mpqc::fenced_now(world);
      array_type Q_ket;
      Q_ket = compute_Q_ket(C_bra_, D, RJ);
      t1 = mpqc::fenced_now(world);
      t_Qket += mpqc::duration_in_s(t0, t1);

      // compute E(Y, μ_0, ρ)
      t0 = mpqc::fenced_now(world);
      DirectTArray E_ket;
      {
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        auto bs1 = shift_basis_origin(*obs_, vec_RJ);

        auto bs_array = utility::make_array_of_refs(*Y_dfbs_, *obs_, *bs1);
        auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs_, *obs_, *bs1}};

        screen_Qket = std::make_shared<Qmatrix>(
            Qmatrix(world, screen_engine, *obs_, *bs1, screen_norm_op));
        auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
            lcao::gaussian::SchwarzScreen(screen_Qbra, screen_Qket,
                                          screen_thresh));
        auto engine =
            make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);

        E_ket = lcao::gaussian::direct_sparse_integrals(
            world, engine, bs_vector, std::move(screener));
      }
      t1 = mpqc::fenced_now(world);
      t_Eket += mpqc::duration_in_s(t0, t1);

      // compute F(Y, μ_0, ρ) = E(Y, μ_0, ρ) - C(X, μ_0, ρ) M(X, Y)
      t0 = mpqc::fenced_now(world);
      array_type F;
      {
        auto normsp =
            force_normsp(Q_ket.shape().data(), E_ket.array().shape().data());
        auto trange = E_ket.array().trange();
        TA::SparseShape<float> forced_shape(world, normsp, trange);

        F("Y, mu, rho") = (E_ket("Y, mu, rho")).set_shape(forced_shape);
        F.truncate();
        F("Y, mu, rho") -=
            (C("X, mu, rho") * M_("X, Y")).set_shape(forced_shape);
        F.truncate();
      }
      t1 = mpqc::fenced_now(world);
      t_F += mpqc::duration_in_s(t0, t1);

      // compute K(μ_0, ν_R) += F(Y, μ_0, ρ) Q(Y, ν_R, ρ)
      t0 = mpqc::fenced_now(world);
      if (!K.is_initialized()) {
        K("mu, nu") = F("Y, mu, rho") * Q_ket("Y, nu, rho");
      } else {
        K("mu, nu") += F("Y, mu, rho") * Q_ket("Y, nu, rho");
      }
      t1 = mpqc::fenced_now(world);
      t_Kpart1 += mpqc::duration_in_s(t0, t1);

      // compute K(μ_0, ν_R) += E(X, ν_R, σ_Rd) Q(X, μ_0, σ_Rd)
      result_trange_ = K.trange();
      result_pmap_ = K.pmap();
      t0 = mpqc::fenced_now(world);
      auto EQ = compute_contr_EQ(Q_bra, RJ, target_precision);
      K("mu, nu") += EQ("mu, nu");
      t1 = mpqc::fenced_now(world);
      t_Kpart2 += mpqc::duration_in_s(t0, t1);
//      ExEnv::out0() << "\nRJ = " << RJ << ", EQ = \n" << EQ << std::endl;

    }
    auto t1_k = mpqc::fenced_now(world);
    auto t_tot = mpqc::duration_in_s(t0_k, t1_k);
    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K time decomposition:\n"
                    << "\tC block:                  " << t_Cblock << " s\n"
                    << "\tQ_bra:                    " << t_Qbra << " s\n"
                    << "\tQ_ket:                    " << t_Qket << " s\n"
                    << "\tE_ket:                    " << t_Eket << " s\n"
                    << "\tF = E - C M:              " << t_F << " s\n"
                    << "\tK_part1 = F Qket:         " << t_Kpart1 << " s\n"
                    << "\tK_part2 = E Qbra:         " << t_Kpart2 << " s\n"
                    << "\nTotal K builder time:     " << t_tot << " s"
                    << std::endl;
    }

    return K;
  }

  auto force_normsp(TA::Tensor<float> const &in,
                    TA::Tensor<float> const &F_norms) {
    auto const &irange = in.range();
    auto iext = irange.extent_data();

    // Q("Y, nu, rho")
    const auto Y_size = iext[0];
    const auto nu_size = iext[1];
    const auto rho_size = iext[2];

    using SigPair = std::pair<int64_t, int64_t>;
    std::unordered_set<SigPair, boost::hash<SigPair>> Y_rho;
    Y_rho.reserve(Y_size * rho_size);

    for (auto Y = 0ul; Y < Y_size; ++Y) {
      for (auto rho = 0ul; rho < rho_size; ++rho) {
        for (auto nu = 0ul; nu < nu_size; ++nu) {
          const auto val = in(Y, nu, rho);
          if (val > force_shape_threshold_) {
            SigPair Y_rho_pair(Y, rho);
            Y_rho.insert(Y_rho_pair);
            break;
          }
        }
      }
    }

    auto const &orange = F_norms.range();
    TA::Tensor<float> out(orange, 0.0);

    // F("X, nu, sigma")
    auto oext = orange.extent_data();
    assert(Y_size == oext[0]);
    assert(rho_size == oext[2]);
    const auto mu_size = oext[1];

    for (auto mu = 0ul; mu < mu_size; ++mu) {
      for (auto const &Y_rho_pair : Y_rho) {
        const auto Y = Y_rho_pair.first;
        const auto rho = Y_rho_pair.second;
        out(Y, mu, rho) = std::numeric_limits<float>::max();
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

    const auto &C_norms = C_repl.shape().data();

    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };
    auto is_in_lattice_range = [](Vector3i const &in_idx, Vector3i const &range,
                                  Vector3i const &center = {0, 0, 0}) {
      if (in_idx(0) <= center(0) + range(0) &&
          in_idx(0) >= center(0) - range(0) &&
          in_idx(1) <= center(1) + range(1) &&
          in_idx(1) >= center(1) - range(1) &&
          in_idx(2) <= center(2) + range(2) &&
          in_idx(2) >= center(2) - range(2))
        return true;
      else
        return false;
    };

    const auto unshifted_Y_latt_range =
        max_latt_range(R_max_, RJ_max_ + RD_max_);

    const auto shifted_Y_latt_range = RJ_max_;

    const auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
    const auto vec_RJ = direct_vector(RJ_ord, RJ_max_, dcell_);
    const auto bs1 = shift_basis_origin(*obs_, vec_RJ);

    const auto trange = lcao::gaussian::detail::create_trange(
        lcao::gaussian::BasisVector{{*Y_dfbs_, *basisR_, *bs1}});
    const auto tvolume = trange.tiles_range().volume();
    const auto pmap = Policy::default_pmap(world, tvolume);

    // force Q norms
    const size_t ntiles_per_uc = obs_->nclusters();
    const size_t ntiles_Y = Y_dfbs_->nclusters();
    const size_t ntiles_nu = basisR_->nclusters();
    const size_t ntiles_rho = bs1->nclusters();
    const size_t ntiles_sig = ntiles_per_uc * RD_size_;

    TA::Range range;
    range = TA::Range(std::array<size_t, 3>{{ntiles_Y, ntiles_nu, ntiles_rho}});
    TA::Tensor<float> norms(range, 0.0);

    using SigPair = std::pair<size_t, size_t>;
    using TilePairList = std::unordered_map<size_t, std::vector<SigPair>>;
    TilePairList significant_Y_rho;

    for (auto nu = 0ul; nu != ntiles_nu; ++nu) {
      auto R_ord = nu / ntiles_per_uc;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_per_uc;

      significant_Y_rho.insert(std::make_pair(nu, std::vector<SigPair>()));

      for (auto Y = 0ul; Y != ntiles_Y; ++Y) {
        auto RY_ord = Y / ntiles_per_uc;
        auto RY_3D = direct_3D_idx(RY_ord, unshifted_Y_latt_range);

        auto RYmR_3D = RY_3D - R_3D;
        if (!is_in_lattice_range(RYmR_3D, shifted_Y_latt_range)) continue;

        auto RYmR_ord = direct_ord_idx(RYmR_3D, shifted_Y_latt_range);
        auto Y_in_C = Y % ntiles_per_uc + RYmR_ord * ntiles_per_uc;

        for (auto rho = 0ul; rho != ntiles_rho; ++rho) {
          for (auto sig = 0ul; sig != ntiles_sig; ++sig) {
            auto RD_ord = sig / ntiles_per_uc;
            auto RJpRD_3D = RJ_3D + direct_3D_idx(RD_ord, RD_max_);

            if (!(RY_3D == R_3D) && !(RY_3D == RJpRD_3D)) continue;

            auto RJpRDmR_3D = RJpRD_3D - R_3D;
            if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

            auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
            if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) ==
                RJ_list_.end())
              continue;

            auto sig_in_C = sig % ntiles_per_uc + RJpRDmR_ord * ntiles_per_uc;
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
        auto RD_ord = sig / ntiles_per_uc;
        auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
        auto RJpRD_3D = RJ_3D + RD_3D;

        if (!(RY_3D == R_3D) && !(RY_3D == RJpRD_3D)) continue;

        auto RJpRDmR_3D = RJpRD_3D - R_3D;
        if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

        auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
        if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) ==
            RJ_list_.end())
          continue;

        auto sig_in_C = sig % ntiles_per_uc + RJpRDmR_ord * ntiles_per_uc;

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
      auto R_ord = nu / ntiles_per_uc;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_per_uc;

      const auto &rng_nu = trange.dim(1).tile(nu);
      const auto offset_nu = rng_nu.first;
      const auto ext_nu = rng_nu.second - rng_nu.first;

      for (const auto &significant_pair : significant_Y_rho[nu]) {
        auto Y = significant_pair.first;
        auto rho = significant_pair.second;

        std::array<size_t, 3> tile_idx = {{Y, nu, rho}};
        const size_t tile_ord = trange.tiles_range().ordinal(tile_idx);

        if (Q.is_local(tile_idx) && !Q.is_zero(tile_idx)) {
          auto RY_ord = Y / ntiles_per_uc;
          auto RY_3D = direct_3D_idx(RY_ord, unshifted_Y_latt_range);
          auto Y_in_C = Y % ntiles_per_uc +
                        direct_ord_idx(RY_3D - R_3D, shifted_Y_latt_range) *
                            ntiles_per_uc;

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

  array_type compute_contr_EQ(array_type const &Q, int64_t const RJ,
                              double target_precision) {
    auto &world = this->get_world();
    const auto me = world.rank();
    const auto nproc = world.nproc();
    target_precision_ = target_precision;

    // # of tiles per basis
    const auto ntiles_per_uc = obs_->nclusters();
    const auto ntiles_mu = obs_->nclusters();
    const auto ntiles_nu = basisR_->nclusters();
    const auto ntiles_X = X_dfbs_->nclusters();
    const auto ntiles_sig = basisRD_->nclusters();

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;
    const auto RJ_3D = direct_3D_idx(RJ, RJ_max_);

    auto is_in_lattice_range = [](Vector3i const &in_idx, Vector3i const &range,
                                  Vector3i const &center = {0, 0, 0}) {
      if (in_idx(0) <= center(0) + range(0) &&
          in_idx(0) >= center(0) - range(0) &&
          in_idx(1) <= center(1) + range(1) &&
          in_idx(1) >= center(1) - range(1) &&
          in_idx(2) <= center(2) + range(2) &&
          in_idx(2) >= center(2) - range(2))
        return true;
      else
        return false;
    };

    for (auto tile_X = 0ul, task = 0ul; tile_X != ntiles_X; ++tile_X) {
      const auto RX_ord = tile_X / ntiles_per_uc;
      const auto RX_3D = direct_3D_idx(RX_ord, RJ_max_);
      for (auto tile_mu = 0ul; tile_mu != ntiles_mu; ++tile_mu) {
        for (auto tile_sig = 0ul; tile_sig != ntiles_sig; ++tile_sig) {
          std::array<size_t, 3> idx_Q = {{tile_X, tile_mu, tile_sig}};
          if (Q.is_zero(idx_Q)) continue;

          auto Qtile = Q.find(idx_Q);

          const auto RD_ord = tile_sig / ntiles_per_uc;
          const auto RD_3D = direct_3D_idx(RD_ord, RD_max_);

          for (auto tile_nu = 0ul; tile_nu != ntiles_nu; ++tile_nu, ++task) {
            const auto R_ord = tile_nu / ntiles_per_uc;
            const auto R_3D = direct_3D_idx(R_ord, R_max_);

            const auto RJpRDmR_3D = RJ_3D + RD_3D - R_3D;
            if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

            auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
            if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) ==
                RJ_list_.end())
              continue;

            const auto RXmR_3D = RX_3D - R_3D;
            const auto RXmR_ord = direct_ord_idx(RXmR_3D, RJ_max_ + R_max_);

            const auto eri3_tile_X =
                tile_X % ntiles_per_uc + RXmR_ord * ntiles_per_uc;
            const auto eri3_tile_nu = tile_nu % ntiles_per_uc;
            const auto eri3_tile_sig =
                tile_sig % ntiles_per_uc + RJpRDmR_ord * ntiles_per_uc;

            if (task % nproc == me) {
              WorldObject_::task(
                  me, &PeriodicCADFKBuilder_::compute_contr_EQ_task, Qtile, RJ,
                  std::array<size_t, 4>{{tile_mu, tile_nu, tile_X, tile_sig}},
                  std::array<size_t, 3>{{eri3_tile_X,
                                         eri3_tile_nu, eri3_tile_sig}});
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

  void compute_contr_EQ_task(Tile Q, int64_t RJ, std::array<size_t, 4> tile_idx,
                             std::array<size_t, 3> eri3_idx) {
    const auto tile_mu = tile_idx[0];
    const auto tile_nu = tile_idx[1];
    const auto tile_X = tile_idx[2];
    const auto tile_sig = tile_idx[3];

    // get tile indices of eri3
    const auto eri3_tile_X = eri3_idx[0];
    const auto eri3_tile_nu = eri3_idx[1];
    const auto eri3_tile_sig = eri3_idx[2];

    // get reference to basis sets
    const auto &eri3_basis_X = eri3_X_dfbs_;
    const auto &eri3_basis_nu = eri3_bs0_;
    const auto &eri3_basis_sig = eri3_bs1_;

    // shell clusters for this tile
    const auto &eri3_cluster_X = eri3_basis_X->cluster_shells()[eri3_tile_X];
    const auto &eri3_cluster_nu = eri3_basis_nu->cluster_shells()[eri3_tile_nu];
    const auto &eri3_cluster_sig =
        eri3_basis_sig->cluster_shells()[eri3_tile_sig];

    // # of shells in each cluster
    const auto eri3_nshells_X = eri3_cluster_X.size();
    const auto eri3_nshells_nu = eri3_cluster_nu.size();
    const auto eri3_nshells_sig = eri3_cluster_sig.size();

    // 1-d tile ranges
    const auto &tr0 = result_trange_.dim(0);
    const auto &tr1 = result_trange_.dim(1);
    const auto ntiles_nu = tr1.tile_extent();
    const auto &rng_mu = tr0.tile(tile_mu);
    const auto &rng_nu = tr1.tile(tile_nu);
    const auto &rng_X = X_dfbs_trange1_.tile(tile_X);
    const auto &rng_sig = basisRD_trange1_.tile(tile_sig);

    // TODO make eri3-specific trange1 outside of task
    const auto &eri3_rng_X = eri3_X_trange1_.tile(eri3_tile_X);
    const auto &eri3_rng_nu = eri3_bs0_trange1_.tile(eri3_tile_nu);
    const auto &eri3_rng_sig = eri3_bs1_trange1_.tile(eri3_tile_sig);

    // range sizes
    const auto rng_size_mu = rng_mu.second - rng_mu.first;
    const auto rng_size_nu = rng_nu.second - rng_nu.first;
    const auto rng_size_X = rng_X.second - rng_X.first;
    const auto rng_size_sig = rng_sig.second - rng_sig.first;

    // 2-d tile ranges describing the contribution blocks produced by this
    auto result_rng = TA::Range({rng_mu, rng_nu});
    // initialize contribution to the result matrices
    auto result_tile = Tile(std::move(result_rng), 0.0);

    // grab ptrs to tile data to make addressing more efficient
    auto *result_ptr = result_tile.data();
    const auto *Q_ptr = Q.data();
    assert(Q_ptr != nullptr);

    // compute eri3 * Q contribution to all Fock matrices
    {
      // index of first shell in this cluster
      const auto eri3_sh_offset_nu = eri3_bs0_shell_offset_map_[eri3_tile_nu];
      const auto eri3_sh_offset_sig = eri3_bs1_shell_offset_map_[eri3_tile_sig];

      // index of last shell in this cluster
      const auto eri3_sh_max_nu = eri3_sh_offset_nu + eri3_nshells_nu;
      const auto eri3_sh_max_sig = eri3_sh_offset_sig + eri3_nshells_sig;

      auto is_significant = false;
      {
        for (auto eri3_sh_nu = eri3_sh_offset_nu; eri3_sh_nu != eri3_sh_max_nu;
             ++eri3_sh_nu) {
          for (const auto eri3_sh_sig : sig_shellpair_list_[eri3_sh_nu]) {
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
        auto engine = engines_->local();
        const auto engine_precision = target_precision_;
        engine.set_precision(engine_precision);
        const auto &computed_shell_sets = engine.results();

        // compute offset list of cluster_sig
        auto eri3_offset_list_sig =
            compute_func_offset_list(eri3_cluster_sig, eri3_rng_sig.first);

        // this is the index of the first basis functions for each shell *in
        // this shell cluster*
        auto eri3_cf_offset_nu = 0;
        auto cf_offset_nu = 0;
        // this is the index of the first basis functions for each shell *in the
        // basis set*
        auto eri3_bf_offset_nu = eri3_rng_nu.first;
        auto bf_offset_nu = rng_nu.first;

        size_t eri3_cf_offset_sig, eri3_bf_offset_sig;

        // loop over all shell sets
        for (auto eri3_sh_nu = 0; eri3_sh_nu != eri3_nshells_nu; ++eri3_sh_nu) {
          const auto &eri3_shell_nu = eri3_cluster_nu[eri3_sh_nu];
          const auto eri3_nf_nu = eri3_shell_nu.size();

          const auto sh_nu_in_basis = eri3_sh_nu + eri3_sh_offset_nu;
          auto cf_offset_sig = 0;
          auto bf_offset_sig = rng_sig.first;
          for (auto const &sh_sig_in_basis :
               sig_shellpair_list_[sh_nu_in_basis]) {
            if (sh_sig_in_basis < eri3_sh_offset_sig ||
                sh_sig_in_basis >= eri3_sh_max_sig)
              continue;

            const auto eri3_sh_sig = sh_sig_in_basis - eri3_sh_offset_sig;
            std::tie(eri3_cf_offset_sig, eri3_bf_offset_sig) =
                eri3_offset_list_sig[eri3_sh_sig];

            const auto &eri3_shell_sig = eri3_cluster_sig[eri3_sh_sig];
            const auto eri3_nf_sig = eri3_shell_sig.size();

            auto eri3_cf_offset_X = 0;
            auto eri3_bf_offset_X = eri3_rng_X.first;
            auto cf_offset_X = 0;
            auto bf_offset_X = rng_X.first;
            for (auto eri3_sh_X = 0; eri3_sh_X != eri3_nshells_X; ++eri3_sh_X) {
              const auto &eri3_shell_X = eri3_cluster_X[eri3_sh_X];
              const auto eri3_nf_X = eri3_shell_X.size();

              if (screen.skip(eri3_bf_offset_X, eri3_bf_offset_nu,
                              eri3_bf_offset_sig))
                continue;

              // compute shell set
              engine.compute(eri3_shell_X, eri3_shell_nu, eri3_shell_sig);
              const auto &eri3 = computed_shell_sets[0];

              if (eri3 != nullptr) {
                for (auto eri3_f_X = 0, eri3_ord = 0; eri3_f_X != eri3_nf_X;
                     ++eri3_f_X) {
                  const auto cf_X = eri3_f_X + cf_offset_X;
                  for (auto eri3_f_nu = 0; eri3_f_nu != eri3_nf_nu;
                       ++eri3_f_nu) {
                    const auto cf_nu = eri3_f_nu + cf_offset_nu;
                    for (auto eri3_f_sig = 0; eri3_f_sig != eri3_nf_sig;
                         ++eri3_f_sig, ++eri3_ord) {
                      const auto cf_sig = eri3_f_sig + cf_offset_sig;

                      const auto eri3_value = eri3[eri3_ord];

                      for (auto cf_mu = 0; cf_mu != rng_size_mu; ++cf_mu) {
                        const auto cf_mu_nu = cf_mu * rng_size_nu + cf_nu;
                        const auto cf_X_mu_sig =
                            cf_X * rng_size_mu * rng_size_sig +
                            cf_mu * rng_size_sig + cf_sig;
                        result_ptr[cf_mu_nu] += eri3_value * Q_ptr[cf_X_mu_sig];
                      }
                    }
                  }
                }
              }

              eri3_cf_offset_X += eri3_nf_X;
              eri3_bf_offset_X += eri3_nf_X;
              cf_offset_X += eri3_nf_X;
              bf_offset_X += eri3_nf_X;
            }

            cf_offset_sig += eri3_nf_sig;
            bf_offset_sig += eri3_nf_sig;
          }

          eri3_cf_offset_nu += eri3_nf_nu;
          eri3_bf_offset_nu += eri3_nf_nu;
          cf_offset_nu += eri3_nf_nu;
          bf_offset_nu += eri3_nf_nu;
        }
      }
    }

    // accumulate the local contributions
    {
      const auto tile_ord = tile_mu * ntiles_nu + tile_nu;
      PeriodicCADFKBuilder_::accumulate_local_task(result_tile, tile_ord);
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

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
