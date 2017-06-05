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
class PeriodicCADFKBuilder : public madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using DirectTArray = typename Factory::DirectTArray;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;
  using Qmatrix = ::mpqc::lcao::gaussian::Qmatrix;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

  using WorldObject_ = madness::WorldObject<PeriodicCADFKBuilder<Tile, Policy, Factory>>;
  using PeriodicCADFKBuilder_ = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  using Engine = ::mpqc::lcao::gaussian::ShrPool<libint2::Engine>;

  template <int rank>
  using norm_type = std::vector<std::pair<std::array<int, rank>, float>>;

  PeriodicCADFKBuilder(madness::World &world, Factory &ao_factory) : WorldObject_(world), ao_factory_(ao_factory) {

    // WorldObject mandates this is called from the ctor
    WorldObject_::process_pending();

    print_detail_ = ao_factory_.print_detail();
    force_shape_threshold_ = 1.0e-12;
    target_precision_ = std::numeric_limits<double>::epsilon();
    mpqc::time_point t0, t1;

    // by-cluster orbital basis and df basis
    obs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    dfbs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));

    dcell_ = ao_factory_.unitcell().dcell();
    R_max_ = ao_factory_.R_max();
    RJ_max_ = R_max_;  // replace RJ range with R range
    RD_max_ = ao_factory_.RD_max();
    R_size_ = ao_factory_.R_size();
    RJ_size_ = R_size_;  // replace RJ range with R range
    RD_size_ = ao_factory_.RD_size();

    const Vector3i ref_latt_range = {0, 0, 0};
    const Vector3i ref_latt_center = {0, 0, 0};
    const auto natoms_per_uc = ao_factory_.unitcell().natoms();
    Vector3d zero_shift_base(0.0, 0.0, 0.0);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    // compute C(X, μ_0, ρ_Rj)
    t0 = mpqc::fenced_now(world);
    X_dfbs_ = shift_basis_origin(*dfbs_, zero_shift_base, RJ_max_, dcell_);
    {
      C_bra_ = std::vector<array_type>(RJ_size_, array_type());

      ExEnv::out0() << "\nComputing C_bra ...\n";
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*X_dfbs_);

      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        array_type &C = C_bra_[RJ];
        auto RJ_3D = direct_3D_idx(RJ, RJ_max_);
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        auto bs1 = shift_basis_origin(*obs_, vec_RJ);
        C = lcao::cadf_fitting_coefficients<Tile, Policy>(
            M, *obs_, *bs1, *X_dfbs_, natoms_per_uc, ref_latt_range,
            ref_latt_range, RJ_max_, ref_latt_center, RJ_3D, ref_latt_center);
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_C_bra = mpqc::duration_in_s(t0, t1);

    // test new C_bra
    {
      auto t0 = mpqc::fenced_now(world);
      ExEnv::out0() << "\nComputing new C_bra ...\n";
      auto bs1 = shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);

      const auto by_atom_dfbs = lcao::detail::by_center_basis(*X_dfbs_);
      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

      C_bra_new_ = lcao::cadf_fitting_coefficients<Tile, Policy>(
          M, *obs_, *bs1, *X_dfbs_, natoms_per_uc, ref_latt_range, RJ_max_,
          RJ_max_);

      auto t1 = mpqc::fenced_now(world);
      auto dur = mpqc::duration_in_s(t0, t1);
      detail::print_size_info(C_bra_new_, "C bra");
      ExEnv::out0() << "\tnew C_bra:            " << dur << " s\n" << std::endl;

      // compute significant shell pair list
      sig_shellpair_list_ = parallel_compute_shellpair_list(*obs_, *bs1);
      // make a list of significant Rj's as in overlap between μ and ρ_Rj
      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        const auto nshells = obs_->flattened_shells().size();
        const auto shell1_min = nshells * RJ;
        const auto shell1_max = shell1_min + nshells;

        auto is_significant = false;
        for (auto shell0 = 0; shell0 != nshells; ++shell0) {
          for (const auto &shell1: sig_shellpair_list_[shell0]) {
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

    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };

    // compute C(Y, ν_R, σ_(Rj+Rd))
    t0 = mpqc::fenced_now(world);
    auto Y_latt_range = max_latt_range(R_max_, RJ_max_ + RD_max_);
    Y_dfbs_ =
        shift_basis_origin(*dfbs_, zero_shift_base, Y_latt_range, dcell_);
    {
      C_ket_ = std::vector<array_type>(RJ_size_, array_type());

      ExEnv::out0() << "\nComputing C_ket ...\n";
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*Y_dfbs_);

      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);
      auto bs0 = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);
      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        array_type &C = C_ket_[RJ];
        auto RJ_3D = direct_3D_idx(RJ, RJ_max_);
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        auto bs1 = shift_basis_origin(*obs_, vec_RJ, RD_max_, dcell_);
        C = lcao::cadf_fitting_coefficients<Tile, Policy>(
            M, *bs0, *bs1, *Y_dfbs_, natoms_per_uc, R_max_, RD_max_,
            Y_latt_range, ref_latt_center, RJ_3D, ref_latt_center);
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_C_ket = mpqc::duration_in_s(t0, t1);

    // test new C_ket
//    {
//      auto t0 = mpqc::fenced_now(world);
//      ExEnv::out0() << "\nComputing new C_ket ...\n";
//      auto max_range = R_max_ + RJ_max_ + RD_max_;
//      auto bs1 = shift_basis_origin(*obs_, zero_shift_base, max_range, dcell_);

//      auto Y_dfbs =
//          shift_basis_origin(*dfbs_, zero_shift_base, max_range, dcell_);
//      const auto by_atom_dfbs = lcao::detail::by_center_basis(*Y_dfbs);
//      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

//      C_ket_new_ = lcao::cadf_fitting_coefficients<Tile, Policy>(
//          M, *obs_, *bs1, *Y_dfbs, natoms_per_uc, ref_latt_range, max_range,
//          max_range);

//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);
//      detail::print_size_info(C_ket_new_, "C ket");
//      ExEnv::out0() << "\tnew C_ket:            " << dur << " s\n" << std::endl;
//    }

    // compute M(X, Y)
    t0 = mpqc::fenced_now(world);
    M_ = compute_eri2(world, *X_dfbs_, *Y_dfbs_);
    t1 = mpqc::fenced_now(world);
    double t_M = mpqc::duration_in_s(t0, t1);

    auto oper_type = libint2::Operator::coulomb;
    const auto screen_thresh = ao_factory_.screen_threshold();
    auto screen_norm_op = ::mpqc::lcao::gaussian::detail::l2Norm;
    auto screen_engine = make_engine_pool(
        oper_type, utility::make_array_of_refs(*dfbs_, *obs_, *obs_),
        libint2::BraKet::xx_xx);

    // compute E(X, ν_R, σ_(Rj+Rd))
    t0 = mpqc::fenced_now(world);
    {
      if (E_bra_.empty())
        E_bra_ = std::vector<DirectTArray>(RJ_size_, DirectTArray());

      auto bs0 = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);

      std::shared_ptr<Qmatrix> Qbra, Qket;
      Qbra = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *X_dfbs_, screen_norm_op));

      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        DirectTArray &E = E_bra_[RJ];
        if (!E.array().is_initialized()) {
          auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
          auto bs1 = shift_basis_origin(*obs_, vec_RJ, RD_max_, dcell_);

          auto bs_array = utility::make_array_of_refs(*X_dfbs_, *bs0, *bs1);
          auto bs_vector = lcao::gaussian::BasisVector{{*X_dfbs_, *bs0, *bs1}};

          Qket = std::make_shared<Qmatrix>(
              Qmatrix(world, screen_engine, *bs0, *bs1, screen_norm_op));
          auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
              lcao::gaussian::SchwarzScreen(Qbra, Qket, screen_thresh));

          auto engine =
              make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);

          E = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                      std::move(screener));
        }
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_E_bra = mpqc::duration_in_s(t0, t1);

    // test new E_bra
//    {
//      auto t0 = mpqc::fenced_now(world);
//      ExEnv::out0() << "\nComputing new E_bra ...\n";
//      std::shared_ptr<Qmatrix> Qbra, Qket;
//      Qbra = std::make_shared<Qmatrix>(
//          Qmatrix(world, screen_engine, *X_dfbs_, screen_norm_op));

//      auto max_range = R_max_ + RJ_max_ + RD_max_;
//      auto bs1 = shift_basis_origin(*obs_, zero_shift_base, max_range, dcell_);

//      Qket = std::make_shared<Qmatrix>(
//          Qmatrix(world, screen_engine, *obs_, *bs1, screen_norm_op));
//      auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
//          lcao::gaussian::SchwarzScreen(Qbra, Qket, screen_thresh));

//      auto bs_array = utility::make_array_of_refs(*X_dfbs_, *obs_, *bs1);
//      auto bs_vector = lcao::gaussian::BasisVector{{*X_dfbs_, *obs_, *bs1}};
//      auto engine =
//          make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);
//      E_bra_new_ = lcao::gaussian::direct_sparse_integrals(
//          world, engine, bs_vector, std::move(screener));
//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);
//      ExEnv::out0() << "\tnew E_bra:            " << dur << " s\n" << std::endl;
//    }

    // compute E(Y, μ_0, ρ_Rj)
    t0 = mpqc::fenced_now(world);
    {
      if (E_ket_.empty())
        E_ket_ = std::vector<DirectTArray>(RJ_size_, DirectTArray());

      std::shared_ptr<Qmatrix> Qbra, Qket;
      Qbra = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *Y_dfbs_, screen_norm_op));

      for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
        DirectTArray &E = E_ket_[RJ];
        if (!E.array().is_initialized()) {
          auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
          auto bs1 = shift_basis_origin(*obs_, vec_RJ);

          auto bs_array = utility::make_array_of_refs(*Y_dfbs_, *obs_, *bs1);
          auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs_, *obs_, *bs1}};

          Qket = std::make_shared<Qmatrix>(
              Qmatrix(world, screen_engine, *obs_, *bs1, screen_norm_op));
          auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
              lcao::gaussian::SchwarzScreen(Qbra, Qket, screen_thresh));

          auto engine =
              make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);

          E = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                      std::move(screener));
        }
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_E_ket = mpqc::duration_in_s(t0, t1);

    // test new E_ket
//    {
//      auto t0 = mpqc::fenced_now(world);
//      ExEnv::out0() << "\nComputing new E_ket ...\n";
//      std::shared_ptr<Qmatrix> Qbra, Qket;
//      Qbra = std::make_shared<Qmatrix>(
//          Qmatrix(world, screen_engine, *Y_dfbs_, screen_norm_op));

//      auto bs1 = shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);
//      Qket = std::make_shared<Qmatrix>(
//          Qmatrix(world, screen_engine, *obs_, *bs1, screen_norm_op));
//      auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
//          lcao::gaussian::SchwarzScreen(Qbra, Qket, screen_thresh));

//      auto bs_array = utility::make_array_of_refs(*Y_dfbs_, *obs_, *bs1);
//      auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs_, *obs_, *bs1}};
//      auto engine =
//          make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);
//      E_ket_new_ = lcao::gaussian::direct_sparse_integrals(
//          world, engine, bs_vector, std::move(screener));

//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);
//      ExEnv::out0() << "\tnew E_ket:            " << dur << " s\n" << std::endl;
//    }

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K init time decomposition:\n"
                    << "\tC_bra:               " << t_C_bra << " s\n"
                    << "\tC_ket:               " << t_C_ket << " s\n"
                    << "\tM:                   " << t_M << " s\n"
                    << "\tE_bra:               " << t_E_bra << " s\n"
                    << "\tE_ket:               " << t_E_ket << " s"
                    << std::endl;
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
  double target_precision_ = 0.0;
  std::unique_ptr<PTC_Builder> three_center_builder_;
  shellpair_list_t sig_shellpair_list_;
  std::vector<int64_t> RJ_list_;

  std::shared_ptr<Basis> obs_;
  std::shared_ptr<Basis> dfbs_;
  std::shared_ptr<Basis> X_dfbs_;
  std::shared_ptr<Basis> Y_dfbs_;

  Vector3d dcell_;
  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  int64_t R_size_;
  int64_t RJ_size_;
  int64_t RD_size_;

  std::vector<array_type> C_bra_;
  std::vector<array_type> C_ket_;
  std::vector<array_type> F_;
  array_type M_;
  std::vector<DirectTArray> E_bra_;
  std::vector<DirectTArray> E_ket_;

  array_type C_bra_new_;
  array_type C_ket_new_;
  DirectTArray E_bra_new_;
  DirectTArray E_ket_new_;

  std::shared_ptr<lcao::Screener> p_screener_;


 private:
  array_type compute_K(const array_type &D, double target_precision) {
    auto &world = ao_factory_.world();
    //    auto RJ_size = ao_factory_.RJ_size();
    auto R_size = ao_factory_.R_size();
    auto RJ_size = R_size;

    mpqc::time_point t0, t1;
    auto t0_k_builder = mpqc::fenced_now(world);

    // compute Q(X, μ_0, σ_(Rj+Rd)) = C(X, μ_0, ρ_Rj) D(ρ_0, σ_Rd)
    t0 = mpqc::fenced_now(world);
    std::vector<array_type> Q_bra(RJ_size, array_type());
    for (auto RJ = 0; RJ != RJ_size; ++RJ) {
      array_type &C = C_bra_[RJ];
      const double norm = C("X, mu, rho").norm();

      if (norm == 0.0) {
        continue;
      }
      array_type &Q = Q_bra[RJ];
      Q("X, mu, sig") = C("X, mu, rho") * D("rho, sig");
    }
    t1 = mpqc::fenced_now(world);
    double t_Q_bra = mpqc::duration_in_s(t0, t1);

    // test new Q_bra
//    {
//      auto t0 = mpqc::fenced_now(world);
//      ExEnv::out0() << "\nComputing new Q_bra ...\n";

//      array_type Q_bra_new;
//      Q_bra_new = compute_Q_bra(C_bra_new_, D);
//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);

//      detail::print_size_info(Q_bra_new, "new Q_bra");
//      ExEnv::out0() << "\tnew Q_bra:            " << dur << " s\n";
//    }

    t0 = mpqc::fenced_now(world);
    array_type F, K_part1;
    double t_force_shape = 0.0;
    double t_build_f = 0.0;
    double t_contr = 0.0;
    for (auto RJ = 0; RJ != RJ_size; ++RJ) {
      array_type &C = C_ket_[RJ];
      array_type &Q = Q_bra[RJ];
      if (!Q.is_initialized()) {
        continue;
      }
      DirectTArray &E = E_bra_[RJ];

      mpqc::time_point t0, t1;
      // force shape of F
      // auto norms = force_norms(Q_bra[RJ].shape().data(), 1, R_size);
      t0 = mpqc::fenced_now(world);
      auto normsp = force_normsp(Q.shape().data(), E.array().shape().data());
      auto trange = E.array().trange();
      TA::SparseShape<float> forced_shape(world, normsp, trange);
      t1 = mpqc::fenced_now(world);
      t_force_shape += mpqc::duration_in_s(t0, t1);

      // compute F(X, ν_R, σ_(Rj+Rd)) =
      //           E(X, ν_R, σ_(Rj+Rd)) - M(X, Y) C(Y, ν_R, σ_(Rj+Rd))
      t0 = mpqc::fenced_now(world);
      F("X, nu, sig") = (E("X, nu, sig")).set_shape(forced_shape);
      F.truncate();
      F("X, nu, sig") -= (M_("X, Y") * C("Y, nu, sig")).set_shape(forced_shape);
      F.truncate();
      world.gop.fence();

      t1 = mpqc::fenced_now(world);
      t_build_f += mpqc::duration_in_s(t0, t1);

      t0 = mpqc::fenced_now(world);
      if (!K_part1.is_initialized()) {
        K_part1("mu, nu") = Q("X, mu, sig") * F("X, nu, sig");
      } else {
        K_part1("mu, nu") += Q("X, mu, sig") * F("X, nu, sig");
      }
      t1 = mpqc::fenced_now(world);
      t_contr += mpqc::duration_in_s(t0, t1);
    }
    t1 = mpqc::fenced_now(world);
    double t_F_and_K1 = mpqc::duration_in_s(t0, t1);

    // compute Q(Y, ν_R, ρ_Rj) = C(Y, ν_R, σ_(Rj+Rd)) D(ρ_0, σ_Rd)
    t0 = mpqc::fenced_now(world);
    std::vector<array_type> Q_ket(RJ_size, array_type());
    for (auto RJ = 0; RJ != RJ_size; ++RJ) {
      array_type &C = C_ket_[RJ];
      const double norm = C("Y, nu, sig").norm();
      if (norm == 0.0) {
        continue;
      }
      array_type &Q = Q_ket[RJ];
      Q("Y, nu, rho") = C("Y, nu, sig") * D("rho, sig");
    }
    t1 = mpqc::fenced_now(world);
    double t_Q_ket = mpqc::duration_in_s(t0, t1);

    // test new Q_ket
//    array_type Q_ket_new;
//    {
//      auto t0 = mpqc::fenced_now(world);
//      ExEnv::out0() << "\nComputing new Q_ket ...\n";

//      Q_ket_new = compute_Q_ket(C_ket_new_, D);
//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);

//      detail::print_size_info(Q_ket_new, "new Q_ket");
//      ExEnv::out0() << "\tnew Q_ket:            " << dur << " s\n";
//    }

    // test new F
//    array_type F_new;
//    {
//      auto t0 = mpqc::fenced_now(world);
//      ExEnv::out0() << "\nComputing new F ...\n";

//      auto normsp = force_normsp(Q_ket_new.shape().data(),
//                                 E_ket_new_.array().shape().data());
//      auto trange = E_ket_new_.array().trange();
//      TA::SparseShape<float> forced_shape(world, normsp, trange);

//      F_new("Y, mu, rho") = (E_ket_new_("Y, mu, rho")).set_shape(forced_shape);
//      F_new.truncate();
//      F_new("Y, mu, rho") -=
//          (M_("X, Y") * C_bra_new_("X, mu, rho")).set_shape(forced_shape);
//      F_new.truncate();
//      world.gop.fence();

//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);

//      detail::print_size_info(F_new, "new F");
//      ExEnv::out0() << "\tnew F:            " << dur << " s\n";
//    }

    // test new K_part1
//    {
//      auto t0 = mpqc::fenced_now(world);
//      array_type K_part1_new;
//      ExEnv::out0() << "\nComputing new K part1 = F Q_ket ...\n";
//      K_part1_new("mu, nu") = F_new("Y, mu, rho") * Q_ket_new("Y, nu, rho");
//      auto t1 = mpqc::fenced_now(world);
//      auto dur = mpqc::duration_in_s(t0, t1);
//      detail::print_size_info(F_new, "new K part1");
//      ExEnv::out0() << "\tnew K part1:            " << dur << " s\n";
//    }

    // compute K_part2
    t0 = mpqc::fenced_now(world);
    array_type K_part2;
    for (auto RJ = 0; RJ != RJ_size; ++RJ) {
      auto const &Q = Q_ket[RJ];
      auto const &E = E_ket_[RJ];
      if (Q_ket[RJ].is_initialized()) {
        auto normsp = force_normsp(Q.shape().data(), E.array().shape().data());

        auto trange = E.array().trange();
        TA::SparseShape<float> forced_shape(world, normsp, trange);

        array_type F2;
        F2("Y, mu, rho") = (E("Y, mu, rho")).set_shape(forced_shape);
        if (!K_part2.is_initialized()) {
          K_part2("mu, nu") = F2("Y, mu, rho") * Q("Y, nu, rho");
        } else {
          K_part2("mu, nu") += F2("Y, mu, rho") * Q("Y, nu, rho");
        }
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_K_part2 = mpqc::duration_in_s(t0, t1);

    // compute exchange
    t0 = mpqc::fenced_now(world);
    array_type K;
    K("mu, nu") = K_part1("mu, nu") + K_part2("mu, nu");
    t1 = mpqc::fenced_now(world);
    double t_K = mpqc::duration_in_s(t0, t1);

    auto t1_k_builder = mpqc::fenced_now(world);
    double t_tot = mpqc::duration_in_s(t0_k_builder, t1_k_builder);

    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K time decomposition:\n"
                    << "\tQ_bra:                " << t_Q_bra << " s\n"
                    << "\tF and Q F:            " << t_F_and_K1 << " s\n"
                    << "\t\tforce shape:        " << t_force_shape << " s\n"
                    << "\t\tbuild f:            " << t_build_f << " s\n"
                    << "\t\tcontraction:        " << t_contr << " s\n"
                    << "\tQ_ket:                " << t_Q_ket << " s\n"
                    << "\tK_part2 = E Q:        " << t_K_part2 << " s\n"
                    << "\tK = K1 + K2:          " << t_K << " s\n"
                    << "\nTotal K builder time: " << t_tot << " s" << std::endl;
    }

    // ***************** try something different *****************
    auto t0_new_k = mpqc::fenced_now(world);
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::make_engine_pool;

    array_type K_new;

    auto oper_type = libint2::Operator::coulomb;
    const auto screen_thresh = ao_factory_.screen_threshold();
    auto screen_norm_op = ::mpqc::lcao::gaussian::detail::l2Norm;
    auto screen_engine = make_engine_pool(
        oper_type, utility::make_array_of_refs(*dfbs_, *obs_, *obs_),
        libint2::BraKet::xx_xx);

    const auto &trange_C = C_bra_new_.trange();
    const auto &X_range = trange_C.dim(0).tiles_range();
    const auto &mu_range = trange_C.dim(1).tiles_range();
    const size_t ntiles_per_uc = obs_->nclusters();

    std::shared_ptr<Qmatrix> screen_Qbra, screen_Qket;
    screen_Qbra = std::make_shared<Qmatrix>(
        Qmatrix(world, screen_engine, *Y_dfbs_, screen_norm_op));

    time_point t00, t11;
    double t_Cblock = 0.0;
    double t_Qbra_new = 0.0;
    double t_Qket_new = 0.0;
    double t_Eket_new = 0.0;
    double t_F_new = 0.0;
    double t_Kpart1_new = 0.0;

    for (const auto RJ : RJ_list_) {
      t00 = mpqc::fenced_now(world);
      std::vector<size_t> C_low{X_range.first, mu_range.first, RJ * ntiles_per_uc};
      std::vector<size_t> C_up{X_range.second, mu_range.second, (RJ + 1) * ntiles_per_uc};
      array_type C;
      C("X, mu, rho") = C_bra_new_("X, mu, rho").block(C_low, C_up);
      t11 = mpqc::fenced_now(world);
      t_Cblock += mpqc::duration_in_s(t00, t11);

      const double norm_C = C("X, mu, rho").norm();
      if (norm_C == 0.0) continue;

      t00 = mpqc::fenced_now(world);
      array_type Q_bra;
      Q_bra("X, mu, sig") = C("X, mu, rho") * D("rho, sig");
      t11 = mpqc::fenced_now(world);
      t_Qbra_new += mpqc::duration_in_s(t00, t11);

      detail::print_size_info(Q_bra, "Q_bra_new");

      t00 = mpqc::fenced_now(world);
      array_type Q_ket;
      Q_ket = compute_Q_ket(C_bra_new_, D, RJ);
      t11 = mpqc::fenced_now(world);
      t_Qket_new += mpqc::duration_in_s(t00, t11);
      detail::print_size_info(Q_ket, "Q_ket_new");

      detail::print_size_info(D, "D");

      t00 = mpqc::fenced_now(world);
      DirectTArray E_ket;
      {
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        auto bs1 = shift_basis_origin(*obs_, vec_RJ);

        auto bs_array = utility::make_array_of_refs(*Y_dfbs_, *obs_, *bs1);
        auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs_, *obs_, *bs1}};

        screen_Qket = std::make_shared<Qmatrix>(
            Qmatrix(world, screen_engine, *obs_, *bs1, screen_norm_op));
        auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
            lcao::gaussian::SchwarzScreen(screen_Qbra, screen_Qket, screen_thresh));
        auto engine =
            make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);

        E_ket = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                    std::move(screener));
      }
      t11 = mpqc::fenced_now(world);
      t_Eket_new += mpqc::duration_in_s(t00, t11);

      t00 = mpqc::fenced_now(world);
      array_type F;
      {
        auto normsp = force_normsp(Q_ket.shape().data(), E_ket.array().shape().data());
        auto trange = E_ket.array().trange();
        TA::SparseShape<float> forced_shape(world, normsp, trange);

        F("Y, mu, rho") = (E_ket("Y, mu, rho")).set_shape(forced_shape);
        F.truncate();
        F("Y, mu, rho") -= (C("X, mu, rho") * M_("X, Y")).set_shape(forced_shape);
        F.truncate();
      }
      t11 = mpqc::fenced_now(world);
      t_F_new += mpqc::duration_in_s(t00, t11);

      t00 = mpqc::fenced_now(world);
      if (!K_new.is_initialized()) {
        K_new("mu, nu") = F("Y, mu, rho") * Q_ket("Y, nu, rho");
      } else {
        K_new("mu, nu") += F("Y, mu, rho") * Q_ket("Y, nu, rho");
      }
      t11 = mpqc::fenced_now(world);
      t_Kpart1_new += mpqc::duration_in_s(t00, t11);




    }
    auto t1_new_k = mpqc::fenced_now(world);
    auto t_tot_new = mpqc::duration_in_s(t0_new_k, t1_new_k);
    if (print_detail_) {
      ExEnv::out0() << "\nCADF-K new time decomposition:\n"
                    << "\tC block:                  " << t_Cblock << " s\n"
                    << "\tQ_bra:                    " << t_Qbra_new << " s\n"
                    << "\tQ_ket:                    " << t_Qket_new << " s\n"
                    << "\tE_ket:                    " << t_Eket_new << " s\n"
                    << "\tF = E - C M:              " << t_F_new << " s\n"
                    << "\tK_part1 = F Qket:         " << t_Kpart1_new << " s\n"
                    << "\nTotal K builder new time: " << t_tot_new << " s" << std::endl;
    }

    return K;
  }

  auto force_normsp(TA::Tensor<float> const &in,
                    TA::Tensor<float> const &F_norms) {
    auto const &irange = in.range();
    auto iext = irange.extent_data();

    // Q("X, mu, sigma")
    const auto X_size = iext[0];
    const auto mu_size = iext[1];
    const auto sig_size = iext[2];

    using SigPair = std::pair<int64_t, int64_t>;
    std::unordered_set<SigPair, boost::hash<SigPair>> Xsig;
    Xsig.reserve(X_size * sig_size);

    for (auto X = 0ul; X < X_size; ++X) {
      for (auto sig = 0ul; sig < sig_size; ++sig) {
        for (auto mu = 0ul; mu < mu_size; ++mu) {
          const auto val = in(X, mu, sig);
          if (val > force_shape_threshold_) {
            SigPair Xs_pair(X, sig);
            Xsig.insert(Xs_pair);
            break;
          }
        }
      }
    }

    auto const &orange = F_norms.range();
    TA::Tensor<float> out(orange, 0.0);

    // F("X, nu, sigma")
    auto oext = orange.extent_data();
    assert(X_size == oext[0]);
    assert(sig_size == oext[2]);
    const auto nu_size = oext[1];

    for (auto nu = 0ul; nu < nu_size; ++nu) {
      for (auto const &Xspair : Xsig) {
        const auto X = Xspair.first;
        const auto sig = Xspair.second;
        out(X, nu, sig) = std::numeric_limits<float>::max();
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

  array_type compute_Q_bra(const array_type &C, const array_type &D) {
    auto &world = C.world();

    array_type C_repl, D_repl;
    C_repl("X, mu, rho") = C("X, mu, rho");
    D_repl("rho, sig") = D("rho, sig");
    C_repl.make_replicated();
    D_repl.make_replicated();

    const auto &C_norms = C_repl.shape().data();

    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::detail::direct_3D_idx;

    Basis bs1((Basis()));
    for (auto RJ = 0; RJ != RJ_size_; ++RJ) {
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      auto tmp_bs = shift_basis_origin(*obs_, vec_RJ, RD_max_, dcell_);
      bs1 = lcao::gaussian::merge(bs1, *tmp_bs);
    }

    auto trange = lcao::gaussian::detail::create_trange(
        lcao::gaussian::BasisVector{{*X_dfbs_, *obs_, bs1}});
    auto tvolume = trange.tiles_range().volume();
    auto pmap = Policy::default_pmap(world, tvolume);

    // force Q norms
    auto ntiles_df = dfbs_->nclusters();
    auto ntiles_X_df = X_dfbs_->nclusters();
    auto ntiles0 = obs_->nclusters();
    auto ntiles1 = bs1.nclusters();
    auto nsig_per_RJ = ntiles0 * RD_size_;

    TA::Range range;
    range = TA::Range(std::array<int64_t, 3>{{ntiles_X_df, ntiles0, ntiles1}});
    TA::Tensor<float> norms(range, 0.0);
    for (auto X = 0ul; X != ntiles_X_df; ++X) {
      auto RX_ord = X / ntiles_df;
      auto RX_3D = direct_3D_idx(RX_ord, RJ_max_);

      for (auto sig = 0ul; sig != ntiles1; ++sig) {
        auto RJ_ord = sig / nsig_per_RJ;
        auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
        if (!(RX_3D == RJ_3D) && !(RX_3D == Vector3i{0, 0, 0})) continue;

        for (auto mu = 0ul; mu != ntiles0; ++mu) {
          auto significant = false;
          for (auto rho = 0ul; rho != ntiles0; ++rho) {
            auto val = C_norms(X, mu, RJ_ord * ntiles0 + rho);
            if (val >= force_shape_threshold_) {
              significant = true;
              break;
            }
          }
          if (significant)
            norms(X, mu, sig) = std::numeric_limits<float>::max();
        }
      }
    }

    TA::SparseShape<float> shape(world, norms, trange);

    array_type Q(world, trange, shape, pmap);

    auto create_task_Q_bra_tile = [&](
        array_type *Q, array_type *C_array, array_type *D_array, int64_t ord,
        int64_t RJ, std::array<int64_t, 3> tile_idx,
        std::array<size_t, 3> out_ext) {
      const auto tile_C_df = tile_idx[0];
      const auto tile_C_0 = tile_idx[1];
      const auto tile_D_1 = tile_idx[2];

      const auto ntiles_obs = obs_->nclusters();

      const auto D_0_stride = ntiles_obs * RD_size_;
      const auto C_0_stride = ntiles_obs * RJ_size_;
      const auto C_1_offset = (tile_C_df * ntiles_obs + tile_C_0) * C_0_stride;

      RowMatrixXd out_eig(out_ext[0] * out_ext[1], out_ext[2]);
      out_eig.setZero();

      for (auto tile_sum = 0; tile_sum != ntiles_obs; ++tile_sum) {
        const auto ord_D = tile_sum * D_0_stride + tile_D_1;
        const auto ord_C = C_1_offset + ntiles_obs * RJ + tile_sum;
        if (D_array->is_zero(ord_D) || C_array->is_zero(ord_C)) continue;

        Tile D = D_array->find(ord_D);
        Tile C = C_array->find(ord_C);

        const auto C_ext = C.range().extent();
        const auto D_ext = D.range().extent();

        assert(C_ext[2] == D_ext[0]);

        RowMatrixXd C_eig = TA::eigen_map(C, C_ext[0] * C_ext[1], C_ext[2]);
        RowMatrixXd D_eig = TA::eigen_map(D, D_ext[0], D_ext[1]);

        out_eig += C_eig * D_eig;
      }

      TA::Range out_range;
      out_range = TA::Range(
          std::array<size_t, 3>{{out_ext[0], out_ext[1], out_ext[2]}});
      TA::TensorD out_tile(out_range, 0.0);

      TA::eigen_map(out_tile, out_ext[0] * out_ext[1], out_ext[2]) = out_eig;

      Q->set(ord, out_tile);
    };

    for (auto X = 0; X != ntiles_X_df; ++X) {
      auto RX_ord = X / ntiles_df;
      auto RX_3D = direct_3D_idx(RX_ord, RJ_max_);
      const auto X_times_stride = X * ntiles0 * ntiles1;

      const auto &rng_X = trange.dim(0).tile(X);
      const size_t ext_X = rng_X.second - rng_X.first;

      for (auto sig = 0; sig != ntiles1; ++sig) {
        auto RJ_ord = sig / nsig_per_RJ;
        auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);

        if (!(RX_3D == RJ_3D) && !(RX_3D == Vector3i{0, 0, 0})) continue;

        const auto &rng_sig = trange.dim(2).tile(sig);
        const size_t ext_sig = rng_sig.second - rng_sig.first;

        for (auto mu = 0; mu != ntiles0; ++mu) {
          const auto tile_idx = X_times_stride + mu * ntiles1 + sig;
          if (Q.is_local(tile_idx) && !Q.is_zero(tile_idx)) {
            const auto &rng_mu = trange.dim(1).tile(mu);
            const size_t ext_mu = rng_mu.second - rng_mu.first;
            world.taskq.add(create_task_Q_bra_tile, &Q, &C_repl, &D_repl,
                            tile_idx, RJ_ord,
                            std::array<int64_t, 3>{{X, mu, sig % nsig_per_RJ}},
                            std::array<size_t, 3>{{ext_X, ext_mu, ext_sig}});
          }
        }
      }
    }
    world.gop.fence();
    Q.truncate();

    return Q;
  }

  array_type compute_Q_ket(const array_type &C, const array_type &D, const int64_t RJ_ord) {
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

    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };
    auto is_in_lattice_range = [](Vector3i const &in_idx, Vector3i const &range, Vector3i const &center = {0, 0, 0}) {
      if (in_idx(0) <= center(0) + range(0) && in_idx(0) >= center(0) - range(0) &&
          in_idx(1) <= center(1) + range(1) && in_idx(1) >= center(1) - range(1) &&
          in_idx(2) <= center(2) + range(2) && in_idx(2) >= center(2) - range(2))
        return true;
      else
        return false;
    };

    auto unshifted_Y_latt_range =
        max_latt_range(R_max_, RJ_max_ + RD_max_);

    auto shifted_Y_latt_range = RJ_max_;

    auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
    auto vec_RJ = direct_vector(RJ_ord, RJ_max_, dcell_);
    auto bs0 = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);
    auto bs1 = shift_basis_origin(*obs_, vec_RJ);

    auto trange = lcao::gaussian::detail::create_trange(
        lcao::gaussian::BasisVector{{*Y_dfbs_, *bs0, *bs1}});
    auto tvolume = trange.tiles_range().volume();
    auto pmap = Policy::default_pmap(world, tvolume);

    // force Q norms
    auto ntiles_df = dfbs_->nclusters();
    auto ntiles_obs = obs_->nclusters();

    auto ntilesY = Y_dfbs_->nclusters();
    auto ntiles0 = bs0->nclusters();
    auto ntiles1 = bs1->nclusters();
    auto ntiles_sum = ntiles_obs * RD_size_;

    TA::Range range;
    range = TA::Range(std::array<int64_t, 3>{{ntilesY, ntiles0, ntiles1}});
    TA::Tensor<float> norms(range, 0.0);

    using SigPair = std::pair<int64_t, int64_t>;
    using TilePairList = std::unordered_map<int64_t, std::vector<SigPair>>;
    TilePairList significant_Y_rho;

    for (auto nu = 0; nu != ntiles0; ++nu) {
      auto R_ord = nu / ntiles_obs;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_obs;
      significant_Y_rho.insert(std::make_pair(nu, std::vector<SigPair>()));

      for (auto Y = 0; Y != ntilesY; ++Y) {
        auto RY_ord = Y / ntiles_df;
        auto RY_3D = direct_3D_idx(RY_ord, unshifted_Y_latt_range);

        auto RYmR_3D = RY_3D - R_3D;
        if (!is_in_lattice_range(RYmR_3D, shifted_Y_latt_range)) continue;

        auto RYmR_ord = direct_ord_idx(RYmR_3D, shifted_Y_latt_range);
        auto Y_in_C = Y % ntiles_df + RYmR_ord * ntiles_df;

        for (auto rho = 0; rho != ntiles1; ++rho) {

          for (auto sig = 0; sig != ntiles_sum; ++sig) {
            auto RD_ord = sig / ntiles_obs;
            auto RJpRD_3D = RJ_3D + direct_3D_idx(RD_ord, RD_max_);

            if (!(RY_3D == R_3D) && !(RY_3D == RJpRD_3D)) continue;

            auto RJpRDmR_3D = RJpRD_3D - R_3D;
            if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

            auto RJpRDmR_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);
            if (std::find(RJ_list_.begin(), RJ_list_.end(), RJpRDmR_ord) == RJ_list_.end()) continue;

            auto sig_in_C = sig % ntiles_obs + RJpRDmR_ord * ntiles_obs;
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

    auto create_task_Q_ket_tile = [&](
        array_type *Q, array_type *C_array, array_type *D_array, int64_t ord,
        Vector3i RJmR, std::array<int64_t, 3> external_tile_idx,
        std::array<size_t, 3> external_offset,
        std::array<size_t, 3> external_extent) {
      const auto tile_C_df = external_tile_idx[0];
      const auto tile_C_0 = external_tile_idx[1];
      const auto tile_D_0 = external_tile_idx[2];

      const auto offset_df = external_offset[0];
      const auto offset_0 = external_offset[1];
      const auto offset_1 = external_offset[2];

      const auto ext_df = external_extent[0];
      const auto ext_0 = external_extent[1];
      const auto ext_1 = external_extent[2];

      const auto C_0_stride = ntiles_obs * RJ_size_;
      const auto C_df_stride = ntiles_obs * C_0_stride;
      const auto C_1_offset = tile_C_df * C_df_stride + tile_C_0 * C_0_stride;

      RowMatrixXd out_eig(ext_df * ext_0, ext_1);
      out_eig.setZero();

      for (auto tile_sum = 0; tile_sum != ntiles_sum; ++tile_sum) {
        const auto ord_D = tile_D_0 * ntiles_sum + tile_sum;
        auto RD_ord = tile_sum / ntiles_obs;
        auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
        auto RJpRDmR_3D = RJmR + RD_3D;
        if (!is_in_lattice_range(RJpRDmR_3D, RJ_max_)) continue;

        auto RJmRpRD_ord = direct_ord_idx(RJpRDmR_3D, RJ_max_);

        const auto ord_C =
            C_1_offset + RJmRpRD_ord * ntiles_obs + tile_sum % ntiles_obs;
        if (D_array->is_zero(ord_D) || C_array->is_zero(ord_C)) continue;

        Tile D = D_array->find(ord_D);
        Tile C = C_array->find(ord_C);

        const auto C_ext = C.range().extent();
        const auto D_ext = D.range().extent();

        assert(C_ext[2] == D_ext[1]);

        RowMatrixXd C_eig = TA::eigen_map(C, C_ext[0] * C_ext[1], C_ext[2]);
        RowMatrixXd D_eig = TA::eigen_map(D, D_ext[0], D_ext[1]);

        out_eig += C_eig * (D_eig.transpose());
      }

      TA::Range out_range;
      std::array<size_t, 3> lb{{offset_df, offset_0, offset_1}};
      std::array<size_t, 3> ub{
          {offset_df + ext_df, offset_0 + ext_0, offset_1 + ext_1}};
      out_range = TA::Range(lb, ub);

      TA::TensorD out_tile(out_range, 0.0);

      TA::eigen_map(out_tile, external_extent[0] * external_extent[1],
                    external_extent[2]) = out_eig;

      Q->set(ord, out_tile);
    };

    const auto Y_stride = ntiles0 * ntiles1;
    for (auto nu = 0; nu != ntiles0; ++nu) {
      auto R_ord = nu / ntiles_obs;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_obs;

      const auto &rng_nu = trange.dim(1).tile(nu);
      const auto offset_nu = rng_nu.first;
      const auto ext_nu = rng_nu.second - rng_nu.first;

      const auto nu_times_stride = nu * ntiles1;
      for (auto &significant_pair : significant_Y_rho[nu]) {
        auto Y = significant_pair.first;
        auto rho = significant_pair.second;

        const auto tile_idx = Y * Y_stride + nu_times_stride + rho;

        if (Q.is_local(tile_idx) && !Q.is_zero(tile_idx)) {
          auto RY_ord = Y / ntiles_df;
          auto RY_3D = direct_3D_idx(RY_ord, unshifted_Y_latt_range);
          auto Y_in_C =
              Y % ntiles_df +
              direct_ord_idx(RY_3D - R_3D, shifted_Y_latt_range) * ntiles_df;

          const auto &rng_Y = trange.dim(0).tile(Y);
          const auto offset_Y = rng_Y.first;
          const auto ext_Y = rng_Y.second - rng_Y.first;

          const auto &rng_rho = trange.dim(2).tile(rho);
          const auto offset_rho = rng_rho.first;
          const auto ext_rho = rng_rho.second - rng_rho.first;

          world.taskq.add(
              create_task_Q_ket_tile, &Q, &C_repl, &D_repl, tile_idx,
              RJ_3D - R_3D, std::array<int64_t, 3>{{Y_in_C, nu_in_C, rho}},
              std::array<size_t, 3>{{offset_Y, offset_nu, offset_rho}},
              std::array<size_t, 3>{{ext_Y, ext_nu, ext_rho}});
        }
      }
    }

    world.gop.fence();
    Q.truncate();

    return Q;
  }

  array_type compute_Q_ket(const array_type &C, const array_type &D) {
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

    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };

    auto unshifted_sig_latt_range = RJ_max_ + RD_max_;
    auto unshifted_Y_latt_range =
        max_latt_range(R_max_, unshifted_sig_latt_range);
    auto unshifted_Y_dfbs = shift_basis_origin(*dfbs_, zero_shift_base,
                                               unshifted_Y_latt_range, dcell_);

    auto shifted_sig_latt_range = RJ_max_ + RD_max_ + R_max_;
    auto shifted_Y_latt_range = shifted_sig_latt_range;

    auto bs0 = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);
    auto bs1 = shift_basis_origin(*obs_, zero_shift_base, RJ_max_, dcell_);

    auto trange = lcao::gaussian::detail::create_trange(
        lcao::gaussian::BasisVector{{*unshifted_Y_dfbs, *bs0, *bs1}});
    auto tvolume = trange.tiles_range().volume();
    auto pmap = Policy::default_pmap(world, tvolume);

    // force Q norms
    auto ntiles_df = dfbs_->nclusters();
    auto ntiles_obs = obs_->nclusters();

    auto ntiles_Y_df = unshifted_Y_dfbs->nclusters();
    auto ntiles0 = bs0->nclusters();
    auto ntiles1 = bs1->nclusters();
    auto ntiles_per_RJ = ntiles_obs * RD_size_;

    TA::Range range;
    range = TA::Range(std::array<int64_t, 3>{{ntiles_Y_df, ntiles0, ntiles1}});
    TA::Tensor<float> norms(range, 0.0);

    using SigPair = std::pair<int64_t, int64_t>;
    using TilePairList = std::unordered_map<int64_t, std::vector<SigPair>>;
    TilePairList significant_Y_rho;

    for (auto nu = 0; nu != ntiles0; ++nu) {
      auto R_ord = nu / ntiles_obs;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_obs;
      significant_Y_rho.insert(std::make_pair(nu, std::vector<SigPair>()));

      for (auto Y = 0; Y != ntiles_Y_df; ++Y) {
        auto RY_ord = Y / ntiles_df;
        auto RY_3D = direct_3D_idx(RY_ord, unshifted_Y_latt_range);

        auto RYmR_ord = direct_ord_idx(RY_3D - R_3D, shifted_Y_latt_range);
        auto Y_in_C = Y % ntiles_df + RYmR_ord * ntiles_df;

        for (auto rho = 0; rho != ntiles1; ++rho) {
          auto RJ_ord = rho / ntiles_obs;
          auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);

          for (auto sig = 0; sig != ntiles_per_RJ; ++sig) {
            auto RD_ord = sig / ntiles_obs;
            auto RJpRD_3D = RJ_3D + direct_3D_idx(RD_ord, RD_max_);

            if (!(RY_3D == R_3D) && !(RY_3D == RJpRD_3D)) continue;

            auto RJpRDmR_ord =
                direct_ord_idx(RJpRD_3D - R_3D, shifted_sig_latt_range);
            auto sig_in_C = sig % ntiles_obs + RJpRDmR_ord * ntiles_obs;
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

    auto create_task_Q_ket_tile = [&](
        array_type *Q, array_type *C_array, array_type *D_array, int64_t ord,
        Vector3i RJmR, std::array<int64_t, 3> external_tile_idx,
        std::array<size_t, 3> external_offset,
        std::array<size_t, 3> external_extent) {
      const auto tile_C_df = external_tile_idx[0];
      const auto tile_C_0 = external_tile_idx[1];
      const auto tile_D_0 = external_tile_idx[2];

      const auto offset_df = external_offset[0];
      const auto offset_0 = external_offset[1];
      const auto offset_1 = external_offset[2];

      const auto ext_df = external_extent[0];
      const auto ext_0 = external_extent[1];
      const auto ext_1 = external_extent[2];

      const auto ntiles_obs = obs_->nclusters();
      const auto D_0_stride = ntiles_obs * RD_size_;

      const auto shifted_sig_latt_range = RJ_max_ + RD_max_ + R_max_;

      const auto shifted_sig_latt_size =
          1 + direct_ord_idx(shifted_sig_latt_range, shifted_sig_latt_range);
      const auto C_0_stride = ntiles_obs * shifted_sig_latt_size;
      const auto C_df_stride = ntiles_obs * C_0_stride;
      const auto C_1_offset = tile_C_df * C_df_stride + tile_C_0 * C_0_stride;

      RowMatrixXd out_eig(ext_df * ext_0, ext_1);
      out_eig.setZero();

      for (auto tile_sum = 0; tile_sum != D_0_stride; ++tile_sum) {
        const auto ord_D = tile_D_0 * D_0_stride + tile_sum;
        auto RD_ord = tile_sum / ntiles_obs;
        auto RD_3D = direct_3D_idx(RD_ord, RD_max_);
        auto RJmRpRD_ord = direct_ord_idx(RJmR + RD_3D, shifted_sig_latt_range);

        const auto ord_C =
            C_1_offset + RJmRpRD_ord * ntiles_obs + tile_sum % ntiles_obs;
        if (D_array->is_zero(ord_D) || C_array->is_zero(ord_C)) continue;

        Tile D = D_array->find(ord_D);
        Tile C = C_array->find(ord_C);

        const auto C_ext = C.range().extent();
        const auto D_ext = D.range().extent();

        assert(C_ext[2] == D_ext[1]);

        RowMatrixXd C_eig = TA::eigen_map(C, C_ext[0] * C_ext[1], C_ext[2]);
        RowMatrixXd D_eig = TA::eigen_map(D, D_ext[0], D_ext[1]);

        out_eig += C_eig * (D_eig.transpose());
      }

      TA::Range out_range;
      std::array<size_t, 3> lb{{offset_df, offset_0, offset_1}};
      std::array<size_t, 3> ub{
          {offset_df + ext_df, offset_0 + ext_0, offset_1 + ext_1}};
      out_range = TA::Range(lb, ub);

      TA::TensorD out_tile(out_range, 0.0);

      TA::eigen_map(out_tile, external_extent[0] * external_extent[1],
                    external_extent[2]) = out_eig;

      Q->set(ord, out_tile);
    };

    const auto Y_stride = ntiles0 * ntiles1;
    for (auto nu = 0; nu != ntiles0; ++nu) {
      auto R_ord = nu / ntiles_obs;
      auto R_3D = direct_3D_idx(R_ord, R_max_);
      auto nu_in_C = nu % ntiles_obs;

      const auto &rng_nu = trange.dim(1).tile(nu);
      const auto offset_nu = rng_nu.first;
      const auto ext_nu = rng_nu.second - rng_nu.first;

      const auto nu_times_stride = nu * ntiles1;
      for (auto &significant_pair : significant_Y_rho[nu]) {
        auto Y = significant_pair.first;
        auto rho = significant_pair.second;

        const auto tile_idx = Y * Y_stride + nu_times_stride + rho;

        if (Q.is_local(tile_idx) && !Q.is_zero(tile_idx)) {
          auto RY_ord = Y / ntiles_df;
          auto RY_3D = direct_3D_idx(RY_ord, unshifted_Y_latt_range);
          auto Y_in_C =
              Y % ntiles_df +
              direct_ord_idx(RY_3D - R_3D, shifted_Y_latt_range) * ntiles_df;

          auto RJ_ord = rho / ntiles_obs;
          auto RJ_3D = direct_3D_idx(RJ_ord, RJ_max_);
          auto rho_in_D = rho % ntiles_obs;

          const auto &rng_Y = trange.dim(0).tile(Y);
          const auto offset_Y = rng_Y.first;
          const auto ext_Y = rng_Y.second - rng_Y.first;

          const auto &rng_rho = trange.dim(2).tile(rho);
          const auto offset_rho = rng_rho.first;
          const auto ext_rho = rng_rho.second - rng_rho.first;

          world.taskq.add(
              create_task_Q_ket_tile, &Q, &C_repl, &D_repl, tile_idx,
              RJ_3D - R_3D, std::array<int64_t, 3>{{Y_in_C, nu_in_C, rho_in_D}},
              std::array<size_t, 3>{{offset_Y, offset_nu, offset_rho}},
              std::array<size_t, 3>{{ext_Y, ext_nu, ext_rho}});
        }
      }
    }

    world.gop.fence();
    Q.truncate();

    return Q;
  }

  array_type compute_contr_EQ(array_type const &Q, double target_precision) {

  }

  TA::Tensor<float> force_norms(TA::Tensor<float> const &in_norms,
                                int64_t const in_latt_range_size,
                                int64_t const out_latt_range_size) {
    auto &in_range = in_norms.range();
    auto in_extent = in_range.extent_data();

    TA::Range out_range;
    out_range = TA::Range(std::array<size_t, 3>{
        {in_extent[0], in_extent[1] / in_latt_range_size * out_latt_range_size,
         in_extent[2]}});
    auto out_extent = out_range.extent_data();
    TA::Tensor<float> out_norms(out_range, 0.0);

    std::unordered_set<int64_t> sig_sigma;
    std::unordered_set<int64_t> sig_X;

    sig_sigma.clear();
    sig_X.clear();
    for (auto X = 0ul; X < in_extent[0]; ++X) {
      // For every ν, we determine save a list of important X and sigma
      for (auto sigma = 0ul; sigma < in_extent[2]; ++sigma) {
        for (auto mu = 0ul; mu < in_extent[1]; ++mu) {
          const auto val = in_norms(X, mu, sigma);

          if (val >= force_shape_threshold_) {
            sig_X.insert(X);
            sig_sigma.insert(sigma);
            break;
          }
        }
      }
    }

    for (auto nu = 0ul; nu < out_extent[1]; ++nu) {
      // Then for every X we mark all important sigma. The output will have
      // significant X sigma pairs that the input did not.
      for (auto const &X : sig_X) {
        for (auto const &sigma : sig_sigma) {
          out_norms(X, nu, sigma) = std::numeric_limits<float>::max();
        }
      }
    }

    return out_norms;
  }

  norm_type<3> compute_eri3_norm(
      TA::Tensor<float> const &in, size_t const &natoms_per_uc,
      double thresh = 1.0e-12,
      Vector3i const &lattice_range0 = Vector3i({0, 0, 0}),
      Vector3i const &lattice_range1 = Vector3i({0, 0, 0}),
      Vector3i const &lattice_range_df = Vector3i({0, 0, 0}),
      Vector3i const &lattice_center0 = Vector3i({0, 0, 0}),
      Vector3i const &lattice_center1 = Vector3i({0, 0, 0}),
      Vector3i const &lattice_center_df = Vector3i({0, 0, 0})) {
    auto const &range = in.range();
    auto ext = range.extent_data();

    norm_type<3> norms;
    norms.reserve(ext[1] * ext[2]);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;

    using idx_type = std::array<int, 3>;
    auto val_max = std::numeric_limits<float>::max();
    for (auto a = 0; a < ext[1]; ++a) {
      const auto uc0_ord = a / natoms_per_uc;
      const auto uc0_3D =
          direct_3D_idx(uc0_ord, lattice_range0) + lattice_center0;
      const auto uc0_ord_in_df =
          direct_ord_idx(uc0_3D - lattice_center_df, lattice_range_df);
      const auto a_in_df = uc0_ord_in_df * natoms_per_uc + a % natoms_per_uc;

      for (auto b = 0; b < ext[2]; ++b) {
        auto uc1_ord = b / natoms_per_uc;
        auto uc1_3D = direct_3D_idx(uc1_ord, lattice_range1) + lattice_center1;
        auto uc1_ord_in_df =
            direct_ord_idx(uc1_3D - lattice_center_df, lattice_range_df);
        auto b_in_df = uc1_ord_in_df * natoms_per_uc + b % natoms_per_uc;

        auto in_val = std::max(in(a_in_df, a, b), in(b_in_df, a, b));

        if (in_val >= thresh) {
          norms.emplace_back(
              std::make_pair(idx_type{{int(a_in_df), a, b}}, val_max));
          norms.emplace_back(
              std::make_pair(idx_type{{int(b_in_df), a, b}}, val_max));
        }
      }
    }

    return norms;
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

};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
