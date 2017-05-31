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
class PeriodicCADFKBuilder {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using DirectTArray = typename Factory::DirectTArray;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;
  using Qmatrix = ::mpqc::lcao::gaussian::Qmatrix;
  using Basis = ::mpqc::lcao::gaussian::Basis;

  template <int rank>
  using norm_type = std::vector<std::pair<std::array<int, rank>, float>>;

  PeriodicCADFKBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();
    print_detail_ = ao_factory_.print_detail();
    force_shape_threshold_ = 1.0e-12;
    mpqc::time_point t0, t1;

    // by-cluster orbital basis and df basis
    auto obs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    auto dfbs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));

    auto dcell = ao_factory_.unitcell().dcell();
    auto R_max = ao_factory_.R_max();
    auto RJ_max = R_max;  // replace RJ range with R range
    auto RD_max = ao_factory_.RD_max();
    auto R_size = ao_factory_.R_size();
    auto RJ_size = R_size;  // replace RJ range with R range

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
    auto X_dfbs = shift_basis_origin(*dfbs, zero_shift_base, RJ_max, dcell);
    {
      C_bra_ = std::vector<array_type>(RJ_size, array_type());

      ExEnv::out0() << "\nComputing C_bra ...\n";
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*X_dfbs);

      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        array_type &C = C_bra_[RJ];
        auto RJ_3D = direct_3D_idx(RJ, RJ_max);
        auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
        auto bs1 = shift_basis_origin(*obs, vec_RJ);
        C = lcao::cadf_fitting_coefficients<Tile, Policy>(
            M, *obs, *bs1, *X_dfbs, natoms_per_uc, ref_latt_range,
            ref_latt_range, RJ_max, ref_latt_center, RJ_3D, ref_latt_center);
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_C_bra = mpqc::duration_in_s(t0, t1);

    // test new C_bra
    {
      auto t0 = mpqc::fenced_now(world);
      ExEnv::out0() << "\nComputing new C_bra ...\n";
      auto X_dfbs = shift_basis_origin(*dfbs, zero_shift_base, RJ_max, dcell);
      auto bs1 = shift_basis_origin(*obs, zero_shift_base, RJ_max, dcell);

      const auto by_atom_dfbs = lcao::detail::by_center_basis(*X_dfbs);
      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

      array_type C = lcao::cadf_fitting_coefficients<Tile, Policy>(
          M, *obs, *bs1, *X_dfbs, natoms_per_uc, ref_latt_range, RJ_max,
          RJ_max);

      detail::print_size_info(C, "C bra");
      auto t1 = mpqc::fenced_now(world);
      auto dur = mpqc::duration_in_s(t0, t1);
      ExEnv::out0() << "\tnew C_bra:            " << dur << " s\n" << std::endl;

      dump_shape_0_12(C, "C_bra");
    }

    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };

    // compute C(Y, ν_R, σ_(Rj+Rd))
    t0 = mpqc::fenced_now(world);
    auto Y_latt_range = max_latt_range(R_max, RJ_max + RD_max);
    auto Y_dfbs =
        shift_basis_origin(*dfbs, zero_shift_base, Y_latt_range, dcell);
    {
      C_ket_ = std::vector<array_type>(RJ_size, array_type());

      ExEnv::out0() << "\nComputing C_ket ...\n";
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*Y_dfbs);

      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);
      auto bs0 = shift_basis_origin(*obs, zero_shift_base, R_max, dcell);
      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        array_type &C = C_ket_[RJ];
        auto RJ_3D = direct_3D_idx(RJ, RJ_max);
        auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
        auto bs1 = shift_basis_origin(*obs, vec_RJ, RD_max, dcell);
        C = lcao::cadf_fitting_coefficients<Tile, Policy>(
            M, *bs0, *bs1, *Y_dfbs, natoms_per_uc, R_max, RD_max, Y_latt_range,
            ref_latt_center, RJ_3D, ref_latt_center);

//        ExEnv::out0() << "\nRJ = " << RJ << ":\n" << std::endl;
//        detail::print_size_info(C, "C ket");
      }
    }
    t1 = mpqc::fenced_now(world);
    double t_C_ket = mpqc::duration_in_s(t0, t1);

    // test new C_ket
    {
      auto t0 = mpqc::fenced_now(world);
      ExEnv::out0() << "\nComputing new C_ket ...\n";
      auto max_range = R_max + RJ_max + RD_max;
      auto bs1 = shift_basis_origin(*obs, zero_shift_base, max_range, dcell);

      auto Y_dfbs =
          shift_basis_origin(*dfbs, zero_shift_base, max_range, dcell);
      const auto by_atom_dfbs = lcao::detail::by_center_basis(*Y_dfbs);
      auto M = compute_eri2(world, by_atom_dfbs, by_atom_dfbs);

      array_type C = lcao::cadf_fitting_coefficients<Tile, Policy>(
          M, *obs, *bs1, *Y_dfbs, natoms_per_uc, ref_latt_range, max_range,
          max_range);

      detail::print_size_info(C, "C ket");
      auto t1 = mpqc::fenced_now(world);

      auto dur = mpqc::duration_in_s(t0, t1);
      ExEnv::out0() << "\tnew C_ket:            " << dur << " s\n" << std::endl;

      dump_shape_0_12(C, "C_ket");
    }

    // compute M(X, Y)
    t0 = mpqc::fenced_now(world);
    M_ = compute_eri2(world, *X_dfbs, *Y_dfbs);
    t1 = mpqc::fenced_now(world);
    double t_M = mpqc::duration_in_s(t0, t1);

    auto oper_type = libint2::Operator::coulomb;
    const auto screen_thresh = ao_factory_.screen_threshold();
    auto screen_norm_op = ::mpqc::lcao::gaussian::detail::l2Norm;
    auto screen_engine = make_engine_pool(
        oper_type, utility::make_array_of_refs(*dfbs, *obs, *obs),
        libint2::BraKet::xx_xx);

    // compute E(X, ν_R, σ_(Rj+Rd))
    t0 = mpqc::fenced_now(world);
    {
      if (E_bra_.empty())
        E_bra_ = std::vector<DirectTArray>(RJ_size, DirectTArray());

      auto bs0 = shift_basis_origin(*obs, zero_shift_base, R_max, dcell);

      std::shared_ptr<Qmatrix> Qbra, Qket;
      Qbra = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *X_dfbs, screen_norm_op));

      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        DirectTArray &E = E_bra_[RJ];
        if (!E.array().is_initialized()) {
          auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
          auto bs1 = shift_basis_origin(*obs, vec_RJ, RD_max, dcell);

          auto bs_array = utility::make_array_of_refs(*X_dfbs, *bs0, *bs1);
          auto bs_vector = lcao::gaussian::BasisVector{{*X_dfbs, *bs0, *bs1}};

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
    {
      auto t0 = mpqc::fenced_now(world);
      ExEnv::out0() << "\nComputing new E_bra ...\n";
      std::shared_ptr<Qmatrix> Qbra, Qket;
      Qbra = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *X_dfbs, screen_norm_op));

      auto max_range = R_max + RJ_max + RD_max;
      auto bs1 = shift_basis_origin(*obs, zero_shift_base, max_range, dcell);

      Qket = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *obs, *bs1, screen_norm_op));
      auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
          lcao::gaussian::SchwarzScreen(Qbra, Qket, screen_thresh));

      auto bs_array = utility::make_array_of_refs(*X_dfbs, *obs, *bs1);
      auto bs_vector = lcao::gaussian::BasisVector{{*X_dfbs, *obs, *bs1}};
      auto engine =
          make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);
      auto E = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                  std::move(screener));
      auto t1 = mpqc::fenced_now(world);
      auto dur = mpqc::duration_in_s(t0, t1);
      ExEnv::out0() << "\tnew E_bra:            " << dur << " s\n" << std::endl;
    }

    // compute E(Y, μ_0, ρ_Rj)
    t0 = mpqc::fenced_now(world);
    {
      if (E_ket_.empty())
        E_ket_ = std::vector<DirectTArray>(RJ_size, DirectTArray());

      std::shared_ptr<Qmatrix> Qbra, Qket;
      Qbra = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *Y_dfbs, screen_norm_op));

      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        DirectTArray &E = E_ket_[RJ];
        if (!E.array().is_initialized()) {
          auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
          auto bs1 = shift_basis_origin(*obs, vec_RJ);

          auto bs_array = utility::make_array_of_refs(*Y_dfbs, *obs, *bs1);
          auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs, *obs, *bs1}};

          Qket = std::make_shared<Qmatrix>(
              Qmatrix(world, screen_engine, *obs, *bs1, screen_norm_op));
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

    // test new ket
    {
      auto t0 = mpqc::fenced_now(world);
      ExEnv::out0() << "\nComputing new E_ket ...\n";
      std::shared_ptr<Qmatrix> Qbra, Qket;
      Qbra = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *Y_dfbs, screen_norm_op));

      auto bs1 = shift_basis_origin(*obs, zero_shift_base, RJ_max, dcell);
      Qket = std::make_shared<Qmatrix>(
          Qmatrix(world, screen_engine, *obs, *bs1, screen_norm_op));
      auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
          lcao::gaussian::SchwarzScreen(Qbra, Qket, screen_thresh));

      auto bs_array = utility::make_array_of_refs(*Y_dfbs, *obs, *bs1);
      auto bs_vector = lcao::gaussian::BasisVector{{*Y_dfbs, *obs, *bs1}};
      auto engine =
          make_engine_pool(oper_type, bs_array, libint2::BraKet::xs_xx);
      auto E = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                  std::move(screener));

      auto t1 = mpqc::fenced_now(world);
      auto dur = mpqc::duration_in_s(t0, t1);
      ExEnv::out0() << "\tnew E_ket:            " << dur << " s\n" << std::endl;
    }

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
  std::unique_ptr<PTC_Builder> three_center_builder_;

  std::vector<array_type> C_bra_;
  std::vector<array_type> C_ket_;
  std::vector<array_type> F_;
  array_type M_;
  std::vector<DirectTArray> E_bra_;
  std::vector<DirectTArray> E_ket_;

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
      const double norm = C("i,j").norm();
      if (norm == 0.0) {
        continue;
      }
      array_type &Q = Q_bra[RJ];
      Q("X, mu, sig") = C("X, mu, rho") * D("rho, sig");
    }
    t1 = mpqc::fenced_now(world);
    double t_Q_bra = mpqc::duration_in_s(t0, t1);

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
      //      detail::print_size_info(F, "F");
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
      const double norm = C("i,j").norm();
      if (norm == 0.0) {
        continue;
      }
      array_type &Q = Q_ket[RJ];
      Q("Y, nu, rho") = C("Y, nu, sig") * D("rho, sig");
    }
    t1 = mpqc::fenced_now(world);
    double t_Q_ket = mpqc::duration_in_s(t0, t1);

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

  void dump_shape_02_1(array_type const &Q, std::string const &name, int rj = 0) {
    std::cout << "Dumping shape for 3-d tensor (02-1) to file " << name << "_" << rj << ".csv ...\n";
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

    std::ofstream outfile(name + std::to_string(rj) + ".csv");
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

  void dump_shape_01_2(array_type const &C, std::string const &name, int rj = 0) {
    std::cout << "Dumping shape for 3-d tensor (01-2) to file " << name << "_" << rj << ".csv ...\n";
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

    std::ofstream outfile(name + std::to_string(rj) + ".csv");
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

  void dump_shape_0_12(array_type const &C, std::string const &name, int rj = 0) {
    std::cout << "Dumping shape for 3-d tensor (0-12) to file " << name << "_" << rj << ".csv ...\n";
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

    std::ofstream outfile(name + std::to_string(rj) + ".csv");
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

    std::ofstream outfile(name + std::to_string(rj) + ".csv");
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

};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
