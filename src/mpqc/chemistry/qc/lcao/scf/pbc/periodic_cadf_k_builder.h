#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_

#include "mpqc/chemistry/qc/lcao/integrals/density_fitting/cadf_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_three_center_contraction_builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicCADFKBuilder {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using DirectTArray = typename Factory::DirectTArray;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;

  template <int rank>
  using norm_type = std::vector<std::pair<std::array<int, rank>, float>>;

  PeriodicCADFKBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // by-cluster orbital basis and df basis
    auto obs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    auto dfbs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));

    auto dcell = ao_factory_.unitcell().dcell();
    auto R_max = ao_factory_.R_max();
    auto RJ_max = ao_factory_.RJ_max();
    auto RD_max = ao_factory_.RD_max();
    auto R_size = ao_factory_.R_size();
    auto RJ_size = ao_factory_.RJ_size();
    auto RD_size = ao_factory_.RD_size();

    const Vector3i ref_latt_range = {0, 0, 0};
    const Vector3i ref_latt_center = {0, 0, 0};
    const auto natoms_per_uc = ao_factory_.unitcell().natoms();
    Vector3d zero_shift_base(0.0, 0.0, 0.0);

    using ::mpqc::lcao::detail::direct_3D_idx;
    using ::mpqc::lcao::detail::direct_ord_idx;
    using ::mpqc::lcao::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

    // compute C(X, μ_0, ρ_Rj)
    auto X_dfbs = shift_basis_origin(*dfbs, zero_shift_base, RJ_max, dcell);
    {
      if (C_bra_.empty())
        C_bra_ = std::vector<array_type>(RJ_size, array_type());

      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        array_type &C = C_bra_[RJ];
        auto RJ_3D = direct_3D_idx(RJ, RJ_max);
        auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
        auto bs1 = shift_basis_origin(*obs, vec_RJ);
        C = lcao::cadf_fitting_coefficients<Tile, Policy>(
            world, *obs, *bs1, *X_dfbs, natoms_per_uc, ref_latt_range,
            ref_latt_range, RJ_max, ref_latt_center, RJ_3D, ref_latt_center);
      }
    }

    //    auto bs1 = shift_basis_origin(*obs, zero_shift_base, RJ_max, dcell);
    //    auto dfbs1 = shift_basis_origin(*dfbs, zero_shift_base, RJ_max,
    //    dcell);
    //    C_ = lcao::cadf_fitting_coefficients<Tile, Policy>(
    //        world, *obs, *bs1, *dfbs1, natoms_per_uc, ref_latt_range, RJ_max,
    //        RJ_max);
    //    ExEnv::out0() << "\nC_ = \n" << C_ << std::endl;

    auto max_latt_range = [](Vector3i const &l, Vector3i const &r) {
      auto x = std::max(l(0), r(0));
      auto y = std::max(l(1), r(1));
      auto z = std::max(l(2), r(2));
      return Vector3i({x, y, z});
    };

    // compute C(Y, ν_R, σ_(Rj+Rd))
    auto Y_latt_range = max_latt_range(R_max, RJ_max + RD_max);
    auto Y_dfbs =
        shift_basis_origin(*dfbs, zero_shift_base, Y_latt_range, dcell);
    {
      if (C_ket_.empty())
        C_ket_ = std::vector<array_type>(RJ_size, array_type());

      auto bs0 = shift_basis_origin(*obs, zero_shift_base, R_max, dcell);
      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        array_type &C = C_ket_[RJ];
        auto RJ_3D = direct_3D_idx(RJ, RJ_max);
        auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
        auto bs1 = shift_basis_origin(*obs, vec_RJ, RD_max, dcell);
        C = lcao::cadf_fitting_coefficients<Tile, Policy>(
            world, *bs0, *bs1, *Y_dfbs, natoms_per_uc, R_max, RD_max,
            Y_latt_range, ref_latt_center, RJ_3D, ref_latt_center);
      }
    }

    // compute M(X, Y)
    {
      const auto dfbs_vector = lcao::gaussian::BasisVector{{*X_dfbs, *Y_dfbs}};
      auto engine = lcao::gaussian::make_engine_pool(
          libint2::Operator::coulomb,
          utility::make_array_of_refs(*X_dfbs, *Y_dfbs),
          libint2::BraKet::xs_xs);
      M_ = lcao::gaussian::sparse_integrals(world, engine, dfbs_vector);
    }

    const auto screen_thresh = ao_factory_.screen_threshold();
    // compute E(X, ν_R, σ_(Rj+Rd))
    {
      if (E_bra_.empty())
        E_bra_ = std::vector<DirectTArray>(RJ_size, DirectTArray());

      auto bs0 = shift_basis_origin(*obs, zero_shift_base, R_max, dcell);

      for (auto RJ = 0; RJ != RJ_size; ++RJ) {
        DirectTArray &E = E_bra_[RJ];
        if (!E.array().is_initialized()) {
          auto vec_RJ = direct_vector(RJ, RJ_max, dcell);
          auto bs1 = shift_basis_origin(*obs, vec_RJ, RD_max, dcell);

          auto bs_array = utility::make_array_of_refs(*X_dfbs, *bs0, *bs1);
          auto bs_vector = lcao::gaussian::BasisVector{{*X_dfbs, *bs0, *bs1}};

          auto screen_eng = lcao::gaussian::make_engine_pool(
              libint2::Operator::coulomb, bs_array, libint2::BraKet::xx_xx);
          auto screener = std::make_shared<lcao::gaussian::SchwarzScreen>(
              lcao::gaussian::create_schwarz_screener(
                  world, screen_eng, bs_vector, screen_thresh));

          auto engine = lcao::gaussian::make_engine_pool(
              libint2::Operator::coulomb, bs_array, libint2::BraKet::xs_xx);

          E = lcao::gaussian::direct_sparse_integrals(world, engine, bs_vector,
                                                      std::move(screener));
        }
      }
    }
  }

  ~PeriodicCADFKBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    array_type K;
    // to be implemented
    K("mu, nu") = 1.0 * D("mu, nu");
    return K;
  }

 private:
  Factory &ao_factory_;
  bool print_detail_;
  std::unique_ptr<PTC_Builder> three_center_builder_;

  //  array_type C_;
  std::vector<array_type> C_bra_;
  std::vector<array_type> C_ket_;
  std::vector<array_type> F_;
  array_type M_;
  std::vector<DirectTArray> E_bra_;
  std::vector<DirectTArray> E_ket_;

 private:
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
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
