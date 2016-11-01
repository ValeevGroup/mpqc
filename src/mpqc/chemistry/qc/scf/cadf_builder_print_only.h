#pragma once
#ifndef MPQC_SCF_CADFBUILDERPRINTONLY_H
#define MPQC_SCF_CADFBUILDERPRINTONLY_H

#include <tiledarray.h>


#include "mpqc/util/misc/time.h"
#include "mpqc/math/external/tiledarray/array_info.h"
#include "mpqc/chemistry/qc/scf/util.h"

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/tile.h"
#include "mpqc/math/tensor/clr/tensor_transforms.h"

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/minimize_storage.h"

#include <mpqc/chemistry/qc/scf/builder.h>
#include <mpqc/chemistry/qc/integrals/make_engine.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral.h>
#include <mpqc/chemistry/qc/scf/cadf_fitting_coeffs.h>
#include <mpqc/chemistry/qc/scf/cadf_helper_functions.h>

#include <mpqc/math/external/tiledarray/tensor_store.h>

#include <vector>
#include <iostream>
#include <unordered_set>

namespace mpqc {
namespace scf {

class PrintOnlyCADFFockBuilder : public FockBuilder {
 public:
  using TileType = TA::TensorD;
  using ArrayType = FockBuilder::array_type;

 private:
  ArrayType E_;  // <Κ |G| κ λ > three center two electron coulomb integrals
  ArrayType M_;  // <Κ |G| Λ > two center two electron coulomb integrals
  ArrayType Mchol_inv_;  // Chol(<Κ |G| Λ >)^-1
  ArrayType C_df_;       // CADF fitting coeffs

  bool use_forced_shape_ = false;
  float force_threshold_ = 0.0;
  double lcao_chop_threshold_ = 0.0;

  int iteration = 0;

 public:
  PrintOnlyCADFFockBuilder(
      Molecule const &clustered_mol,
      Molecule const &df_clustered_mol,
      basis::BasisSet const &obs_set, basis::BasisSet const &dfbs_set,
      integrals::AtomicIntegral<TileType, TA::SparsePolicy> &ao_ints,
      bool use_forced_shape, double force_threshold,
      double lcao_chop_threshold = 0.0)
      : PrintOnlyCADFFockBuilder(clustered_mol, df_clustered_mol, obs_set,
                                 dfbs_set, ao_ints) {
    use_forced_shape_ = use_forced_shape;
    force_threshold_ = force_threshold;
    lcao_chop_threshold_ = lcao_chop_threshold;
  }

  PrintOnlyCADFFockBuilder(
      Molecule const &clustered_mol,
      Molecule const &df_clustered_mol,
      basis::BasisSet const &obs_set, basis::BasisSet const &dfbs_set,
      integrals::AtomicIntegral<TileType, TA::SparsePolicy> &ao_ints)
      : FockBuilder() {
    // Grab needed ao integrals
    E_ = ao_ints.compute(L"( Κ | G|κ λ)");
    util::write_shape_tuple3D(
        E_, std::string("E_shape.txt"));
    M_ = ao_ints.compute(L"( Κ | G| Λ )");

    // Form L^{-1} for M
    auto M_eig = array_ops::array_to_eigen(M_);
    using MatType = decltype(M_eig);
    MatType L_inv_eig = MatType(Eigen::LLT<MatType>(M_eig).matrixL()).inverse();

    auto trange1_M = M_.trange().data()[0];  // Assumes symmetric blocking
    Mchol_inv_ = array_ops::eigen_to_array<TA::TensorD>(
        M_.world(), L_inv_eig, trange1_M, trange1_M);

    std::unordered_map<std::size_t, std::size_t> obs_atom_to_cluster_map;
    std::unordered_map<std::size_t, std::size_t> dfbs_atom_to_cluster_map;

    basis::Basis obs = ao_ints.orbital_basis_registry().retrieve(L"κ");
    basis::Basis dfbs = ao_ints.orbital_basis_registry().retrieve(L"Κ");

    auto eng_pool = integrals::make_engine_pool(
        libint2::Operator::coulomb, utility::make_array_of_refs(dfbs, dfbs),
        libint2::BraKet::xs_xs);

    ArrayType C_df_temp = scf::compute_atomic_fitting_coeffs(
        M_.world(), clustered_mol, df_clustered_mol, obs_set, dfbs_set,
        eng_pool, obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);

    auto by_cluster_trange =
        integrals::detail::create_trange(utility::make_array(dfbs, obs, obs));

    C_df_ =
        scf::reblock_from_atoms(C_df_temp, obs_atom_to_cluster_map,
                                dfbs_atom_to_cluster_map, by_cluster_trange);

    util::write_shape_tuple3D(
        C_df_, std::string("C_df_shape.txt"));
  }

  ~PrintOnlyCADFFockBuilder() = default;

  void register_fock(const TA::TSpArrayD &fock,
                     FormulaRegistry<TA::TSpArrayD> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)[df]"), fock);
  }

  ArrayType operator()(ArrayType const &D, ArrayType const &C) override {
    ++iteration;

    ArrayType E_mo;  // Temp array shared by J and K
    E_mo("X, i, mu") = E_("X, mu, nu") * C("nu, i");
    E_mo.truncate();

    ArrayType G;
    G("m, n") = 2 * compute_J(C, E_mo)("m, n") - compute_K(C, E_mo)("m, n");
    return G;
  }

  void print_iter(std::string const &leader) override {}

  rapidjson::Value results(rapidjson::Document &d) override {}

 private:
  ArrayType compute_J(ArrayType const &C, ArrayType const &E_mo) {
    ArrayType J;
    J("mu, nu") = E_("X, mu, nu") *
                  (Mchol_inv_("Z, X") *
                   (Mchol_inv_("Z, Y") * (E_mo("Y, i, rho") * C("rho, i"))));

    return J;
  }

  array_type compute_K(ArrayType const &C_in, ArrayType const &E_mo) {
    ArrayType L, K;        // Matrices
    ArrayType C_mo, F_df;  // Tensors

    // Deep copy C for chopping
    ArrayType C;
    C("mu, i") = C_in("mu, i");

    if (lcao_chop_threshold_ != 0.0) {
      TA::foreach_inplace(C, [&](TA::Tensor<double> &t) {
        const auto norm = t.norm();
        if (norm > lcao_chop_threshold_) {
          return norm;
        } else {
          return 0.0;
        }
      });
    }

    // Contract C_df with orbitals
    C_mo("X, i, mu") = C_df_("X, mu, nu") * C("nu, i");
    C_mo.truncate();

    util::write_shape_tuple3D(
        C_mo, std::string("C_mo_shape") + std::to_string(iteration) + ".txt");
    

    // Get forced output shape
    TA::SparseShape<float> forced_shape;
    if (use_forced_shape_) {
      auto cadf_df_k_shape = [&](TA::Tensor<float> const &input) {
        auto &range = input.range();
        auto extent = range.extent_data();

        TA::Tensor<float> t(range, 0.0);

        std::unordered_set<int64_t> sig_mu;
        std::unordered_set<int64_t> sig_X;

        for (auto i = 0ul; i < extent[1]; ++i) {
          sig_mu.clear();
          sig_X.clear();

          // For every I we determine save a list of important X and mu
          for (auto X = 0ul; X < extent[0]; ++X) {
            for (auto mu = 0ul; mu < extent[2]; ++mu) {
              const auto val = input(X, i, mu);

              if (val >= force_threshold_) {
                sig_X.insert(X);
                sig_mu.insert(mu);
              }
            }
          }

          // Then for every X we mark all important mu. The output will have
          // significant mu X pairs that the input did not.
          for (auto const &X : sig_X) {
            for (auto const &mu : sig_mu) {
              t(X, i, mu) = std::numeric_limits<float>::max();
            }
          }
        }

        return t;
      };

      forced_shape = C_mo.shape().transform(cadf_df_k_shape);
    }

    // Construct F_df
    if (!use_forced_shape_) {
      F_df("X, i, mu") = E_mo("X, i, mu") - 0.5 * M_("X,Y") * C_mo("Y, i, mu");
    } else {
      array_type E_mo_forced;
      E_mo_forced("X,i,mu") =
          (E_("X, mu, nu") * C("nu,i")).set_shape(forced_shape);
      E_mo_forced.truncate();

      array_type R_df;
      R_df("X, i, mu") =
          (M_("X, Y") * C_mo("Y, i, mu")).set_shape(forced_shape);

      F_df("X, i, mu") = E_mo_forced("X, i, mu") - 0.5 * R_df("X, i, mu");
    }
    F_df.truncate();
    util::write_shape_tuple3D(
        F_df, std::string("F_df_shape") + std::to_string(iteration) + ".txt");

    L("mu, nu") = C_mo("X, i, mu") * F_df("X, i, nu");
    L.truncate();

    K("mu, nu") = L("mu, nu") + L("nu, mu");

    return K;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC_SCF_CADFBUILDERPRINTONLY_H
