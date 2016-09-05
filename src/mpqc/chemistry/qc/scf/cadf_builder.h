#pragma once
#ifndef MPQC_SCF_CADFBUILDER_H
#define MPQC_SCF_CADFBUILDER_H

#include "../../../../../common/namespaces.h"
#include "../../../../../include/tiledarray.h"
#include "../../../../../utility/time.h"
#include "../../../../../utility/array_info.h"
#include "../../../../../utility/vector_functions.h"

#include "../../../../../tensor/decomposed_tensor.h"
#include "../../../../../tensor/mpqc_tile.h"
#include "../../../../../tensor/tensor_transforms.h"

#include "../../../../../ta_routines/array_to_eigen.h"
#include "../../../../../ta_routines/minimize_storage.h"

#include <mpqc/chemistry/qc/scf/builder.h>
#include <mpqc/chemistry/qc/integrals/make_engine.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral.h>
#include <mpqc/chemistry/qc/scf/cadf_fitting_coeffs.h>
#include <mpqc/chemistry/qc/scf/cadf_helper_functions.h>

#include <vector>
#include <iostream>
#include <unordered_set>

namespace mpqc {
namespace scf {

class CADFFockBuilder : public FockBuilder {
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

  // Vectors to store timings
  std::vector<double> e_mo_times_;   // E_ * C times
  std::vector<double> j_times_;      // Time J in a single step times
  std::vector<double> c_mo_times_;   // C_df_ * C times
  std::vector<double> shape_times_;  // Time to compute the forced shape
  std::vector<double> f_df_times_;   // E_mo_ + M * C_mo times
  std::vector<double> l_times_;      // (C_mo)^T * F_df times
  std::vector<double> k_times_;      // L^T + L times

 public:
  CADFFockBuilder(
      molecule::Molecule const &clustered_mol,
      molecule::Molecule const &df_clustered_mol,
      basis::BasisSet const &obs_set, basis::BasisSet const &dfbs_set,
      integrals::AtomicIntegral<TileType, TA::SparsePolicy> &ao_ints,
      bool use_forced_shape, double force_threshold)
      : CADFFockBuilder(clustered_mol, df_clustered_mol, obs_set, dfbs_set,
                        ao_ints) {
    use_forced_shape_ = use_forced_shape;
    force_threshold_ = force_threshold;
  }

  CADFFockBuilder(
      molecule::Molecule const &clustered_mol,
      molecule::Molecule const &df_clustered_mol,
      basis::BasisSet const &obs_set, basis::BasisSet const &dfbs_set,
      integrals::AtomicIntegral<TileType, TA::SparsePolicy> &ao_ints)
      : FockBuilder() {
    // Grab needed ao integrals
    E_ = ao_ints.compute(L"( Κ | G|κ λ)");
    M_ = ao_ints.compute(L"( Κ | G| Λ )");

    // Form L^{-1} for M
    auto M_eig = array_ops::array_to_eigen(M_);
    using MatType = decltype(M_eig);
    MatType L_inv_eig = MatType(Eig::LLT<MatType>(M_eig).matrixL()).inverse();

    auto trange1_M = M_.trange().data()[0];  // Assumes symmetric blocking
    Mchol_inv_ = array_ops::eigen_to_array<TA::TensorD>(
        M_.get_world(), L_inv_eig, trange1_M, trange1_M);

    std::unordered_map<std::size_t, std::size_t> obs_atom_to_cluster_map;
    std::unordered_map<std::size_t, std::size_t> dfbs_atom_to_cluster_map;

    basis::Basis obs = ao_ints.orbital_basis_registry()->retrieve(L"κ");
    basis::Basis dfbs = ao_ints.orbital_basis_registry()->retrieve(L"Κ");

    auto ref_array = utility::make_array_of_refs(dfbs, dfbs);
    auto eng_pool = integrals::make_engine_pool(
        libint2::Operator::coulomb, utility::make_array_of_refs(dfbs, dfbs),
        libint2::BraKet::xs_xs);

    ArrayType C_df_temp = scf::compute_atomic_fitting_coeffs(
        M_.get_world(), clustered_mol, df_clustered_mol, obs_set, dfbs_set,
        eng_pool, obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);

    auto by_cluster_trange =
        integrals::detail::create_trange(utility::make_array(dfbs, obs, obs));

    C_df_ =
        scf::reblock_from_atoms(C_df_temp, obs_atom_to_cluster_map,
                                dfbs_atom_to_cluster_map, by_cluster_trange);
  }

  ~CADFFockBuilder() = default;

  ArrayType operator()(ArrayType const &D, ArrayType const &C) override {
    auto &world = D.get_world();

    auto e_mo0 = mpqc_time::fenced_now(world);
    ArrayType E_mo;  // Temp array shared by J and K
    E_mo("X, i, mu") = E_("X, mu, nu") * C("nu, i");
    E_mo.truncate();
    auto e_mo1 = mpqc_time::fenced_now(world);
    e_mo_times_.push_back(mpqc_time::duration_in_s(e_mo0, e_mo1));

    ArrayType G;
    G("m, n") = 2 * compute_J(C, E_mo)("m, n") - compute_K(C, E_mo)("m, n");
    return G;
  }

  void print_iter(std::string const &leader) override {
    if (E_.get_world().rank() == 0) {
      auto et = e_mo_times_.back();
      auto jt = j_times_.back();
      auto ct = c_mo_times_.back();
      auto ft = f_df_times_.back();
      auto lt = l_times_.back();
      auto kt = k_times_.back();  // L^T + L

      auto shape_time = 0.0;
      if (!shape_times_.empty()) {
        shape_time = shape_times_.back();
      }
      auto ktotal = et + ct + ft + lt + kt + shape_time;
      std::cout << leader << "E_mo time: " << et << "\n";
      std::cout << leader << "J time   : " << jt << "\n";
      if (!shape_times_.empty()) {
        std::cout << leader << "Shape Time: " << shape_time << "\n";
      }
      std::cout << leader << "F_df time: " << ft << "\n";
      std::cout << leader << "L time   : " << lt << "\n";
      std::cout << leader << "K time   : " << kt << "\n";
      std::cout << leader << "Exchange time   : " << ktotal << "\n";
    }
  }

  rapidjson::Value results(rapidjson::Document &d) override {
    rapidjson::Value fock_builder(rapidjson::kObjectType);
    fock_builder.AddMember("Type", "CADFFockBuilder", d.GetAllocator());

    auto e_mo_build = utility::vec_avg(e_mo_times_);
    auto j_build = utility::vec_avg(j_times_);
    auto c_mo_build = utility::vec_avg(c_mo_times_);
    auto f_df_build = utility::vec_avg(f_df_times_);
    auto l_build = utility::vec_avg(l_times_);
    auto k_build = utility::vec_avg(k_times_);
    auto shape_time = 0.0;
    if (!shape_times_.empty()) {
      shape_time = utility::vec_avg(shape_times_);
    }
    auto ktotal =
        e_mo_build + c_mo_build + f_df_build + l_build + k_build + shape_time;

    fock_builder.AddMember("Forced Sahpe", use_forced_shape_, d.GetAllocator());
    if (use_forced_shape_) {
      fock_builder.AddMember("Forced Shape Threshold", force_threshold_,
                             d.GetAllocator());
      fock_builder.AddMember("Avg Shape Time", shape_time, d.GetAllocator());
    }
    fock_builder.AddMember("Avg Emo Time", e_mo_build, d.GetAllocator());
    fock_builder.AddMember("Avg J Time", j_build, d.GetAllocator());
    fock_builder.AddMember("Avg Cmo Time", c_mo_build, d.GetAllocator());
    fock_builder.AddMember("Avg Fdf Time", f_df_build, d.GetAllocator());
    fock_builder.AddMember("Avg L Time", l_build, d.GetAllocator());
    fock_builder.AddMember("Avg K Time", k_build, d.GetAllocator());
    fock_builder.AddMember("Avg Total Exchange Time", ktotal, d.GetAllocator());

    return fock_builder;
  }

 private:
  ArrayType compute_J(ArrayType const &C, ArrayType const &E_mo) {
    auto &world = C.get_world();
    auto j0 = mpqc_time::fenced_now(world);
    ArrayType J;
    J("mu, nu") = E_("X, mu, nu") *
                  (Mchol_inv_("Z, X") *
                   (Mchol_inv_("Z, Y") * (E_mo("Y, i, rho") * C("rho, i"))));
    auto j1 = mpqc_time::fenced_now(world);
    j_times_.push_back(mpqc_time::duration_in_s(j0, j1));

    return J;
  }

  array_type compute_K(ArrayType const &C, ArrayType const &E_mo) {
    auto &world = M_.get_world();
    ArrayType L, K;        // Matrices
    ArrayType C_mo, F_df;  // Tensors

    // Contract C_df with orbitals
    auto c_mo0 = mpqc_time::fenced_now(world);
    C_mo("X, i, mu") = C_df_("X, mu, nu") * C("nu, i");
    C_mo.truncate();
    auto c_mo1 = mpqc_time::fenced_now(world);
    c_mo_times_.push_back(mpqc_time::duration_in_s(c_mo0, c_mo1));

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

      auto shape_time0 = mpqc_time::fenced_now(world);
      forced_shape = C_mo.get_shape().transform(cadf_df_k_shape);
      auto shape_time1 = mpqc_time::fenced_now(world);
      shape_times_.push_back(
          mpqc_time::duration_in_s(shape_time0, shape_time1));
    }

    // Construct F_df
    auto f_df0 = mpqc_time::fenced_now(world);
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
    auto f_df1 = mpqc_time::fenced_now(world);
    f_df_times_.push_back(mpqc_time::duration_in_s(f_df0, f_df1));

    // Construct L
    auto l0 = mpqc_time::fenced_now(world);
    L("mu, nu") = C_mo("X, i, mu") * F_df("X, i, nu");
    L.truncate();
    auto l1 = mpqc_time::fenced_now(world);
    l_times_.push_back(mpqc_time::duration_in_s(l0, l1));

    auto k0 = mpqc_time::fenced_now(world);
    l_times_.push_back(mpqc_time::duration_in_s(l0, l1));
    K("mu, nu") = L("mu, nu") + L("nu, mu");
    K.truncate();
    auto k1 = mpqc_time::fenced_now(world);
    k_times_.push_back(mpqc_time::duration_in_s(k0, k1));

    return K;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC_SCF_CADFBUILDER_H
