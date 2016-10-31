#pragma once
#ifndef MPQC_SCF_CLRCADFBUILDER_H
#define MPQC_SCF_CLRCADFBUILDER_H

#include "../../../../../common/namespaces.h"
#include <tiledarray.h>
#include "../../../../../utility/time.h"
#include "../../../../../utility/array_info.h"
#include "../../../../../utility/vector_functions.h"

#include "../../../../../tensor/decomposed_tensor.h"
#include "../../../../../tensor/mpqc_tile.h"
#include "../../../../../tensor/tensor_transforms.h"

#include "mpqc/math/external/eigen/eigen.h"
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

class ClrCADFFockBuilder : public FockBuilder {
 public:
  using TileType = TA::TensorD;
  using ArrayType = FockBuilder::array_type;
  using DArrayType =
      TA::DistArray<tensor::Tile<tensor::DecomposedTensor<double>>,
                    TA::SparsePolicy>;

 private:
  ArrayType E_;  // <Κ |G| κ λ > three center two electron coulomb integrals
  ArrayType Mchol_inv_;  // Chol(<Κ |G| Λ >)^-1
  DArrayType C_df_;      // CADF fitting coeffs
  DArrayType dE_;
  DArrayType M_;

  bool use_forced_shape_ = false;
  float force_threshold_ = 0.0;
  double lcao_chop_threshold_ = 0.0;
  double clr_threshold_ = 0.0;

  // Vectors to store timings
  std::vector<double> e_mo_times_;       // E_ * C times
  std::vector<double> r_df_times_;       // M * C_mo
  std::vector<double> j_times_;          // Time J in a single step times
  std::vector<double> lcao_chop_times_;  // Time spent chopping orbitals
  std::vector<double> c_mo_times_;       // C_df_ * C times
  std::vector<double> shape_times_;      // Time to compute the forced shape
  std::vector<double> f_df_times_;       // E_mo_ + M * C_mo times
  std::vector<double> l_times_;          // (C_mo)^T * F_df times
  std::vector<double> k_times_;          // L^T + L times

  std::array<double, 3> E_J_sizes_;
  std::array<double, 3> E_clr_sizes_;
  std::vector<std::array<double, 3>> e_mo_sizes_;
  std::vector<std::array<double, 3>> lcao_sizes_;
  std::vector<std::array<double, 3>> lcao_chopped_sizes_;
  std::vector<std::array<double, 3>> c_mo_sizes_;
  std::vector<std::array<double, 3>> f_df_sizes_;

  DArrayType make_three_center_integrals(basis::Basis const &obs,
                                         basis::Basis const &dfbs) {
    auto &world = E_.world();

    auto basis_array = std::vector<basis::Basis>{{dfbs, obs, obs}};

    auto eng_pool = integrals::make_engine_pool(
        libint2::Operator::coulomb, utility::make_array_of_refs(dfbs, obs, obs),
        libint2::BraKet::xs_xx);

    auto p_screen =
        std::make_shared<integrals::Screener>(integrals::Screener{});

    auto my_class = tensor::TaToDecompTensor(clr_threshold_);
    using clr_type = decltype(my_class(std::declval<TA::TensorD>()));
    auto func = std::function<clr_type(TA::TensorD &&)>(my_class);

    return integrals::sparse_integrals(world, eng_pool, basis_array,
                                           p_screen, func);
  }

 public:
  ClrCADFFockBuilder(
      molecule::Molecule const &clustered_mol,
      molecule::Molecule const &df_clustered_mol,
      basis::BasisSet const &obs_set, basis::BasisSet const &dfbs_set,
      integrals::AtomicIntegral<TileType, TA::SparsePolicy> &ao_ints,
      bool use_forced_shape, double force_threshold, double lcao_chop_threshold,
      double clr_thresh)
      : ClrCADFFockBuilder(clustered_mol, df_clustered_mol, obs_set, dfbs_set,
                           ao_ints, clr_thresh) {
    use_forced_shape_ = use_forced_shape;
    force_threshold_ = force_threshold;
    lcao_chop_threshold_ = lcao_chop_threshold;

    basis::Basis obs = ao_ints.orbital_basis_registry().retrieve(L"κ");
    basis::Basis dfbs = ao_ints.orbital_basis_registry().retrieve(L"Κ");

    dE_ = make_three_center_integrals(obs, dfbs);

    E_clr_sizes_ = utility::array_storage(dE_);
    if (dE_.world().rank() == 0) {
      std::cout << "E with clr storage:\n"
                << "\tDense  " << E_clr_sizes_[0] << "\n"
                << "\tSparse " << E_clr_sizes_[1] << "\n"
                << "\tCLR    " << E_clr_sizes_[2] << std::endl;
    }
  }

  ClrCADFFockBuilder(
      molecule::Molecule const &clustered_mol,
      molecule::Molecule const &df_clustered_mol,
      basis::BasisSet const &obs_set, basis::BasisSet const &dfbs_set,
      integrals::AtomicIntegral<TileType, TA::SparsePolicy> &ao_ints,
      double clr_threshold)
      : FockBuilder(), clr_threshold_(clr_threshold) {
    // Grab needed ao integrals
    E_ = ao_ints.compute(L"( Κ | G|κ λ)");
    E_J_sizes_ = utility::array_storage(E_);
    if (E_.world().rank() == 0) {
      std::cout << "E for J storage:\n"
                << "\tDense  " << E_J_sizes_[0] << "\n"
                << "\tSparse " << E_J_sizes_[1] << std::endl;
    }
    ArrayType M = ao_ints.compute(L"( Κ | G| Λ )");

    bool compress_M = false;
    M_ = TA::to_new_tile_type(
        M, tensor::TaToDecompTensor(clr_threshold_, compress_M));
    M_.world().gop.fence();

    auto m_store = utility::array_storage(M_);
    if (E_.world().rank() == 0) {
      std::cout << "M storage:\n"
                << "\tDense  " << m_store[0] << "\n"
                << "\tSparse " << m_store[1] << "\n"
                << "\tCLR    " << m_store[2] << std::endl;
    }

    // Form L^{-1} for M
    auto M_eig = array_ops::array_to_eigen(M);
    using MatType = decltype(M_eig);
    MatType L_inv_eig = MatType(Eig::LLT<MatType>(M_eig).matrixL()).inverse();

    auto trange1_M = M.trange().data()[0];  // Assumes symmetric blocking
    Mchol_inv_ = array_ops::eigen_to_array<TA::TensorD>(
        M.world(), L_inv_eig, trange1_M, trange1_M);

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

    ArrayType C_df =
        scf::reblock_from_atoms(C_df_temp, obs_atom_to_cluster_map,
                                dfbs_atom_to_cluster_map, by_cluster_trange);
    C_df_ =
        TA::to_new_tile_type(C_df, tensor::TaToDecompTensor(clr_threshold_));
    E_.world().gop.fence();
    auto c_df_store = utility::array_storage(C_df_);
    if (E_.world().rank() == 0) {
      std::cout << "C_df storage:\n"
                << "\tDense  " << c_df_store[0] << "\n"
                << "\tSparse " << c_df_store[1] << "\n"
                << "\tCLR    " << c_df_store[2] << std::endl;
    }
  }

  ~ClrCADFFockBuilder() = default;

  void register_fock(const TA::TSpArrayD &fock,
                     FormulaRegistry<TA::TSpArrayD> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)[df]"), fock);
  }

  ArrayType operator()(ArrayType const &D, ArrayType const &C) override {
    ArrayType G;
    G("m, n") = 2 * compute_J(D)("m, n") - compute_K(C)("m, n");
    return G;
  }

  void print_iter(std::string const &leader) override {
    if (E_.world().rank() == 0) {
      auto jt = j_times_.back();
      auto ct = c_mo_times_.back();
      auto ft = f_df_times_.back();
      auto et = e_mo_times_.back();
      auto rt = r_df_times_.back();
      auto lt = l_times_.back();
      auto kt = k_times_.back();  // L^T + L

      auto shape_time = 0.0;
      if (!shape_times_.empty()) {
        shape_time = shape_times_.back();
      }

      auto cut_time = 0.0;
      if (!lcao_chop_times_.empty()) {
        cut_time = lcao_chop_times_.back();
      }

      auto ktotal = 0.0;
      ktotal = ct + ft + lt + kt + shape_time + cut_time;

      std::cout << leader << "J time   : " << jt << "\n";
      if (!shape_times_.empty()) {
        std::cout << leader << "Shape Time: " << shape_time << "\n";
      }
      std::cout << leader << "C_mo time: " << ct << "\n";
      std::cout << leader << "F_df time: " << ft << "\n";
      std::cout << leader << "\tE_mo time: " << et << "\n";
      std::cout << leader << "\tR_df time: " << rt << "\n";
      std::cout << leader << "L time   : " << lt << "\n";
      std::cout << leader << "K time   : " << kt << "\n";
      std::cout << leader << "Exchange time   : " << ktotal << "\n";
      std::cout << leader << "Storages(Dense, Sparse, CLR):\n";
      std::cout << leader << "\tE_mo(" << e_mo_sizes_.back()[0] << ", "
                << e_mo_sizes_.back()[1] << ", " << e_mo_sizes_.back()[2]
                << ")\n";
      std::cout << leader << "\tlcao(" << lcao_sizes_.back()[0] << ", "
                << lcao_sizes_.back()[1] << ")\n";
      if (!lcao_chopped_sizes_.empty()) {
        std::cout << leader << "\tlcao chopped("
                  << lcao_chopped_sizes_.back()[0] << ", "
                  << lcao_chopped_sizes_.back()[1] << ")\n";
      }
      std::cout << leader << "\tC_mo(" << c_mo_sizes_.back()[0] << ", "
                << c_mo_sizes_.back()[1] << ", " << c_mo_sizes_.back()[2]
                << ")\n";
      std::cout << leader << "\tF_df(" << f_df_sizes_.back()[0] << ", "
                << f_df_sizes_.back()[1] << ", " << f_df_sizes_.back()[2] << ")"
                << std::endl;
    }
  }

  rapidjson::Value results(rapidjson::Document &d) override {
    rapidjson::Value fock_builder(rapidjson::kObjectType);
    fock_builder.AddMember("Type", "ClrCADFFockBuilder", d.GetAllocator());
    fock_builder.AddMember("clr threshold", clr_threshold_, d.GetAllocator());

    auto e_mo_build = utility::vec_avg(e_mo_times_);
    auto r_df_build = utility::vec_avg(r_df_times_);
    auto j_build = utility::vec_avg(j_times_);
    auto c_mo_build = utility::vec_avg(c_mo_times_);
    auto f_df_build = utility::vec_avg(f_df_times_);
    auto l_build = utility::vec_avg(l_times_);
    auto k_build = utility::vec_avg(k_times_);

    auto cut_time = 0.0;
    if (lcao_chop_threshold_ != 0.0) {
      auto cut_time = utility::vec_avg(lcao_chop_times_);
      fock_builder.AddMember("Avg LCAO Cut Time", cut_time, d.GetAllocator());
    }

    auto shape_time = 0.0;
    if (!shape_times_.empty()) {
      shape_time = utility::vec_avg(shape_times_);
    }

    auto ktotal =
        c_mo_build + f_df_build + l_build + k_build + shape_time + cut_time;

    fock_builder.AddMember("Forced Shape", use_forced_shape_, d.GetAllocator());
    if (use_forced_shape_) {
      fock_builder.AddMember("Forced Shape Threshold", force_threshold_,
                             d.GetAllocator());
      fock_builder.AddMember("Avg Shape Time", shape_time, d.GetAllocator());
    }
    if (lcao_chop_threshold_ != 0.0) {
      fock_builder.AddMember("LCAO Cut Threshold", lcao_chop_threshold_,
                             d.GetAllocator());
      fock_builder.AddMember("Avg LCAO Cut Time", cut_time, d.GetAllocator());
    }
    fock_builder.AddMember("Avg J Time", j_build, d.GetAllocator());
    fock_builder.AddMember("Avg C_mo Time", c_mo_build, d.GetAllocator());
    fock_builder.AddMember("Avg E_mo Time", e_mo_build, d.GetAllocator());
    fock_builder.AddMember("Avg R_df Time", r_df_build, d.GetAllocator());
    fock_builder.AddMember("Avg F_df Time", f_df_build, d.GetAllocator());
    fock_builder.AddMember("Avg L Time", l_build, d.GetAllocator());
    fock_builder.AddMember("Avg K Time", k_build, d.GetAllocator());
    fock_builder.AddMember("Avg Total Exchange Time", ktotal, d.GetAllocator());

    auto average_storage = [](std::vector<std::array<double, 3>> const &store) {
      auto dense = 0.0;
      auto sparse = 0.0;
      auto clr = 0.0;

      for (auto const &s : store) {
        dense += s[0];
        sparse += s[1];
        clr += s[2];
      }
      dense /= double(store.size());
      sparse /= double(store.size());
      clr /= double(store.size());

      return std::array<double, 3>{{dense, sparse, clr}};
    };

    fock_builder.AddMember("E Dense", E_clr_sizes_[0], d.GetAllocator());
    fock_builder.AddMember("E Sparse", E_clr_sizes_[1], d.GetAllocator());
    fock_builder.AddMember("E CLR", E_clr_sizes_[2], d.GetAllocator());

    auto store = average_storage(e_mo_sizes_);
    fock_builder.AddMember("Avg E_mo Dense", store[0], d.GetAllocator());
    fock_builder.AddMember("Avg E_mo Sparse", store[1], d.GetAllocator());

    store = average_storage(lcao_sizes_);
    fock_builder.AddMember("Avg C Dense", store[0], d.GetAllocator());
    fock_builder.AddMember("Avg C Sparse", store[1], d.GetAllocator());

    if (!lcao_chopped_sizes_.empty()) {
      store = average_storage(lcao_chopped_sizes_);
      fock_builder.AddMember("Avg C(chopped) Dense", store[0],
                             d.GetAllocator());
      fock_builder.AddMember("Avg C(chopped) Sparse", store[1],
                             d.GetAllocator());
    }

    store = average_storage(c_mo_sizes_);
    fock_builder.AddMember("Avg C_mo Dense", store[0], d.GetAllocator());
    fock_builder.AddMember("Avg C_mo Sparse", store[1], d.GetAllocator());

    store = average_storage(f_df_sizes_);
    fock_builder.AddMember("Avg F_df Dense", store[0], d.GetAllocator());
    fock_builder.AddMember("Avg F_df Sparse", store[1], d.GetAllocator());

    return fock_builder;
  }

 private:
  ArrayType compute_J(ArrayType const &D) {
    auto &world = D.world();
    auto j0 = mpqc_time::fenced_now(world);
    ArrayType J;
    J("mu, nu") =
        E_("X, mu, nu") *
        (Mchol_inv_("Z, X") *
         (Mchol_inv_("Z, Y") * (E_("Y, rho, sigma") * D("rho, sigma"))));
    auto j1 = mpqc_time::fenced_now(world);
    j_times_.push_back(mpqc_time::duration_in_s(j0, j1));

    return J;
  }

  array_type compute_K(ArrayType const &C_in) {
    auto &world = M_.world();
    ArrayType K;            // Matrices
    DArrayType dC, dL;      // Matrices
    DArrayType C_mo, F_df;  // Tensors

    // Deep copy C for chopping
    ArrayType C;
    C("mu, i") = C_in("mu, i");
    // Capture C sizes
    lcao_sizes_.push_back(utility::array_storage(C));

    if (lcao_chop_threshold_ != 0.0) {
      auto chop0 = mpqc_time::fenced_now(world);
      TA::foreach_inplace(C, [&](TA::Tensor<double> &t) {
        const auto norm = t.norm();
        if (norm > lcao_chop_threshold_) {
          return norm;
        } else {
          return 0.0;
        }
      });
      auto chop1 = mpqc_time::fenced_now(world);
      lcao_chop_times_.push_back(mpqc_time::duration_in_s(chop0, chop1));

      lcao_chopped_sizes_.push_back(utility::array_storage(C));
    }

    bool compress_C = false;
    dC = TA::to_new_tile_type(
        C, tensor::TaToDecompTensor(clr_threshold_, compress_C));

    // Contract C_df with orbitals
    auto c_mo0 = mpqc_time::fenced_now(world);
    C_mo("X, i, mu") = C_df_("X, mu, nu") * dC("nu, i");
    C_mo.truncate();
    auto c_mo1 = mpqc_time::fenced_now(world);
    c_mo_times_.push_back(mpqc_time::duration_in_s(c_mo0, c_mo1));
    c_mo_sizes_.push_back(utility::array_storage(C_mo));

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
      forced_shape = C_mo.shape().transform(cadf_df_k_shape);
      auto shape_time1 = mpqc_time::fenced_now(world);
      shape_times_.push_back(
          mpqc_time::duration_in_s(shape_time0, shape_time1));
    }

    // Construct F_df
    auto f_df0 = mpqc_time::fenced_now(world);
    if (!use_forced_shape_) {
      auto emo0 = mpqc_time::fenced_now(world);
      DArrayType E_mo;
      E_mo("X,i,mu") = dE_("X, mu, nu") * dC("nu,i");
      E_mo.truncate();
      auto emo1 = mpqc_time::fenced_now(world);
      e_mo_times_.push_back(mpqc_time::duration_in_s(emo0, emo1));
      e_mo_sizes_.push_back(utility::array_storage(E_mo));

      auto rdf0 = mpqc_time::fenced_now(world);
      DArrayType R_df;
      R_df("X, i, mu") = M_("X, Y") * C_mo("Y, i, mu");
      auto rdf1 = mpqc_time::fenced_now(world);
      r_df_times_.push_back(mpqc_time::duration_in_s(rdf0, rdf1));
      F_df("X, i, mu") = E_mo("X, i, mu") - 0.5 * R_df("X, i, mu");
    } else {
      auto emo0 = mpqc_time::fenced_now(world);
      DArrayType E_mo;
      E_mo("X,i,mu") = (dE_("X, mu, nu") * dC("nu,i")).set_shape(forced_shape);
      E_mo.truncate();
      auto emo1 = mpqc_time::fenced_now(world);
      e_mo_times_.push_back(mpqc_time::duration_in_s(emo0, emo1));
      e_mo_sizes_.push_back(utility::array_storage(E_mo));

      auto rdf0 = mpqc_time::fenced_now(world);
      DArrayType R_df;
      R_df("X, i, mu") =
          (M_("X, Y") * C_mo("Y, i, mu")).set_shape(forced_shape);
      auto rdf1 = mpqc_time::fenced_now(world);
      r_df_times_.push_back(mpqc_time::duration_in_s(rdf0, rdf1));

      F_df("X, i, mu") = E_mo("X, i, mu") - 0.5 * R_df("X, i, mu");
    }
    F_df.truncate();
    auto f_df1 = mpqc_time::fenced_now(world);
    f_df_times_.push_back(mpqc_time::duration_in_s(f_df0, f_df1));
    f_df_sizes_.push_back(utility::array_storage(F_df));

    // Construct L
    auto l0 = mpqc_time::fenced_now(world);
    dL("mu, nu") = C_mo("X, i, mu") * F_df("X, i, nu");
    dL.truncate();
    auto L = TA::to_new_tile_type(dL, tensor::DecompToTaTensor{});
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

#endif  // MPQC_SCF_CLRCADFBUILDER_H
