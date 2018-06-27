#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CADF_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CADF_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/factory/ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/util.h"
#include "mpqc/math/external/tiledarray/array_info.h"
#include "mpqc/util/misc/time.h"

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/four_center_exact_k_diagonal_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/traditional_four_center_fock_builder.h"

#include <tiledarray.h>

#include <iostream>
#include <unordered_set>
#include <vector>

namespace mpqc {
namespace lcao {
namespace scf {

template<typename Tile, typename Policy, typename DirectArray>
class CADFFockBuilder : public FockBuilder<Tile, Policy> {
 public:
  using ArrayType = TA::DistArray<Tile, Policy>;

 private:
  DirectArray E_;  // <Κ |G| κ λ > three center two electron direct integrals
  ArrayType M_;    // <Κ |G| Λ > two center two electron coulomb integrals
  ArrayType Mchol_inv_;  // Chol(<Κ |G| Λ >)^-1
  ArrayType C_;          // CADF fitting coeffs

  // seCADF builder
  std::unique_ptr<ExactKDiagonalBuilder < Tile, Policy>> exactK_;

  float force_threshold_ = 0.0;
  double LMO_chop_threshold_ = 0.0;
  bool secadf_ = false;

  double energy_J_;
  double energy_K_;

  // Vectors to store timings
  std::vector<double> j_times_;         // Time J in a single step times
  std::vector<double> LMO_chop_times_;  // Time spent chopping orbitals
  std::vector<double> Z_times_;         // C_df_ * C times
  std::vector<double> shape_times_;     // Time to compute the forced shape
  std::vector<double> F_times_;         // E_mo_ + M * C_mo times
  std::vector<double> L_times_;         // (C_mo)^T * F_df times
  std::vector<double> K_times_;         // L^T + L times
  std::vector<double> se_times_;        // L^T + L times
  std::vector<double> Exch_times_;      // L^T + L times

  std::vector<std::array<double, 2>> LMO_sizes_;
  std::vector<std::array<double, 2>> LMO_chopped_sizes_;
  std::vector<std::array<double, 2>> Z_sizes_;
  std::vector<std::array<double, 2>> E_LMO_sizes_;
  std::vector<std::array<double, 2>> F_sizes_;

  TA::SparseShape<float> K_diagonal_;      // For SECADF
  TA::SparseShape<float> K_off_diagonal_;  // For SECADF

  inline TA::Tensor<float> diagonal_shape(TA::Tensor<float> const &in) {
    auto range = in.range();
    TA::Tensor<float> out(range, 0.0);
    auto lo = range.lobound_data();
    auto up = range.upbound_data();

    for (auto i = lo[0]; i < up[0]; ++i) {
      out(i, i) = in(i, i);  // Assume K is always symmetric
    }

    return out;
  }

  inline TA::Tensor<float> off_diagonal_shape(
      TA::Tensor<float> const &in) {
    auto range = in.range();
    TA::Tensor<float> out = in.clone();
    auto lo = range.lobound_data();
    auto up = range.upbound_data();

    for (auto i = lo[0]; i < up[0]; ++i) {
      out(i, i) = 0.0;
    }

    return out;
  }

 public:
  using BasisFactory = lcao::gaussian::Basis::Factory;

  template<typename Factory>
  CADFFockBuilder(Factory &ao_factory, double force_threshold,
                  double lmo_chop_threshold, bool do_secadf)
      : force_threshold_(force_threshold),
        LMO_chop_threshold_(lmo_chop_threshold),
        secadf_(do_secadf) {
    force_threshold_ = force_threshold;
    LMO_chop_threshold_ = lmo_chop_threshold;

    C_ = ao_factory.compute(L"( Κ | CADF|κ λ)");
    E_ = ao_factory.compute_direct(L"( Κ | G|κ λ)");
    M_ = ao_factory.compute(L"( Κ | G| Λ )");

    if (secadf_) {  // init four center builder
      auto &world = C_.world();
      auto &ao_factory_ref = ::mpqc::lcao::gaussian::to_ao_factory(ao_factory);
      auto screen = ao_factory_ref.screen();
      auto screen_threshold = ao_factory_ref.screen_threshold();
      auto obs = ao_factory.basis_registry()->retrieve(L"κ");

      exactK_ = std::make_unique<ExactKDiagonalBuilder < Tile, Policy>>
      (
          world, obs, obs, obs, screen, screen_threshold);
    }

    // Form L^{-1} for M
    auto M_eig = math::array_to_eigen(M_);
    using MatType = decltype(M_eig);
    MatType L_inv_eig = MatType(Eigen::LLT<MatType>(M_eig).matrixL()).inverse();

    auto trange1_M = M_.trange().data()[0];  // Assumes symmetric blocking
    Mchol_inv_ = math::eigen_to_array<Tile, Policy>(M_.world(), L_inv_eig,
                                                         trange1_M, trange1_M);
  }

  ~CADFFockBuilder() {}

  void register_fock(const ArrayType &fock,
                     FormulaRegistry<ArrayType> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)[df]"), fock);
  }

  ArrayType operator()(ArrayType const &D, ArrayType const &LMO,
                       double target_precision =
                       std::numeric_limits<double>::epsilon()) override {
    ArrayType G;
    G("m, n") =
        2 * compute_J(D)("m, n") - compute_K(LMO, D, target_precision)("m, n");
    return G;
  }

  void print_iter(std::string const &leader) override {
    ExEnv::out0() << indent << "J time: " << j_times_.back()
                  << ", Energy J: " << energy_J_ << "\n";
    ExEnv::out0() << indent << "K time: " << Exch_times_.back()
                  << ", Energy K: " << energy_K_ << "\n";
    if (!se_times_.empty()) {
      ExEnv::out0() << indent << indent
                    << "Semi-exact correction time: " << se_times_.back()
                    << "\n";
    }
  }

 private:
  ArrayType compute_J(ArrayType const &D) {
    auto &world = D.world();
    auto j0 = mpqc::fenced_now(world);
    ArrayType J;
    J("mu, nu") = E_("X, mu, nu") *
        (Mchol_inv_("Z, X") *
            (Mchol_inv_("Z, Y") * (E_("Y, sig, rho") * D("sig, rho"))));
    auto j1 = mpqc::fenced_now(world);
    j_times_.push_back(mpqc::duration_in_s(j0, j1));

    energy_J_ = D("i,j").dot(J("i,j"));

    return J;
  }

  ArrayType compute_K(ArrayType const &LMO_in, ArrayType const &D,
                      double target_precision) {
    using ::mpqc::detail::array_storage;

    auto &world = M_.world();
    auto Exch0 = mpqc::fenced_now(world);
    ArrayType L, K;  // Matrices
    ArrayType Z, F;  // Tensors

    // Deep copy LMO for chopping
    ArrayType LMO;
    LMO("mu, i") = LMO_in("mu, i");

    if (LMO_chop_threshold_ != 0.0) {
      auto chop0 = mpqc::fenced_now(world);
      TA::foreach_inplace(LMO, [&](TA::Tensor<double> &t) {
        const auto norm = t.norm();
        if (norm > LMO_chop_threshold_) {
          return norm;
        } else {
          return 0.0;
        }
      });
      auto chop1 = mpqc::fenced_now(world);
      LMO_chop_times_.push_back(mpqc::duration_in_s(chop0, chop1));

      auto chop_sizes = array_storage(LMO);
      LMO_chopped_sizes_.push_back(
          std::array<double, 2>{{chop_sizes[0], chop_sizes[1]}});
    }

    // Contract C with orbitals
    auto z0 = mpqc::fenced_now(world);
    Z("X, i, mu") = C_("X, mu, nu") * LMO("nu, i");
    Z.truncate();
    auto z1 = mpqc::fenced_now(world);
    Z_times_.push_back(mpqc::duration_in_s(z0, z1));

    auto Z_sizes = array_storage(Z);
    Z_sizes_.push_back(std::array<double, 2>{{Z_sizes[0], Z_sizes[1]}});

    // Get forced output shape
    TA::SparseShape<float> forced_shape;
    if (force_threshold_ > 0.0) {
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
      auto shape_time0 = mpqc::fenced_now(world);
      forced_shape = Z.shape().transform(cadf_df_k_shape);
      auto shape_time1 = mpqc::fenced_now(world);
      shape_times_.push_back(mpqc::duration_in_s(shape_time0, shape_time1));
    }

    // Construct F
    auto f0 = mpqc::fenced_now(world);
    if (force_threshold_ == 0.0) {
      F("X, i, mu") =
          E_("X, mu, nu") * LMO("nu, i") - 0.5 * M_("X,Y") * Z("Y, i, mu");
    } else {
      F("X, i, mu") = (E_("X, mu, nu") * LMO("nu, i")).set_shape(forced_shape);
      F.truncate();

      F("X, i, mu") +=
          (-0.5 * M_("X, Y") * Z("Y, i, mu")).set_shape(forced_shape);
    }
    F.truncate();
    auto f1 = mpqc::fenced_now(world);
    F_times_.push_back(mpqc::duration_in_s(f0, f1));

    auto F_sizes = array_storage(F);
    F_sizes_.push_back(std::array<double, 2>{{F_sizes[0], F_sizes[1]}});

    // Construct L
    auto l0 = mpqc::fenced_now(world);
    L("mu, nu") = Z("X, i, mu") * F("X, i, nu");
    L.truncate();
    auto l1 = mpqc::fenced_now(world);
    L_times_.push_back(mpqc::duration_in_s(l0, l1));

    auto k0 = mpqc::fenced_now(world);
    if (!secadf_) {
      K("mu, nu") = L("mu, nu") + L("nu, mu");
    } else {
      if (K_off_diagonal_.empty()) {
        K("mu, nu") = L("mu, nu") + L("nu, mu");
        world.gop.fence();
        auto off = [this](TA::Tensor<float> const &t) {
          return this->off_diagonal_shape(t);
        };
        K_off_diagonal_ = K.shape().transform(off);
      }
      K("mu, nu") = (L("mu, nu") + L("nu, mu")).set_shape(K_off_diagonal_);
    }
    K.truncate();

    // seCADF correction
    if (secadf_) {
      auto se0 = mpqc::fenced_now(world);
      auto K_se = exactK_->operator()(D, LMO, target_precision);

      if (K_diagonal_.empty()) {
        auto diag = [this](TA::Tensor<float> const &t) {
          return this->diagonal_shape(t);
        };
        K_diagonal_ = K_se.shape().transform(diag);
      }
      K_se("mu, nu") = K_se("mu, nu").set_shape(K_diagonal_);

      // Use - because exactK_ returns negative K
      K("mu, nu") = K("mu, nu") - K_se("mu, nu");
      auto se1 = mpqc::fenced_now(world);
      se_times_.push_back(mpqc::duration_in_s(se0, se1));
    }
    auto k1 = mpqc::fenced_now(world);
    K_times_.push_back(mpqc::duration_in_s(k0, k1));
    Exch_times_.push_back(mpqc::duration_in_s(Exch0, k1));

    energy_K_ = D("i,j").dot(K("i,j"));

    return K;
  }
};

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CADF_BUILDER_H_
