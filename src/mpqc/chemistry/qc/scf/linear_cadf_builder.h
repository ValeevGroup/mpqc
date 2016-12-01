
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINEAR_CADF_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINEAR_CADF_BUILDER_H_

#include "mpqc/chemistry/qc/scf/util.h"
#include "mpqc/math/external/tiledarray/array_info.h"
#include "mpqc/util/misc/time.h"
#include <tiledarray.h>

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/tensor_transforms.h"
#include "mpqc/math/tensor/clr/tile.h"

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/sqrt_inv.h"
#include "mpqc/math/tensor/clr/minimize_storage.h"

#include "mpqc/chemistry/qc/scf/builder.h"
#include "mpqc/chemistry/qc/scf/ta_shape_tracker.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_set>
#include <vector>

namespace mpqc {
namespace scf {

template <typename Integral>
class ONCADFFockBuilder : public FockBuilder {
 public:
  using dtile_type = tensor::Tile<tensor::DecomposedTensor<double>>;
  using darray_type = TA::DistArray<dtile_type, TA::SparsePolicy>;

 private:
  Integral E_;              // Direct three center integrals
  darray_type M_;           // Fitting Metric for K
  darray_type M_sqrt_inv_;  // Sqrt inverse of M
  darray_type C_df_;        // CADF fitting coeffs

  double clr_thresh_;

  std::vector<double> j_times_;
  std::vector<double> c_mo_times_;
  std::vector<double> e_df_times_;
  std::vector<double> r_df_times_;

  std::vector<double> dl_times_;
  std::vector<double> dl_to_k_times_;

  std::vector<std::array<double, 3>> u_mo_storages_;
  std::vector<std::array<double, 3>> c_mo_storages_;
  std::vector<std::array<double, 3>> e_df_storages_;
  std::vector<std::array<double, 3>> f_df_storages_;

  bool force_shape_;
  double force_cut_thresh_;
  double mo_cut_thresh_;

  //  ShapeTracker shape_tracker_;
  //  std::vector<ShapeTracker> shape_tracker_iters_;

 public:
  ONCADFFockBuilder(darray_type const &Mj, darray_type const &Mk,
                    Integral const &eri3, darray_type const &C_df,
                    double clr_thresh, double force_cut_thresh,
                    double mo_cut_thresh, bool force_shape)
      : FockBuilder(),
        E_(eri3),
        M_(Mk),
        C_df_(C_df),
        clr_thresh_(clr_thresh),
        force_cut_thresh_(force_cut_thresh),
        mo_cut_thresh_(mo_cut_thresh),
        force_shape_(force_shape) {
    auto &world = C_df_.world();
    bool compress_M = false;

    auto l0 = mpqc::fenced_now(world);
    M_sqrt_inv_ = array_ops::inverse_sqrt(Mj);
    auto l1 = mpqc::fenced_now(world);
    auto l_inv_time_ = mpqc::duration_in_s(l0, l1);

    if (world.rank() == 0) {
      std::cout << "Sqrt Inv of Metric time: " << l_inv_time_ << std::endl;
    }
  }

  ~ONCADFFockBuilder() = default;

  array_type operator()(array_type const &D, array_type const &C) override {
    array_type G;
    G("mu, nu") = 2 * compute_J(D)("mu, nu") - compute_K(D, C)("mu, nu");
    return G;
  }

  void print_iter(std::string const &leader) override {
    auto &world = C_df_.world();
    if (world.rank() == 0) {
      std::cout
          << leader << "CADF Builder:\n"
          << leader << "\tStorages:\n"
          << leader << "\t\tFully Dense: " << c_mo_storages_.back()[0]
          << " GB\n"
          << leader << "\t\tU mo dense storage: " << u_mo_storages_.back()[0]
          << " GB\n"
          << leader << "\t\tU mo sparse storage: " << u_mo_storages_.back()[1]
          << " GB\n"
          << leader << "\t\tC mo sparse storage: " << c_mo_storages_.back()[1]
          << " GB\n"
          << leader << "\t\tC mo clr storage: " << c_mo_storages_.back()[2]
          << " GB\n"
          << leader << "\t\tE df storage: " << e_df_storages_.back()[1]
          << " GB\n"
          << leader << "\t\tE df clr storage: " << e_df_storages_.back()[2]
          << " GB\n"
          << leader << "\t\tF df storage: " << f_df_storages_.back()[1]
          << " GB\n"
          << leader << "\t\tF df clr storage: " << f_df_storages_.back()[2]
          << " GB\n"
          << leader << "\tJ time: " << j_times_.back() << "\n"
          << leader << "\tC mo time: " << c_mo_times_.back() << "\n"
          << leader << "\tE df time: " << e_df_times_.back() << "\n"
          << leader << "\tR df time: " << r_df_times_.back() << "\n"
          << leader << "\tdL time: " << dl_times_.back() << "\n"
          << leader << "\tdL to K time: " << dl_to_k_times_.back() << "\n";
      auto total_k_time = c_mo_times_.back() + e_df_times_.back() +
                          r_df_times_.back() + dl_times_.back() +
                          dl_to_k_times_.back();
      std::cout << leader << "\ttotal K time: " << total_k_time << "\n";
    }
  }

  rapidjson::Value results(rapidjson::Document &d) override {
    rapidjson::Value fock_builder(rapidjson::kObjectType);
    fock_builder.AddMember("Type", "ONCADFFockBuilder", d.GetAllocator());
    fock_builder.AddMember("Force shape", force_shape_, d.GetAllocator());

    auto j_build = utility::vec_avg(j_times_);
    auto c_mo_build = utility::vec_avg(c_mo_times_);
    auto e_df_build = utility::vec_avg(e_df_times_);
    auto r_df_build = utility::vec_avg(r_df_times_);
    auto dl_build = utility::vec_avg(dl_times_);
    auto k_build = utility::vec_avg(dl_to_k_times_);

    fock_builder.AddMember("Avg J Build Time", j_build, d.GetAllocator());
    fock_builder.AddMember("Avg C_mo Build Time", c_mo_build, d.GetAllocator());
    fock_builder.AddMember("Avg E_df Build Time", e_df_build, d.GetAllocator());
    fock_builder.AddMember("Avg R_df Build Time", r_df_build, d.GetAllocator());
    fock_builder.AddMember("Avg L Build Time", dl_build, d.GetAllocator());
    fock_builder.AddMember("Avg K from L Build Time", k_build,
                           d.GetAllocator());
    fock_builder.AddMember(
        "Avg Total K Build Time",
        k_build + dl_build + e_df_build + r_df_build + c_mo_build,
        d.GetAllocator());

    auto c_mo_storage = utility::vec_avg(c_mo_storages_);
    fock_builder.AddMember("C_mo Dense Storage", c_mo_storage[0],
                           d.GetAllocator());
    fock_builder.AddMember("C_mo Avg Sparse Storage", c_mo_storage[1],
                           d.GetAllocator());
    fock_builder.AddMember("C_mo Avg CLR Storage", c_mo_storage[2],
                           d.GetAllocator());
    auto e_df_storage = utility::vec_avg(e_df_storages_);
    fock_builder.AddMember("E_df Dense Storage", e_df_storage[0],
                           d.GetAllocator());
    fock_builder.AddMember("E_df Avg Sparse Storage", e_df_storage[1],
                           d.GetAllocator());
    fock_builder.AddMember("E_df Avg CLR Storage", e_df_storage[2],
                           d.GetAllocator());
    auto f_df_storage = utility::vec_avg(f_df_storages_);
    fock_builder.AddMember("F_df Dense Storage", f_df_storage[0],
                           d.GetAllocator());
    fock_builder.AddMember("F_df Avg Sparse Storage", f_df_storage[1],
                           d.GetAllocator());
    fock_builder.AddMember("F_df Avg CLR Storage", f_df_storage[2],
                           d.GetAllocator());

    return fock_builder;
  }

 private:
  array_type compute_J(array_type const &D) {
    auto &world = C_df_.world();

    constexpr bool compress_D = false;
    auto dD = TA::to_new_tile_type(
        D, tensor::TaToDecompTensor(clr_thresh_, compress_D));

    darray_type dJ;
    auto j0 = mpqc::fenced_now(world);

    dJ("mu, nu") =
        E_("X,mu,nu") * (M_sqrt_inv_("X,Z") *
                         (M_sqrt_inv_("Z, Y") * (E_("Y,r,s") * dD("r,s"))));

    auto j1 = mpqc::fenced_now(world);
    j_times_.push_back(mpqc::duration_in_s(j0, j1));

    auto J = TA::to_new_tile_type(dJ, tensor::DecompToTaTensor{});

    return J;
  }

  struct Rdf_shape {
    double th_;
    Rdf_shape(double th) : th_(th) {
      std::cout << "Thresh = " << th_ << std::endl;
    }
    TA::Tensor<float> operator()(TA::Tensor<float> const &norms) {
      auto &range = norms.range();
      auto lo = range.lobound_data();
      auto up = range.upbound_data();

      // int32_t for size at the moment, may need to make bigger in the future.
      std::vector<std::unordered_set<int32_t>> mu(range.extent_data()[1]);
      std::vector<std::unordered_set<int32_t>> Y(range.extent_data()[1]);

      for (auto i = lo[1]; i != up[1]; ++i) {
        for (auto X = lo[0]; X != up[0]; ++X) {
          for (auto n = lo[2]; n != up[2]; ++n) {
            if (norms(X, i, n) > th_) {
              Y[i].insert(X);
              mu[i].insert(n);
            }
          }
        }
      }

      TA::Tensor<float> t(range, 0.0);

      for (auto i = lo[1]; i != up[1]; ++i) {
        for (auto ye : Y[i]) {
          for (auto n : mu[i]) {
            t(ye, i, n) = std::numeric_limits<float>::max();
          }
        }
      }

      return t;
    }
  };

  array_type compute_K(array_type const &D, array_type const &C) {
    auto &world = C_df_.world();

    constexpr bool compress_C = false;
    auto dC = TA::to_new_tile_type(
        C, tensor::TaToDecompTensor(clr_thresh_, compress_C));
    dC.truncate();
    TA::foreach_inplace(dC, [&](typename decltype(dC)::value_type &tile_t) {

      auto &t = tile_t.tile().tensor(0);
      auto range = t.range();
      auto lo = range.lobound_data();
      auto up = range.upbound_data();

      if (t.norm() < 1e-3) {
        for (auto mu = lo[0]; mu < up[0]; ++mu) {
          for (auto i = lo[1]; i < up[1]; ++i) {
            if (t(mu, i) < mo_cut_thresh_) {
              t(mu, i) = 0;
            }
          }
        }
      }

      return t.norm();
    });

    auto Umo_storage = detail::array_storage(dC);
    u_mo_storages_.push_back(detail::array_storage(dC));

    darray_type C_mo, dL, F_df;

    auto cmo0 = mpqc::fenced_now(world);
    C_mo("X, i, mu") = C_df_("X, mu, nu") * dC("nu, i");
    C_mo.truncate();
    auto cmo1 = mpqc::fenced_now(world);
    c_mo_storages_.push_back(detail::array_storage(C_mo));

    TA::SparseShape<float> forced_shape;
    if (force_shape_) {
      forced_shape = C_mo.shape().transform(Rdf_shape{force_cut_thresh_});
    }

    auto edf0 = mpqc::fenced_now(world);
    if (force_shape_) {
      F_df("X, i, mu") = (E_("X, mu, nu") * dC("nu,i")).set_shape(forced_shape);
    } else {
      F_df("X, i, mu") = E_("X, mu, nu") * dC("nu,i");
    }
    F_df.truncate();
    auto edf1 = mpqc::fenced_now(world);
    e_df_storages_.push_back(detail::array_storage(F_df));

    auto rdf0 = mpqc::fenced_now(world);
    if (force_shape_) {
      F_df("X, i, mu") =
          F_df("X, i, mu") -
          (0.5 * M_("X,Y") * C_mo("Y, i, mu")).set_shape(forced_shape);
    } else {
      F_df("X, i, mu") = F_df("X, i , mu") - 0.5 * M_("X,Y") * C_mo("Y, i, mu");
    }
    F_df.truncate();
    auto rdf1 = mpqc::fenced_now(world);
    f_df_storages_.push_back(detail::array_storage(F_df));

    auto l0 = mpqc::fenced_now(world);
    dL("mu, nu") = C_mo("X, i, mu") * F_df("X, i, nu");
    dL.truncate();
    auto l1 = mpqc::fenced_now(world);

    array_type K;
    auto k0 = mpqc::fenced_now(world);
    auto L = TA::to_new_tile_type(dL, tensor::DecompToTaTensor{});
    K("mu, nu") = L("mu, nu") + L("nu, mu");
    auto k1 = mpqc::fenced_now(world);

    c_mo_times_.push_back(mpqc::duration_in_s(cmo0, cmo1));
    e_df_times_.push_back(mpqc::duration_in_s(edf0, edf1));
    r_df_times_.push_back(mpqc::duration_in_s(rdf0, rdf1));
    dl_times_.push_back(mpqc::duration_in_s(l0, l1));
    dl_to_k_times_.push_back(mpqc::duration_in_s(k0, k1));

    return K;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINEAR_CADF_BUILDER_H_
