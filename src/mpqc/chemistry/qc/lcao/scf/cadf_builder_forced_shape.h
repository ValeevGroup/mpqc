
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CADF_BUILDER_FORCED_SHAPE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CADF_BUILDER_FORCED_SHAPE_H_

#include <tiledarray.h>

#include "mpqc/math/external/tiledarray/array_info.h"
#include "mpqc/util/misc/time.h"

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/tensor_transforms.h"
#include "mpqc/math/tensor/clr/tile.h"

#include "mpqc/math/external/eigen/eigen.h"

#include "mpqc/chemistry/qc/lcao/scf/builder.h"

#include <iostream>
#include <vector>

namespace mpqc {
namespace lcao {
namespace scf {

class CADFForcedShapeFockBuilder : public FockBuilder {
 public:
  using dtile_type = Tile<DecomposedTensor<double>>;
  using darray_type = TA::DistArray<dtile_type, TA::SparsePolicy>;

 private:
  darray_type B_;     // CADF fitting coeffs
  darray_type C_df_;  // CADF fitting coeffs
  darray_type G_df_;  // Exchange Temp

  double clr_thresh_;

  double l_inv_time_ = 0;
  double B_time_ = 0;
  std::array<double, 3> B_storages_;

  std::vector<double> j_times_;
  std::vector<double> c_mo_times_;
  std::vector<double> f_df_times_;
  std::vector<double> dl_times_;
  std::vector<double> dl_to_k_times_;

  std::vector<std::array<double, 3>> c_mo_storages_;
  std::vector<std::array<double, 3>> f_df_storages_;

 public:
  template<typename Integral>
  CADFForcedShapeFockBuilder(array_type const &M, Integral const &eri3,
                             darray_type const &C_df, darray_type const &G_df,
                             double clr_thresh, double j_clr_thresh)
      : FockBuilder(), C_df_(C_df), G_df_(G_df), clr_thresh_(clr_thresh) {
    auto &world = C_df_.world();

    auto l0 = mpqc::fenced_now(world);
    auto M_eig = math::array_to_eigen(M);

    RowMatrixXd L_inv_eig =
        RowMatrixXd(Eigen::LLT<RowMatrixXd>(M_eig).matrixL()).inverse();

    auto tr_M = M.trange().data()[0];

    auto L_inv = math::eigen_to_array<TA::TensorD>(M.world(), L_inv_eig,
                                                        tr_M, tr_M);

    constexpr auto compress_L = false;
    auto dL_inv = TA::to_new_tile_type(
        L_inv, TaToDecompTensor(j_clr_thresh, compress_L));

    auto l1 = mpqc::fenced_now(world);
    l_inv_time_ = mpqc::duration_in_s(l0, l1);

    if (world.rank() == 0) {
      std::cout << "L_inv of Metric time: " << l_inv_time_ << std::endl;
    }

    auto B0 = mpqc::fenced_now(world);
    const auto old_compress = detail::recompress;
    detail::recompress = true;

    B_("X, mu, nu") = dL_inv("X, Y") * eri3("Y, mu, nu");
    minimize_storage(B_, j_clr_thresh);

    auto B1 = mpqc::fenced_now(world);
    detail::recompress = old_compress;
    B_time_ = mpqc::duration_in_s(B0, B1);

    B_storages_ = detail::array_storage(B_);
    if (world.rank() == 0) {
      std::cout << "B time: " << B_time_ << std::endl;
      std::cout << "B storage:\n"
                << "Full    = " << B_storages_[0] << "\n"
                << "Sparse  = " << B_storages_[1] << "\n"
                << "CLR     = " << B_storages_[2] << std::endl;
    }
  }

  ~CADFForcedShapeFockBuilder() {}

  array_type operator()(array_type const &D, array_type const &C) override {
    array_type G;
    G("mu, nu") = 2 * compute_J(D)("mu, nu") - compute_K(C)("mu, nu");
    return G;
  }

  void print_iter(std::string const &leader) override {
    auto &world = C_df_.world();
    if (world.rank() == 0) {
      std::cout << leader << "CADF Forced Shape Builder:\n"
                << leader << "\tJ time: " << j_times_.back() << "\n"
                << leader << "\tC mo time: " << c_mo_times_.back() << "\n"
                << leader
                << "\tC mo sparse storage: " << c_mo_storages_.back()[1]
                << " GB\n"
                << leader << "\tC mo clr storage: " << c_mo_storages_.back()[2]
                << " GB\n"
                << leader << "\tF df time: " << f_df_times_.back() << "\n"
                << leader << "\tF df storage: " << f_df_storages_.back()[1]
                << " GB\n"
                << leader << "\tF df clr storage: " << f_df_storages_.back()[2]
                << " GB\n"
                << leader << "\tdL time: " << dl_times_.back() << "\n"
                << leader << "\tdL to K time: " << dl_to_k_times_.back()
                << "\n";
      auto total_k_time = c_mo_times_.back() + f_df_times_.back() +
          dl_times_.back() + dl_to_k_times_.back();
      std::cout << leader << "\ttotal K time: " << total_k_time << "\n";
    }
  }

 private:
  array_type compute_J(array_type const &D) {
    auto &world = C_df_.world();

    constexpr bool compress_D = false;
    auto dD = TA::to_new_tile_type(
        D, TaToDecompTensor(clr_thresh_, compress_D));

    darray_type dJ;
    auto j0 = mpqc::fenced_now(world);

    dJ("mu, nu") = B_("X,mu,nu") * (B_("X,r,s") * dD("r,s"));

    auto j1 = mpqc::fenced_now(world);
    j_times_.push_back(mpqc::duration_in_s(j0, j1));

    auto J = TA::to_new_tile_type(dJ, DecompToTaTensor{});

    return J;
  }

  array_type compute_K(array_type const &C) {
    auto &world = C_df_.world();

    constexpr bool compress_C = false;
    auto dC = TA::to_new_tile_type(
        C, TaToDecompTensor(clr_thresh_, compress_C));

    darray_type C_mo, dL, F_df;
    auto cmo0 = mpqc::fenced_now(world);
    C_mo("X,  i, mu") = C_df_("X, mu, nu") * dC("nu, i");
    C_mo.truncate();
    auto cmo1 = mpqc::fenced_now(world);
    c_mo_storages_.push_back(detail::array_storage(C_mo));

    auto fdf0 = mpqc::fenced_now(world);
    F_df("X, i, mu") = (G_df_("X,mu, nu") * dC("nu,i")).set_shape(C_mo.shape());
    F_df.truncate();
    auto fdf1 = mpqc::fenced_now(world);
    f_df_storages_.push_back(detail::array_storage(F_df));

    auto l0 = mpqc::fenced_now(world);
    dL("mu, nu") = C_mo("X, i, mu") * F_df("X, i, nu");
    dL.truncate();
    auto l1 = mpqc::fenced_now(world);

    array_type K;
    auto k0 = mpqc::fenced_now(world);
    auto L = TA::to_new_tile_type(dL, DecompToTaTensor{});
    K("mu, nu") = L("mu, nu") + L("nu, mu");
    auto k1 = mpqc::fenced_now(world);

    c_mo_times_.push_back(mpqc::duration_in_s(cmo0, cmo1));
    f_df_times_.push_back(mpqc::duration_in_s(fdf0, fdf1));
    dl_times_.push_back(mpqc::duration_in_s(l0, l1));
    dl_to_k_times_.push_back(mpqc::duration_in_s(k0, k1));

    return K;
  }
};

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CADF_BUILDER_FORCED_SHAPE_H_
