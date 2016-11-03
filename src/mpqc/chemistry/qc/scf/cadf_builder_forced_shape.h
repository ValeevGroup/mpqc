
#ifndef MPQC_SCF_CADFBUILDERFOCEDSHAPE_H
#define MPQC_SCF_CADFBUILDERFOCEDSHAPE_H

#include <tiledarray.h>


#include "mpqc/util/misc/time.h"
#include "mpqc/math/external/tiledarray/array_info.h"

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/tile.h"
#include "mpqc/math/tensor/clr/tensor_transforms.h"

#include "mpqc/math/external/eigen/eigen.h"

#include <mpqc/chemistry/qc/scf/builder.h>

#include <vector>
#include <iostream>

namespace mpqc {
namespace scf {

class CADFForcedShapeFockBuilder : public FockBuilder {
  public:
    using dtile_type = tensor::Tile<tensor::DecomposedTensor<double>>;
    using darray_type = TA::DistArray<dtile_type, TA::SparsePolicy>;

  private:
    darray_type B_;    // CADF fitting coeffs
    darray_type C_df_; // CADF fitting coeffs
    darray_type G_df_; // Exchange Temp

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
    template <typename Integral>
    CADFForcedShapeFockBuilder(array_type const &M, Integral const &eri3,
                               darray_type const &C_df, darray_type const &G_df,
                               double clr_thresh, double j_clr_thresh)
            : FockBuilder(), C_df_(C_df), G_df_(G_df), clr_thresh_(clr_thresh) {
        auto &world = C_df_.world();

        auto l0 = mpqc::fenced_now(world);
        auto M_eig = array_ops::array_to_eigen(M);

        RowMatrixXd L_inv_eig
              = RowMatrixXd(Eigen::LLT<RowMatrixXd>(M_eig).matrixL()).inverse();

        auto tr_M = M.trange().data()[0];

        auto L_inv
              = array_ops::eigen_to_array<TA::TensorD>(M.world(), L_inv_eig,
                                                       tr_M, tr_M);

        constexpr auto compress_L = false;
        auto dL_inv = TA::to_new_tile_type(
              L_inv, tensor::TaToDecompTensor(j_clr_thresh, compress_L));

        auto l1 = mpqc::fenced_now(world);
        l_inv_time_ = mpqc::duration_in_s(l0, l1);

        if (world.rank() == 0) {
            std::cout << "L_inv of Metric time: " << l_inv_time_ << std::endl;
        }

        auto B0 = mpqc::fenced_now(world);
        const auto old_compress = tensor::detail::recompress;
        tensor::detail::recompress = true;

        B_("X, mu, nu") = dL_inv("X, Y") * eri3("Y, mu, nu");
        minimize_storage(B_, j_clr_thresh);

        auto B1 = mpqc::fenced_now(world);
        tensor::detail::recompress = old_compress;
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

    ~CADFForcedShapeFockBuilder() = default;

    array_type operator()(array_type const &D, array_type const &C) override {
        array_type G;
        G("mu, nu") = 2 * compute_J(D)("mu, nu") - compute_K(C)("mu, nu");
        return G;
    }

    void print_iter(std::string const &leader) override {
        auto &world = C_df_.world();
        if (world.rank() == 0) {
            std::cout << leader << "CADF Forced Shape Builder:\n" << leader
                      << "\tJ time: " << j_times_.back() << "\n" << leader
                      << "\tC mo time: " << c_mo_times_.back() << "\n" << leader
                      << "\tC mo sparse storage: " << c_mo_storages_.back()[1]
                      << " GB\n" << leader
                      << "\tC mo clr storage: " << c_mo_storages_.back()[2]
                      << " GB\n" << leader
                      << "\tF df time: " << f_df_times_.back() << "\n" << leader
                      << "\tF df storage: " << f_df_storages_.back()[1]
                      << " GB\n" << leader
                      << "\tF df clr storage: " << f_df_storages_.back()[2]
                      << " GB\n" << leader << "\tdL time: " << dl_times_.back()
                      << "\n" << leader
                      << "\tdL to K time: " << dl_to_k_times_.back() << "\n";
            auto total_k_time = c_mo_times_.back() + f_df_times_.back()
                                + dl_times_.back() + dl_to_k_times_.back();
            std::cout << leader << "\ttotal K time: " << total_k_time << "\n";
        }
    }

    rapidjson::Value results(rapidjson::Document &d) override {
        rapidjson::Value fock_builder(rapidjson::kObjectType);
        fock_builder.AddMember("Type", "CADFForcedShapeFockBuilder",
                               d.GetAllocator());
        auto j_build = utility::vec_avg(j_times_);
        auto c_mo_build = utility::vec_avg(c_mo_times_);
        auto f_df_build = utility::vec_avg(f_df_times_);
        auto dl_build = utility::vec_avg(dl_times_);
        auto k_build = utility::vec_avg(dl_to_k_times_);

        fock_builder.AddMember("Avg J Build Time", j_build, d.GetAllocator());
        fock_builder.AddMember("Avg C_mo Build Time", c_mo_build,
                               d.GetAllocator());
        fock_builder.AddMember("Avg F_df Build Time", f_df_build,
                               d.GetAllocator());
        fock_builder.AddMember("Avg L Build Time", dl_build, d.GetAllocator());
        fock_builder.AddMember("Avg K from L Build Time", k_build,
                               d.GetAllocator());
        fock_builder.AddMember("Avg Total K Build Time",
                               k_build + dl_build + f_df_build + c_mo_build,
                               d.GetAllocator());

        auto c_mo_storage = utility::vec_avg(c_mo_storages_);
        fock_builder.AddMember("C_mo Dense Storage", c_mo_storage[0],
                               d.GetAllocator());
        fock_builder.AddMember("C_mo Avg Sparse Storage", c_mo_storage[1],
                               d.GetAllocator());
        fock_builder.AddMember("C_mo Avg CLR Storage", c_mo_storage[2],
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
        auto dD = TA::to_new_tile_type(D, tensor::TaToDecompTensor(clr_thresh_,
                                                                   compress_D));

        darray_type dJ;
        auto j0 = mpqc::fenced_now(world);

        dJ("mu, nu") = B_("X,mu,nu") * (B_("X,r,s") * dD("r,s"));

        auto j1 = mpqc::fenced_now(world);
        j_times_.push_back(mpqc::duration_in_s(j0, j1));

        auto J = TA::to_new_tile_type(dJ, tensor::DecompToTaTensor{});

        return J;
    }

    array_type compute_K(array_type const &C) {
        auto &world = C_df_.world();

        constexpr bool compress_C = false;
        auto dC = TA::to_new_tile_type(C, tensor::TaToDecompTensor(clr_thresh_,
                                                                   compress_C));

        darray_type C_mo, dL, F_df;
        auto cmo0 = mpqc::fenced_now(world);
        C_mo("X,  i, mu") = C_df_("X, mu, nu") * dC("nu, i");
        C_mo.truncate();
        auto cmo1 = mpqc::fenced_now(world);
        c_mo_storages_.push_back(detail::array_storage(C_mo));

        auto fdf0 = mpqc::fenced_now(world);
        F_df("X, i, mu")
              = (G_df_("X,mu, nu") * dC("nu,i")).set_shape(C_mo.shape());
        F_df.truncate();
        auto fdf1 = mpqc::fenced_now(world);
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
        f_df_times_.push_back(mpqc::duration_in_s(fdf0, fdf1));
        dl_times_.push_back(mpqc::duration_in_s(l0, l1));
        dl_to_k_times_.push_back(mpqc::duration_in_s(k0, k1));

        return K;
    }
};

} // namespace scf
} // namespace mpqc


#endif // MPQC_SCF_CADFBUILDERFOCEDSHAPE_H
