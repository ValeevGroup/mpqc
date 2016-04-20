#pragma once
#ifndef MPQC_SCF_LINEARCADFBUILDER_H
#define MPQC_SCF_LINEARCADFBUILDER_H

#include "../common/namespaces.h"
#include "../include/tiledarray.h"
#include "../utility/time.h"
#include "../utility/array_info.h"
#include "../utility/vector_functions.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/mpqc_tile.h"
#include "../tensor/tensor_transforms.h"

#include "../ta_routines/array_to_eigen.h"
#include "../ta_routines/minimize_storage.h"
#include "../ta_routines/sqrt_inv.h"

#include "builder.h"
#include "ta_shape_tracker.h"

#include <vector>
#include <iostream>

namespace mpqc {
namespace scf {

template <typename Integral>
class ONCADFFockBuilder : public FockBuilder {
  public:
    using dtile_type = tensor::Tile<tensor::DecomposedTensor<double>>;
    using darray_type = TA::DistArray<dtile_type, SpPolicy>;

  private:
    Integral E_;             // Direct three center integrals
    darray_type M_;          // Fitting Metric
    darray_type M_sqrt_inv_; // Fitting Metric
    darray_type C_df_;       // CADF fitting coeffs

    double clr_thresh_;

    std::vector<double> j_times_;
    std::vector<double> c_mo_times_;
    std::vector<double> e_df_times_;
    std::vector<double> r_df_times_;
    std::vector<double> f_df_times_;

    std::vector<double> dl_times_;
    std::vector<double> dl_to_k_times_;

    std::vector<std::array<double, 3>> c_mo_storages_;
    std::vector<std::array<double, 3>> e_df_storages_;
    std::vector<std::array<double, 3>> r_df_storages_;
    std::vector<std::array<double, 3>> f_df_storages_;

    ShapeTracker shape_tracker_;
    std::vector<ShapeTracker> shape_tracker_iters_;

  public:
    ONCADFFockBuilder(darray_type const &M, Integral const &eri3,
                      darray_type const &C_df, double clr_thresh)
            : FockBuilder(),
              E_(eri3),
              M_(M),
              C_df_(C_df),
              clr_thresh_(clr_thresh) {
        auto &world = C_df_.get_world();
        bool compress_M = false;

        auto l0 = mpqc_time::fenced_now(world);
        M_sqrt_inv_ = array_ops::inverse_sqrt(M_);
        auto l1 = mpqc_time::fenced_now(world);
        auto l_inv_time_ = mpqc_time::duration_in_s(l0, l1);

        if (world.rank() == 0) {
            std::cout << "Sqrt Inv of Metric time: " << l_inv_time_
                      << std::endl;
        }
    }

    ~ONCADFFockBuilder() = default;

    array_type operator()(array_type const &D, array_type const &C) override {
        array_type G;
        G("mu, nu") = 2 * compute_J(D)("mu, nu") - compute_K(C)("mu, nu");
        return G;
    }

    void print_iter(std::string const &leader) override {
        auto &world = C_df_.get_world();
        if (world.rank() == 0) {
            std::cout << leader << "CADF Builder:\n" << leader
                      << "\tStorages:\n" << leader
                      << "\t\tFully Dense: " << c_mo_storages_.back()[0]
                      << " GB\n" << leader
                      << "\t\tC mo sparse storage: " << c_mo_storages_.back()[1]
                      << " GB\n" << leader
                      << "\t\tC mo clr storage: " << c_mo_storages_.back()[2]
                      << " GB\n" << leader
                      << "\t\tE df storage: " << e_df_storages_.back()[1]
                      << " GB\n" << leader
                      << "\t\tE df clr storage: " << e_df_storages_.back()[2]
                      << " GB\n" << leader
                      << "\t\tR df storage: " << r_df_storages_.back()[1]
                      << " GB\n" << leader
                      << "\t\tR df clr storage: " << r_df_storages_.back()[2]
                      << " GB\n" << leader
                      << "\t\tF df storage: " << f_df_storages_.back()[1]
                      << " GB\n" << leader
                      << "\t\tF df clr storage: " << f_df_storages_.back()[2]
                      << " GB\n" << leader << "\tJ time: " << j_times_.back()
                      << "\n" << leader << "\tC mo time: " << c_mo_times_.back()
                      << "\n" << leader << "\tE df time: " << e_df_times_.back()
                      << "\n" << leader << "\tR df time: " << r_df_times_.back()
                      << "\n" << leader << "\tF df time: " << f_df_times_.back()
                      << "\n" << leader << "\tdL time: " << dl_times_.back()
                      << "\n" << leader
                      << "\tdL to K time: " << dl_to_k_times_.back() << "\n";
            auto total_k_time = c_mo_times_.back() + f_df_times_.back()
                                + dl_times_.back() + dl_to_k_times_.back();
            std::cout << leader << "\ttotal K time: " << total_k_time << "\n";
        }
    }

    rapidjson::Value results(rapidjson::Document &d) override {}

  private:
    array_type compute_J(array_type const &D) {
        auto &world = C_df_.get_world();

        constexpr bool compress_D = false;
        auto dD = TA::to_new_tile_type(D, tensor::TaToDecompTensor(clr_thresh_,
                                                                   compress_D));

        darray_type dJ;
        auto j0 = mpqc_time::fenced_now(world);

        dJ("mu, nu") = E_("X,mu,nu")
                       * (M_sqrt_inv_("X,Z")
                          * (M_sqrt_inv_("Z, Y") * (E_("Y,r,s") * dD("r,s"))));

        auto j1 = mpqc_time::fenced_now(world);
        j_times_.push_back(mpqc_time::duration_in_s(j0, j1));

        auto J = TA::to_new_tile_type(dJ, tensor::DecompToTaTensor{});

        return J;
    }

    array_type compute_K(array_type const &C) {
        auto &world = C_df_.get_world();

        constexpr bool compress_C = false;
        auto dC = TA::to_new_tile_type(C, tensor::TaToDecompTensor(clr_thresh_,
                                                                   compress_C));

        darray_type C_mo, dL, E_df, R_df, F_df;

        auto cmo0 = mpqc_time::fenced_now(world);
        C_mo("X, i, mu") = C_df_("X, mu, nu") * dC("nu, i");
        C_mo.truncate();
        auto cmo1 = mpqc_time::fenced_now(world);
        c_mo_storages_.push_back(utility::array_storage(C_mo));

        auto edf0 = mpqc_time::fenced_now(world);
        E_df("X, i, mu")
              = (E_("X, mu, nu") * dC("nu,i")).set_shape(C_mo.get_shape());
        E_df.truncate();
        auto edf1 = mpqc_time::fenced_now(world);
        e_df_storages_.push_back(utility::array_storage(E_df));

        auto rdf0 = mpqc_time::fenced_now(world);
        R_df("X, i, mu") = M_("X,Y") * C_mo("Y, i, mu");
        R_df.truncate();
        auto rdf1 = mpqc_time::fenced_now(world);
        r_df_storages_.push_back(utility::array_storage(R_df));

        auto fdf0 = mpqc_time::fenced_now(world);
        F_df("X, i, mu") = E_df("X, i, mu") - 0.5 * R_df("X, i, mu");
        F_df.truncate();
        auto fdf1 = mpqc_time::fenced_now(world);
        f_df_storages_.push_back(utility::array_storage(F_df));

        // shape_tracker_.update(C_mo.get_shape().data(),
        // F_df.get_shape().data());
        // ShapeTracker this_iter;
        // this_iter.update(C_mo.get_shape().data(), F_df.get_shape().data());
        // shape_tracker_iters_.emplace_back(std::move(this_iter));

        auto l0 = mpqc_time::fenced_now(world);
        dL("mu, nu") = C_mo("X, i, mu") * F_df("X, i, nu");
        dL.truncate();
        auto l1 = mpqc_time::fenced_now(world);

        array_type K;
        auto k0 = mpqc_time::fenced_now(world);
        auto L = TA::to_new_tile_type(dL, tensor::DecompToTaTensor{});
        K("mu, nu") = L("mu, nu") + L("nu, mu");
        auto k1 = mpqc_time::fenced_now(world);

        c_mo_times_.push_back(mpqc_time::duration_in_s(cmo0, cmo1));
        e_df_times_.push_back(mpqc_time::duration_in_s(edf0, edf1));
        r_df_times_.push_back(mpqc_time::duration_in_s(rdf0, rdf1));
        f_df_times_.push_back(mpqc_time::duration_in_s(fdf0, fdf1));
        dl_times_.push_back(mpqc_time::duration_in_s(l0, l1));
        dl_to_k_times_.push_back(mpqc_time::duration_in_s(k0, k1));

        return K;
    }
};

} // namespace scf
} // namespace mpqc


#endif // MPQC_SCF_LINEARCADFBUILDER_H
