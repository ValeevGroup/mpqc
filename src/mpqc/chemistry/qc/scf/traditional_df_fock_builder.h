#pragma once
#ifndef MPQC_SCF_TRADITIONALDFFOCKBUILDER_H
#define MPQC_SCF_TRADITIONALDFFOCKBUILDER_H

#include "../../../../../common/namespaces.h"
#include "../../../../../include/tiledarray.h"

#include "../../../../../ta_routines/array_to_eigen.h"
#include "../../../../../utility/time.h"
#include "../../../../../utility/vector_functions.h"

#include <mpqc/chemistry/qc/scf/builder.h>

#include <vector>

namespace mpqc {
namespace scf {

template <typename Integral>
class DFFockBuilder : public FockBuilder {
  private:
    array_type L_inv_; // Metric Cholesky inverse
    Integral eri3_;

    std::vector<double> j_times_;
    std::vector<double> w_times_;
    std::vector<double> k_times_;

  public:
    /*! \brief DFFockBuilder constructor takes a metric matrix
     *
     * This is to avoid forcing the user to construct the
     * inverse sqrt or cholesky decomposition of the
     * metric and to allow the behavior of the FockBuilder
     * to change without requiring user code to change.
     */
    DFFockBuilder(array_type const &M, Integral const &eri3) : eri3_(eri3) {
        auto M_eig = array_ops::array_to_eigen(M);

        MatrixD L_inv_eig
              = MatrixD(Eig::LLT<MatrixD>(M_eig).matrixL()).inverse();

        auto tr_M = M.trange().data()[0];

        L_inv_ = array_ops::eigen_to_array<TA::TensorD>(M.get_world(),
                                                        L_inv_eig, tr_M, tr_M);
    }


    const array_type &inv() const {
        return L_inv_;
    }

/*! \brief This builder requires the user to compute coefficients
     *
     * Integral is a type that can be used in a TiledArray expression, the
     * template is to allow for Direct Integral wrappers or other options.
     */
    array_type operator()(array_type const &D, array_type const &C) override {
        auto &world = D.get_world();

        auto w0 = mpqc_time::fenced_now(world);
        array_type W;
        W("X, rho, i") = L_inv_("X,Y") * (eri3_("Y, rho, sig") * C("sig, i"));
        auto w1 = mpqc_time::fenced_now(world);

        // Make J
        array_type J;
        J("mu, nu") = eri3_("Z, mu, nu")
                      * (L_inv_("X, Z") * (W("X, rho, i") * C("rho, i")));
        auto j1 = mpqc_time::fenced_now(world);

        // Permute W
        W("X, i, rho") = W("X, rho, i");

        array_type K;
        K("mu, nu") = W("X, i, mu") * W("X, i, nu");
        auto k1 = mpqc_time::fenced_now(world);

        w_times_.push_back(mpqc_time::duration_in_s(w0, w1));
        j_times_.push_back(mpqc_time::duration_in_s(w1, j1));
        k_times_.push_back(mpqc_time::duration_in_s(j1, k1));

        // Make and return G
        array_type G;
        G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");

        return G;
    }

    void print_iter(std::string const &leader) override {
        auto &world = L_inv_.get_world();

        if (world.rank() == 0) {
            std::cout << leader << "DF Fock builder:\n" << leader
                      << "\tW time: " << w_times_.back() << "\n" << leader
                      << "\tJ time: " << j_times_.back() << "\n" << leader
                      << "\tK time: " << k_times_.back() << "\n" << leader
                      << "\tTotal exchange time: "
                      << w_times_.back() + k_times_.back() << "\n";
        }
    }

    rapidjson::Value results(rapidjson::Document &d) override {
        rapidjson::Value fock_builder(rapidjson::kObjectType);
        fock_builder.AddMember("Type", "DFFockBuilder", d.GetAllocator());

        auto w_build = utility::vec_avg(w_times_);
        auto j_build = utility::vec_avg(j_times_);
        auto k_build = utility::vec_avg(k_times_);
        auto X_avg = w_build + k_build;

        fock_builder.AddMember("Avg W Time", w_build, d.GetAllocator());
        fock_builder.AddMember("Avg J Time", j_build, d.GetAllocator());
        fock_builder.AddMember("Avg K Time", k_build, d.GetAllocator());
        fock_builder.AddMember("Avg Total K Time", X_avg, d.GetAllocator());

        return fock_builder;
    }
};

} // namespace scf
} // namespace mpqc


#endif // MPQC_SCF_TRADITIONALDFFOCKBUILDER_H
