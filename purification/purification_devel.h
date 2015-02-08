#pragma once
#ifndef TCC_PURIFICATION_PURIFICATIONDEVEL_H
#define TCC_PURIFICATION_PURIFICATIONDEVEL_H

#include <cstdlib>
#include "eigen_value_estimation.h"
#include <iostream>
#include "../include/eigen.h"

namespace tcc {
namespace pure {

template <typename Array>
Array initial_guess(Array const &H, Array const &S) {
    return create_eval_scaled_guess(H, S);
}

template <typename Array>
class trace_resetting_poly {
  public:
    void operator()(Array &D, Array const &S, std::size_t occ) const {

        Array DS(D.get_world(), D.trange());
        DS("i,j") = D("i,k") * S("k,j");

        auto trace = array_trace(DS);
        occ = occ / 2;
        auto iter = 1;

        while (std::abs(trace - occ) >= 1e-6 && iter <= 1000) {
            if (trace > occ) {
                D("i,j") = DS("i,k") * D("k,j");
                std::cout << "\tshrinking trace" << std::endl;
            } else {
                D("i,j") = 2 * D("i,j") - DS("i,k") * D("k,j");
                std::cout << "\tincreasing trace" << std::endl;
            }

            DS("i,j") = D("i,k") * S("k,j");

            trace = array_trace(DS);

            std::cout << "DS trace = " << trace << " trace and occ diff "
                      << std::abs(trace - occ) << std::endl;
            ++iter;
            D.get_world().gop.fence();
        }
    }

    double array_trace(Array const &A) const { return A("i,j").trace(); }
};


class purifier {
  public:
    explicit purifier(double cut = 1e-07) : cut_{cut} {}

    template <typename Array, typename Polynomial = trace_resetting_poly<Array>>
    Array operator()(Array const &H, Array const &S, std::size_t occupation,
                     Polynomial poly = Polynomial{}) {
        auto P = initial_guess(H, S);
        Array PS;
        PS("i,j") = P("i,k") * S("k,j");
        auto eig_PS = TiledArray::array_to_eigen(PS);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(eig_PS);
        std::cout << "Input evals = " << cut_ << " "
                  << es.eigenvalues().transpose() << std::endl;

        poly(P, S, occupation);

        return P;
    }

  private:
    double cut_;
}; // class purifier

} // namespace pure
} // namespace tcc

#endif /* end of include guard: TCC_PURIFICATION_PURIFICATIONDEVEL_H */
