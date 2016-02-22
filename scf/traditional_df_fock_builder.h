#pragma once
#ifndef MPQC_SCF_TRADITIONALDFFOCKBUILDER_H
#define MPQC_SCF_TRADITIONALDFFOCKBUILDER_H

#include "../common/namespaces.h"
#include "../include/tiledarray.h"

#include "../ta_routines/array_to_eigen.h"

#include "builder.h"

namespace mpqc {
namespace scf {

template <typename Integral>
class DFFockBuilder : public FockBuilder {
  private:
    array_type L_inv_; // Metric Cholesky inverse
    Integral eri3_;

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


    /*! \brief This builder requires the user to compute coefficients
     *
     * Integral is a type that can be used in a TiledArray expression, the
     * template is to allow for Direct Integral wrappers or other options.
     */
    array_type operator()(array_type const &D, array_type const &C) override {
        array_type W;
        W("X, rho, i") = L_inv_("X,Y") * (eri3_("Y, rho, sig") * C("sig, i"));

        // Make J
        array_type J;
        J("mu, nu") = eri3_("Z, mu, nu")
                      * (L_inv_("X, Z") * (W("X, rho, i") * C("rho, i")));

        // Make K

        // Permute W
        W("X, i, rho") = W("X, rho, i");

        array_type K;
        K("mu, nu") = W("X, i, mu") * W("X, i, nu");

        // Make and return G
        array_type G;
        G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");

        return G;
    }

    void print_iter(std::string const &) override {}
};

} // namespace scf
} // namespace mpqc


#endif // MPQC_SCF_TRADITIONALDFFOCKBUILDER_H
