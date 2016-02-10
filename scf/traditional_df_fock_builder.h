#pragma once
#ifndef MPQC_SCF_TRADITIONALDFFOCKBUILDER_H
#define MPQC_SCF_TRADITIONALDFFOCKBUILDER_H

#include "../common/namespaces.h"
#include "../include/tiledarray.h"

namespace mpqc {
namespace scf {

class DFFockBuilder {
  public:
    using array_type = TA::TSpArrayD;

  private:
    array_type L_inv_; // Metric Cholesky inverse

  public:
    /*! \brief DFFockBuilder constructor takes a metric matrix
     *
     * This is to avoid forcing the user to construct the
     * inverse sqrt or cholesky decomposition of the
     * metric and to allow the behavior of the FockBuilder
     * to change without requiring user code to change.
     */
    DFFockBuilder(array_type const &M);

    /*! \brief This builder requires the user to compute coefficients
     *
     * Integral is a type that can be used in a TiledArray expression, the
     * template is to allow for Direct Integral wrappers or other options.
     */
    template <typename Integral>
    array_type operator()(Integral const &eri3, array_type const &C_mo) {
        array_type W;
        W("X, rho, i") = L_inv_("X,Y") * (eri3("Y, rho, sig") * C_mo("sig, i"));

        // Make J
        array_type J;
        J("mu, nu") = eri3("Z, mu, nu")
                      * (L_inv_("X, Z") * (W("X, rho, i") * C_mo("rho, i")));
        
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
};

} // namespace scf
} // namespace mpqc


#endif // MPQC_SCF_TRADITIONALDFFOCKBUILDER_H
