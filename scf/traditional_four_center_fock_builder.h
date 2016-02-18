#pragma once
#ifndef MPQC_SCF_TRADITIONALFOURCENTERFOCKBUILDER_H
#define MPQC_SCF_TRADITIONALFOURCENTERFOCKBUILDER_H

#include "../common/namespaces.h"
#include "../include/tiledarray.h"

namespace mpqc {
namespace scf {

class FourCenterBuilder {
  public:
    using array_type = TA::TSpArrayD;

  public:
    /*! \brief This builder requires the user to compute coefficients
     *
     * Integral is a type that can be used in a TiledArray expression, the
     * template is to allow for Direct Integral wrappers or other options.
     */
    template <typename Integral>
    array_type operator()(Integral const &eri4, array_type const &C_mo) {
        array_type D;
        D("mu, nu") = C_mo("mu, i") * C_mo("nu, i");
        D.truncate();
       
        // Make J
        array_type J;
        J("mu, nu") = eri4("mu, nu, rho, sig") * D("rho, sig");
        
        // Make K
        array_type K;
        K("mu, nu") = eri4("mu, rho, nu, sig") * D("rho, sig");
        
        // Make and return G
        array_type G;
        G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");

        return G;
    }
};

} // namespace scf
} // namespace mpqc
#endif // MPQC_SCF_TRADITIONALFOURCENTERFOCKBUILDER_H
