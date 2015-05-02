#pragma once
#ifndef TCC_PURIFICATION_PURIFICATIONDEVEL_H
#define TCC_PURIFICATION_PURIFICATIONDEVEL_H

#include "../include/eigen.h"
#include "../include/tiledarray.h"
#include "../utility/parallel_print.h"
#include "../utility/time.h"

#include "eigen_value_estimation.h"
#include "diagonal_array.h"

namespace tcc {
namespace pure {

template <typename T, typename TileType, typename Policy>
class OrthTraceResettingPurifier {
  public:
    using ArrayType = TiledArray::Array<T, 2, TileType, Policy>;

    OrthTraceResettingPurifier(ArrayType const &sqrt_inv)
            : sqrt_inv_(sqrt_inv), I(create_diagonal_matrix(sqrt_inv, 1.0)) {}

    double trace(ArrayType const &A) const { return A("i,j").trace(); }

    ArrayType operator()(ArrayType const &Fao, std::size_t occ) const {
        // F = Z F_{ao} Z^{T}
        ArrayType F;
        F("i,j") = sqrt_inv_("i,k") * Fao("k,l") * sqrt_inv_("l,j");

        auto eig_pair = eval_guess(F);
        auto emax = eig_pair[1];
        auto emin = eig_pair[0];
        auto scale = 1.0 / (emax - emin);

        // D_{o} = \frac{emax I - F}{emax - emin}
        ArrayType D;
        D("i,j") = scale * (emax * I("i,j") - F("i,j"));
        auto tr = trace(D);

        occ = occ / 2;
        auto iter = 1;
        ArrayType D2;
        auto error = std::abs(tr - occ);
        while (error >= 1e-12 && iter <= 100) {
            // Compute D2
            D2("i,j") = D("i,k") * D("k,j");
            if (tr > occ) {
                D = D2;
            } else {
                D("i,j") = 2 * D("i,j") - D2("i,j");
            }
            D.truncate();
            tr = trace(D);
            error = std::abs(tr - occ);
            ++iter;
        }

        // D_{ao} = Z^{T} D Z
        D("i,j") = sqrt_inv_("i,k") * D("k,l") * sqrt_inv_("l,j");
        D.truncate(); // necessary since norms blow up otherwise
        return D;
    }


  private:
    const ArrayType sqrt_inv_;
    const ArrayType I;
};

template <typename T, typename TileType, typename Policy>
OrthTraceResettingPurifier<T, TileType, Policy> make_orthogonal_tr_reset_pure(
      TiledArray::Array<T, 2, TileType, Policy> const &sqrt_inv) {
    return OrthTraceResettingPurifier<T, TileType, Policy>{sqrt_inv};
}

} // namespace pure
} // namespace tcc

#endif /* end of include guard: TCC_PURIFICATION_PURIFICATIONDEVEL_H */
