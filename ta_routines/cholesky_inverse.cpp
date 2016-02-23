#include "cholesky_inverse.h"
#include "../include/eigen.h"
#include "array_to_eigen.h"

namespace mpqc {
namespace array_ops {

TA::DistArray<TA::TensorD, SpPolicy>
cholesky_inverse(TA::DistArray<TA::TensorD, SpPolicy> const &A) {
    auto &world = A.get_world();

    auto A_eig = array_to_eigen(A);
    Eig::LLT<decltype(A_eig)> llt(A_eig);

    auto L_inv_eig = decltype(A_eig)(llt.matrixL()).inverse();

    
    auto tr_A0 = A.trange().data()[0];
    auto tr_A1 = A.trange().data()[1];

    return eigen_to_array<TA::TensorD>(world, L_inv_eig, tr_A0, tr_A1);
}

} // namespace ta_routines
} // namespace mpqc
