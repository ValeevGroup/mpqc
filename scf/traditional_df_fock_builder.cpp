#include "traditional_df_fock_builder.h"
#include "../ta_routines/array_to_eigen.h"

namespace mpqc {
namespace scf {

DFFockBuilder::DFFockBuilder(array_type const &M) {
    auto M_eig = array_ops::array_to_eigen(M);

    MatrixD L_inv_eig = MatrixD(Eig::LLT<MatrixD>(M_eig).matrixL()).inverse();

    auto tr_M = M.trange().data()[0];

    L_inv_ = array_ops::eigen_to_array<TA::TensorD>(M.get_world(), L_inv_eig,
                                                    tr_M, tr_M);
}

} // namespace scf
} // namespace mpqc
