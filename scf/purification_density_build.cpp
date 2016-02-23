#include "../integrals/integrals.h"
#include "../common/typedefs.h"

#include "purification_density_build.h"
#include "../ta_routines/array_to_eigen.h"
#include "../ta_routines/cholesky_inverse.h"
#include "../ta_routines/sqrt_inv.h"
#include "../ta_routines/minimize_storage.h"

#include "../tensor/decomposed_tensor_algebra.h"

#include "diagonalize_for_coffs.hpp"
#include "orbital_localization.h"
#include "clusterd_coeffs.h"


namespace mpqc {
namespace scf {

using array_type = PurificationDensityBuilder::array_type;

PurificationDensityBuilder::PurificationDensityBuilder(
      array_type const &S, std::vector<array_type> r_xyz, int64_t occ,
      int64_t nclusters, double TcutC, bool localize)
        : S_(S),
          r_xyz_ints_(r_xyz),
          occ_(occ),
          n_coeff_clusters_(nclusters),
          TcutC_(TcutC),
          localize_(localize) {
    M_inv_ = array_ops::inverse_sqrt(S_);
    I_ = array_ops::create_diagonal_matrix(S_, 1.0);
}

array_type PurificationDensityBuilder::purify(array_type const &F) {
    auto &world = F.get_world();

    array_type Fp, D, D2;
    Fp("i,j") = M_inv_("i,k") * F("k,l") * M_inv_("j,l");

    auto eig_pair = array_ops::eval_guess(Fp);
    auto emax = eig_pair[1];
    auto emin = eig_pair[0];
    auto scale = 1.0 / (emax - emin);

    D("i,j") = scale * (emax * I_("i,j") - Fp("i,j"));

    auto trace = D("i,j").trace().get();

    auto iter = 0;
    auto error = std::abs(trace - occ_);
    while (error >= 1e-13 && iter <= 100) {
        // Compute D2
        D2("i,j") = D("i,k") * D("k,j");
        if (trace > occ_) {
            D = D2;
        } else {
            D("i,j") = 2 * D("i,j") - D2("i,j");
        }
        D.truncate();
        trace = D("i,j").trace().get();
        error = std::abs(trace - occ_);
        ++iter;
    }

    D("i,j") = M_inv_("i,k") * D("k,l") * M_inv_("l,j");
    D.truncate();
    world.gop.fence();

    return D;
}

array_type PurificationDensityBuilder::orbitals(array_type const &D) {
    auto &world = D.get_world();

    auto D_eig = array_ops::array_to_eigen(D);
    tensor::algebra::piv_cholesky(D_eig);

    auto tr_ao = D.trange().data()[0];
    auto tr_occ = scf::tr_occupied(n_coeff_clusters_, occ_);

    auto Cao = array_ops::eigen_to_array<TA::TensorD>(D.get_world(), D_eig,
                                                      tr_ao, tr_occ);

    if (localize_) {
        auto U = mpqc::scf::BoysLocalization{}(Cao, r_xyz_ints_);
        Cao("mu,i") = Cao("mu,k") * U("k,i");

        auto obs_ntiles = Cao.trange().tiles().extent()[0];
        scf::clustered_coeffs(r_xyz_ints_, Cao, obs_ntiles);
    }

    if (TcutC_ != 0) {
        ta_routines::minimize_storage(Cao, TcutC_);
    }

    return Cao;
}

std::pair<array_type, array_type> PurificationDensityBuilder::
operator()(array_type const &F) {
    auto D = purify(F);
    auto C = orbitals(D);
    return std::make_pair(D, C);
}

} // namespace scf
} // namespace mpqc
