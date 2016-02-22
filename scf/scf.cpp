#include "scf.h"
#include "../integrals/integrals.h"
#include "../utility/time.h"
#include "diagonalize_for_coffs.hpp"
#include "../ta_routines/array_to_eigen.h"
#include "../ta_routines/minimize_storage.h"

namespace mpqc {
namespace scf {

double ClosedShellSCF::energy() {
    return repulsion_
           + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
}

bool ClosedShellSCF::solve(int64_t max_iters, double thresh) {
    auto &world = F_.get_world();

    if (world.rank() == 0) {
        std::cout << "Starting SCF:\n"
                  << "\tThreshold: " << thresh << "\n"
                  << "\tMaximum number of iterations: " << max_iters << "\n";
    }

    auto iter = 0;
    auto error = std::numeric_limits<double>::max();
    auto rms_error = std::numeric_limits<double>::max();
    auto old_energy = energy();
    const double volume = F_.trange().elements().volume();

    while (iter < max_iters && (thresh < (error / old_energy)
                                || thresh < (rms_error / volume))) {

        auto s0 = mpqc_time::fenced_now(world);

        build_F();

        auto current_energy = energy();
        error = std::abs(old_energy - current_energy);
        old_energy = current_energy;

        array_type Grad;
        Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                      - S_("i,k") * D_("k,l") * F_("l,j");

        rms_error = Grad("i,j").norm().get();
        diis_.extrapolate(F_, Grad);

        compute_density(occ_);
        auto s1 = mpqc_time::fenced_now(world);

        scf_times_.push_back(mpqc_time::duration_in_s(s0, s1));

        if (world.rank() == 0) {
            std::cout << "\titeration: " << iter << "\n"
                      << "\t\tEnergy: " << old_energy << "\n"
                      << "\t\tabs(Energy Change)/energy: "
                      << (error / std::abs(old_energy)) << "\n"
                      << "\t\t(Gradient Norm)/n^2: " << (rms_error / volume)
                      << "\n";
        }
        builder_->print_iter("\t\t");
        ++iter;
    }

    if(iter > max_iters || (thresh > (error / old_energy)
                                || thresh > (rms_error / volume))) {
        return false;
    } else {
        return true;
    }

}

void ClosedShellSCF::compute_density(int64_t occ) {
    auto F_eig = array_ops::array_to_eigen(F_);
    auto S_eig = array_ops::array_to_eigen(S_);

    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);
    decltype(S_eig) C = es.eigenvectors().leftCols(occ);
    auto tr_ao = S_.trange().data()[0];

    auto occ_nclusters = (occ_ < 10) ? occ_ : 10;
    auto tr_occ = scf::tr_occupied(occ_nclusters, occ_);

    C_ = array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C, tr_ao,
                                                tr_occ);

    if (TcutC_ != 0) {
        ta_routines::minimize_storage(C_, TcutC_);
    }

    D_("i,j") = C_("i,k") * C_("j,k");
}

void ClosedShellSCF::build_F() {
    F_("i,j") = H_("i,j") + builder_->operator()(D_, C_)("i,j");
}


} // namespace scf
} // namespace mpqc
