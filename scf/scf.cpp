//
// Created by Chong Peng on 2/9/16.
//

#include "scf.h"
#include "../utility/time.h"
#include "diagonalize_for_coffs.hpp"


// member function of RHF
void mpqc::RHF::solve(int64_t max_iters, double thresh, const array_type &eri4) {
    auto iter = 0;
    auto error = std::numeric_limits<double>::max();
    auto old_energy = 0.0;

    while (iter < max_iters && thresh < error) {
        auto s0 = mpqc_time::now();
        F_.get_world().gop.fence();
        form_fock(eri4);

        auto current_energy = energy();
        error = std::abs(old_energy - current_energy);
        old_energy = current_energy;

        array_type Grad;
        Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                      - S_("i,k") * D_("k,l") * F_("l,j");

        diis_.extrapolate(F_, Grad);

        // Lastly update density
        compute_density(occ_);

        F_.get_world().gop.fence();
        auto s1 = mpqc_time::now();
        scf_times_.push_back(mpqc_time::duration_in_s(s0, s1));


        std::cout << "Iteration: " << (iter + 1)
        << " energy: " << old_energy << " error: " << error
        << std::endl;
        std::cout << "\tJ time: " << j_times_.back()
        << " s K time: " << k_times_.back()
        << " s iter time: " << scf_times_.back() << std::endl;

        ++iter;
    }
}


void mpqc::RHF::compute_density(int64_t occ) {
    auto F_eig = array_ops::array_to_eigen(F_);
    auto S_eig = array_ops::array_to_eigen(S_);

    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                               S_eig);
    decltype(S_eig) C = es.eigenvectors().leftCols(occ);
    MatrixD D_eig = C * C.transpose();

    auto tr_ao = S_.trange().data()[0];

    D_ = array_ops::eigen_to_array<TA::TensorD>(F_.get_world(), D_eig,
                                                tr_ao, tr_ao);
}

template<typename Integral>
void mpqc::RHF::form_fock(const Integral &eri4) {
    F_("i,j") = H_("i,j") + 2 * compute_j(eri4)("i,j")
                - compute_k(eri4)("i,j");
}

template<typename Integral>
mpqc::RHF::array_type mpqc::RHF::compute_j(const Integral &eri4) {
    array_type J;
    auto &world = eri4.get_world();
    world.gop.fence();
    auto j0 = mpqc_time::now();
    J("i,j") = eri4("i,j,k,l") * D_("k,l");
    world.gop.fence();
    auto j1 = mpqc_time::now();
    j_times_.push_back(mpqc_time::duration_in_s(j0, j1));

    return J;
}

template<typename Integral>
mpqc::RHF::array_type mpqc::RHF::compute_k(const Integral &eri4) {
    array_type K;
    auto &world = eri4.get_world();
    world.gop.fence();
    auto k0 = mpqc_time::now();
    K("i,j") = eri4("i,k,j,l") * D_("k,l");
    world.gop.fence();
    auto k1 = mpqc_time::now();
    k_times_.push_back(mpqc_time::duration_in_s(k0, k1));

    return K;
}


// member function of DFRHF
void mpqc::DFRHF::solve(int64_t max_iters, double thresh, const array_type &eri3) {
    if(F_.get_world().rank() == 0){
        std::cout << "Start SCF" << std::endl;
        std::cout << "Convergence : " << thresh << std::endl;
        std::cout << "Max Iteration : " << max_iters << std::endl;
    }

    auto iter = 0;
    auto error = std::numeric_limits<double>::max();
    auto old_energy = 0.0;

    while (iter < max_iters && thresh < error) {
        auto s0 = mpqc_time::now();
        F_.get_world().gop.fence();
        form_fock(eri3);

        auto current_energy = energy();
        error = std::abs(old_energy - current_energy);
        old_energy = current_energy;

        array_type Grad;
        Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                      - S_("i,k") * D_("k,l") * F_("l,j");

        diis_.extrapolate(F_, Grad);

        // Lastly update density
        compute_density(occ_);

        F_.get_world().gop.fence();
        auto s1 = mpqc_time::now();
        scf_times_.push_back(mpqc_time::duration_in_s(s0, s1));

        if(F_.get_world().rank() == 0){
            std::cout << "Iteration: " << (iter + 1)
            << " energy: " << old_energy << " error: " << error
            << std::endl;
            std::cout << "\tW time: " << w_times_.back() << std::endl;
            std::cout << "\tJ time: " << j_times_.back()
            << " s K time: " << k_times_.back()
            << " s iter time: " << scf_times_.back() << std::endl;
        }
        ++iter;
    }
}

void mpqc::DFRHF::compute_density(int64_t occ) {
    auto F_eig = array_ops::array_to_eigen(F_);
    auto S_eig = array_ops::array_to_eigen(S_);

    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                               S_eig);
    decltype(S_eig) C = es.eigenvectors().leftCols(occ);
    auto tr_ao = S_.trange().data()[0];

    auto occ_nclusters = (occ_ < 10) ? occ_ : 10;
    auto tr_occ = scf::tr_occupied(occ_nclusters, occ_);

    C_ = array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C, tr_ao, tr_occ);

    D_("i,j") = C_("i,k") * C_("j,k");
}

template<typename Integral>
void mpqc::DFRHF::form_fock(const Integral &eri3) {
    auto &world = F_.get_world();

    world.gop.fence();
    auto w0 = mpqc_time::now();
    TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> W;
    W("X, mu, i") = L_invV_("X,Y") * (eri3("Y, mu, nu") * C_("nu, i"));
    world.gop.fence();
    auto w1 = mpqc_time::now();
    w_times_.push_back(mpqc_time::duration_in_s(w0,w1));


    array_type J;
    J("mu, nu") = eri3("X, mu, nu")
                  * (L_invV_("Y, X") * (W("Y, rho, i") * C_("rho, i")));
    world.gop.fence();
    auto j1 = mpqc_time::now();
    j_times_.push_back(mpqc_time::duration_in_s(w1,j1));


    // Permute W
    W("X,i,nu") = W("X,nu,i");
    array_type K;
    K("mu, nu") = W("X, i, mu") * W("X, i, nu");
    world.gop.fence();
    auto k1 = mpqc_time::now();
    k_times_.push_back(mpqc_time::duration_in_s(j1,k1));

    F_("i,j") = H_("i,j") + 2 * J("i,j") - K("i,j");
}
