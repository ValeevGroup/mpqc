//
// Created by Chong Peng on 2/9/16.
//

#ifndef MPQC_SCF_SCF_H
#define MPQC_SCF_SCF_H


#include "../integrals/integrals.h"
#include "../include/tiledarray.h"
#include "../utility/time.h"
#include "../ta_routines/array_to_eigen.h"
#include "diagonalize_for_coffs.hpp"

namespace mpqc {
namespace scf {

template <typename FockBuilder>
class ClosedShellSCF {
  public:
    using array_type = TA::TSpArrayD;

  protected:
    array_type H_;
    array_type S_;
    array_type F_;
    array_type D_;
    array_type C_;
    TiledArray::DIIS<array_type> diis_;

    FockBuilder builder_;

    std::vector<double> scf_times_;

    int64_t occ_;
    double repulsion_;

  public:
    ClosedShellSCF(array_type const &H, array_type const &S, int64_t occ,
                   double rep, FockBuilder builder,
                   array_type const &F_guess = array_type{})
            : H_(H),
              S_(S),
              builder_(std::move(builder)),
              occ_(occ),
              repulsion_(rep) {

        if (F_guess.is_initialized()) {
            F_ = F_guess;
        } else {
            F_ = H_;
        }

        compute_density(occ_);
    }

    array_type const &overlap() const { return S_; }
    array_type const &fock() const { return F_; }
    array_type const &density() const { return D_; }
    array_type const &coefficents() const { return C_; }

    double energy() {
        return repulsion_
               + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
    }


    /*! Function to compute the density to the desired accuracy.
     *
     * Takes some form of integral and does the scf iterations.  The place to
     *specialized is in build_fock.
     *
     * returns true if the calculation converged to the desired threshold in
     *fewer than max_iters
     */
    template <typename... FockBuilderArgs>
    bool solve(int64_t max_iters, double thresh, FockBuilderArgs... args) {
        auto &world = F_.get_world();

        if (world.rank() == 0) {
            std::cout << "Starting SCF:\n"
                      << "\tThreshold: " << thresh << "\n"
                      << "\tMaximum number of iterations: " << max_iters
                      << "\n";
        }

        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto rms_error = std::numeric_limits<double>::max();
        auto old_energy = energy();
        const double volume = F_.trange().elements().volume();

        while (iter < max_iters && (thresh < (error / old_energy)
                                    || thresh < (rms_error / volume))) {

            auto s0 = mpqc_time::fenced_now(world);

            F_("i,j") = H_("i,j")
                        + builder_(std::forward<FockBuilderArgs>(args)...,
                                   C_)("i,j");

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
                          << (error / old_energy) << "\n"
                          << "\t\t(Gradient Norm)/n^2: " << (rms_error / volume)
                          << "\n"
                          << "\t\tTime: " << scf_times_.back() << "\n";
            }
            ++iter;
        }

        if (iter < max_iters && thresh > (error / old_energy)
            && thresh > (rms_error / volume)) {

            if(world.rank() == 0){
                std::cout << "HF Energy:  " << old_energy << std::endl;
            }

            return true;
        } else {
            return false;
        }
    }


  private:
    void compute_density(int64_t occ) {
        auto F_eig = array_ops::array_to_eigen(F_);
        auto S_eig = array_ops::array_to_eigen(S_);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        decltype(S_eig) C = es.eigenvectors().leftCols(occ);
        auto tr_ao = S_.trange().data()[0];

        auto occ_nclusters = (occ_ < 10) ? occ_ : 10;
        auto tr_occ = scf::tr_occupied(occ_nclusters, occ_);

        C_ = array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C, tr_ao,
                                                    tr_occ);

        D_("i,j") = C_("i,k") * C_("j,k");
    }
};

} // end of namespace scf
} // end of namespace mpqc


#endif // MPQC_SCF_SCF_H
