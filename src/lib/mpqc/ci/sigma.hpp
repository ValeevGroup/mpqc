#ifndef MPQC_CI_SIGMA_HPP
#define MPQC_CI_SIGMA_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"
#include "mpqc/ci/sigma2.hpp"
#include "mpqc/ci/sigma3.hpp"
#include "mpqc/ci/vector.hpp"

#include "mpqc/utility/timer.hpp"
#include "mpqc/range.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/mpi.hpp"
#include "mpqc/mpi/task.hpp"

#include "mpqc/utility/profile.hpp"

namespace mpqc {
namespace ci {

    /// Computes sigma 1,2,3 contributions
    template<class Type, class Index>
    void sigma(const CI<Type, Index> &ci,
               const mpqc::Vector &h, const Matrix &V,
               ci::Vector &C, ci::Vector &S) {

        struct { double s1, s2, s3; timer t; } time = { };

        auto &comm = ci.comm;

        size_t no = ci.alpha[0].size();
        mpqc::Vector H = h;
        for (int j = 0; j < no; ++j) {
            for (int i = 0; i <= j; ++i) {
                double v = 0;
                for (int k = 0; k < no; ++k) {
                    v += V(index(i,k), index(k,j));
                }
                H(index(i,j)) -= 0.5*v;
            }
        }

        const std::vector< Subspace<Alpha> > &alpha = ci.subspace.alpha();
        const std::vector< Subspace<Beta> > &beta = ci.subspace.beta();
        const auto &blocks = ci::blocks(alpha, beta);

        std::unique_ptr<MPI::Task> task;

        task.reset(new MPI::Task(comm));
#pragma omp parallel
        while (true) {

            auto next = task->next(blocks.begin(), blocks.end());
            if (next == blocks.end()) break;

            auto Ia = alpha.at(next->alpha);
            auto Ib = beta.at(next->beta);
            if (!ci.test(Ia,Ib)) continue;

            Matrix s = Matrix::Zero(Ia.size(), Ib.size()); //S(Ia, Ib);

            // sigma1
            foreach (auto Jb, beta) {
                MPQC_PROFILE_LINE;
                // only single and double excitations are allowed
                if (!ci.test(Ia,Jb) || ci.diff(Ib,Jb) > 2) continue;
                Matrix c = C(Ia,Jb);
                timer t;
                sigma12(ci, Ib, Jb, H, V, c, s);
#pragma omp master
                time.s1 += t;
            }

            // with ms == 0 symmetry, S is symmetrized in sigma3 step
            if (ci.ms == 0) goto end;

            // sigma2, need to transpose s, c
            s = Matrix(s.transpose());
            foreach (auto Ja, alpha) {
                MPQC_PROFILE_LINE;
                if (!ci.test(Ja,Ib) || ci.diff(Ia,Ja) > 2) continue;
                Matrix c = Matrix(C(Ja,Ib)).transpose();
                timer t;
                sigma12(ci, Ia, Ja, H, V, c, s);
                time.s2 += t;
            }
            s = Matrix(s.transpose());

        end:
            S(Ia,Ib) = s;

        }

        task.reset(new MPI::Task(comm));
#pragma omp parallel
        while (true) {

            timer t;

            auto next = task->next(blocks.begin(), blocks.end());
            if (next == blocks.end()) break;

            if (ci.ms == 0 && next->alpha > next->beta) continue;

            auto Ia = alpha.at(next->alpha);
            auto Ib = beta.at(next->beta);
            if (!ci.test(Ia,Ib)) continue;

            // excitations from Ib into each Jb subspace
            std::vector< Excitations<Beta> > BB;
            foreach (auto Jb, beta) {
                MPQC_PROFILE_LINE;
                BB.push_back(Excitations<Beta>(ci, Ib, Jb));
            }

            // excitations from Ia into each Ja subspace
            std::vector< Excitations<Alpha> > AA;
            foreach (auto Ja, alpha) {
                MPQC_PROFILE_LINE;
                AA.push_back(Excitations<Alpha>(ci, Ia, Ja));
            }
            
            Matrix s = S(Ia,Ib);

            // if symmetric CI, symmetrize diagonal block
            if (ci.ms == 0 && next->alpha == next->beta) 
                s += Matrix(s).transpose();

            for (auto bb = BB.begin(); bb != BB.end(); ++bb) {
                MPQC_PROFILE_LINE;
                if (!bb->size()) continue; // no beta excitations
                auto Jb = bb->J();
                for (auto aa = AA.begin(); aa != AA.end(); ++aa) {
                    if (!aa->size()) continue; // no alpha excitations
                    auto Ja = aa->J();
                    if (!ci.test(Ja,Jb)) continue; // forbidden block
                    Matrix c = C(Ja,Jb);
                    timer t;
                    sigma3(*aa, *bb, V, c, s);
#pragma omp master
                    time.s3 += t;
                }
            }

            // if symmetric CI, symmetrize off-diagonal blocks S(Ia,Ib) and S(Ib,Ia)
            if (ci.ms == 0 && next->alpha != next->beta) {
                MPQC_PROFILE_LINE;
                Matrix t = S(Ib,Ia);
                t += s.transpose();
                s = t.transpose();
                S(Ib,Ia) = t;
            }

            {
                MPQC_PROFILE_LINE;
                S(Ia,Ib) = s;
            }
            
        }

        S.sync();

        sc::ExEnv::out0() << sc::indent << "sigma took " << double(time.t) << std::endl;
        sc::ExEnv::out0() << sc::indent << "  sigma1: " << time.s1 << std::endl;
        sc::ExEnv::out0() << sc::indent << "  sigma2: " << time.s2 << std::endl;
        sc::ExEnv::out0() << sc::indent << "  sigma3: " << time.s3 << std::endl;


    }

}
}

#endif /* MPQC_CI_SIGMA_HPP */

