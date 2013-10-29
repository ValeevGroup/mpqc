#ifndef MPQC_CI_SIGMA_HPP
#define MPQC_CI_SIGMA_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"
#include "mpqc/ci/sigma2.hpp"
#include "mpqc/ci/sigma3.hpp"

#include "mpqc/utility/timer.hpp"
#include "mpqc/range.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/omp.hpp"

#include "mpqc/array.hpp"
#include "mpqc/array/functions.hpp"

namespace mpqc {
namespace ci {

    /// Computes sigma 1,2,3 contributions
    template<class Type, class Index>
    void sigma(const CI<Type, Index> &ci,
               const mpqc::Vector &h, const Matrix &V,
               ci::Vector<Type> &C, ci::Vector<Type> &S) {

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

        const int ms = 0;

        std::vector< Subspace<Alpha> > alpha = split(ci.subspace.alpha(), ci.block);
        std::vector< Subspace<Beta> > beta = split(ci.subspace.beta(), ci.block);

        struct Tuple {
            int a, b;
            Tuple(int a, int b) : a(a), b(b) {}
        };
        std::vector<Tuple> tuples;
        // (a,b) block tuples
        for (int b = 0; b < beta.size(); ++b) {
            for (int a = 0; a < alpha.size(); ++a) {
                if (!ci.test(alpha.at(a), beta.at(b))) continue;
                tuples.push_back(Tuple(a,b));
            }
        }

        std::auto_ptr<MPI::Task> task;

        task.reset(new MPI::Task(comm));
#pragma omp parallel
        while (true) {

            auto next = task->next(tuples.begin(), tuples.end());
            if (next == tuples.end()) break;
            auto Ia = alpha.at(next->a);
            auto Ib = beta.at(next->b);

            Matrix s = Matrix::Zero(Ia.size(), Ib.size()); //S(Ia, Ib);

            // sigma1
            foreach (auto Jb, beta) {
                // only single and double excitations are allowed
                if (!ci.test(Ia,Jb) || ci.diff(Ib,Jb) > 2) continue;
                Matrix c = C(Ia,Jb);
                timer t;
                sigma12(ci, Ib, Jb, H, V, c, s);
#pragma omp master
                time.s1 += t;
            }

            // with ms == 0 symmetry, S is symmetrized in sigma3 step
            if (ms == 0) goto end;

            // sigma2, need to transpose s, c
            s = Matrix(s.transpose());
            foreach (auto Ja, alpha) {
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

            auto next = task->next(tuples.begin(), tuples.end());
            if (next == tuples.end()) break;

            if (ms == 0 && next->a > next->b) continue;

            auto Ia = alpha.at(next->a);
            auto Ib = beta.at(next->b);
            if (!ci.test(Ia,Ib)) continue;

            // excitations from Ib into each Jb subspace
            std::vector< Excitations<Beta> > BB;
            foreach (auto Jb, beta) {
                BB.push_back(Excitations<Beta>(ci, Ib, Jb));
            }

            // excitations from Ia into each Ja subspace
            std::vector< Excitations<Alpha> > AA;
            foreach (auto Ja, alpha) {
                AA.push_back(Excitations<Alpha>(ci, Ia, Ja));
            }
            
            Matrix s = S(Ia,Ib);

            // if symmetric CI, symmetrize diagonal block
            if (ms == 0 && next->a == next->b) 
                s += Matrix(s).transpose();

            for (auto bb = BB.begin(); bb != BB.end(); ++bb) {
                //if (!bb->size()) continue; // no beta excitations
                auto Jb = bb->J();
                for (auto aa = AA.begin(); aa != AA.end(); ++aa) {
                    //if (!aa->size()) continue; // no alpha excitations
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
            if (ms == 0 && next->a != next->b) {
                Matrix t = S(Ib,Ia);
                t += s.transpose();
                s = t.transpose();
                S(Ib,Ia) = t;
            }
            
            S(Ia,Ib) = s;

        }

        S.sync();

        std::cout << "sigma took " << double(time.t) << std::endl;
        std::cout << "  sigma1: " << time.s1 << std::endl;
        std::cout << "  sigma2: " << time.s2 << std::endl;
        std::cout << "  sigma3: " << time.s3 << std::endl;


    }

}
}

#endif /* MPQC_CI_SIGMA_HPP */

