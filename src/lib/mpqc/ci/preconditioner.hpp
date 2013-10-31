#ifndef MPQC_CI_PRECONDITIONER_HPP
#define MPQC_CI_PRECONDITIONER_HPP

#include "mpqc/ci/hamiltonian.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/utility/foreach.hpp"
#include "mpqc/mpi.hpp"
#include "mpqc/mpi/task.hpp"

namespace mpqc {
namespace ci {

    void symmetrize(mpqc::Matrix &a, double phase, double scale) {
        MPQC_CHECK(a.rows() == a.cols());
        for (size_t j = 0; j < a.cols(); ++j) {
            a(j,j) *= scale;
            for (size_t i = 0; i < j; ++i) {
                a(i,j) = scale*a(i,j); // + a(j,i));
                a(j,i) = phase*a(i,j);
            }
        }
    }

    template<class Type, class Index>
    void preconditioner(CI<Type, Index> &ci,
                        const mpqc::Vector &h, const mpqc::Matrix &V,
                        double lambda, ci::BlockVector &D) {

        const auto &comm = ci.comm;
        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;
        const int ms = 0;

        double dd = 0;

        const std::vector< Subspace<Alpha> > &A = ci.subspace.alpha();
        const std::vector< Subspace<Beta> > &B = ci.subspace.beta();
        const auto &blocks = ci::blocks(A, B);

        std::unique_ptr<MPI::Task> task;

        task.reset(new MPI::Task(comm));
#pragma omp parallel
        while (true) {

            auto next = task->next(blocks.begin(), blocks.end());
            if (next == blocks.end()) break;

            auto Ia = A.at(next->alpha);
            auto Ib = B.at(next->beta);
            if (!ci.test(Ia,Ib)) continue;

            mpqc::Matrix d = D(Ia,Ib);
            mpqc::Vector aa(Ia.size());
            for (int a = 0; a < Ia.size(); ++a) {
                aa(a) = diagonal(alpha[*Ia.begin() + a], h, V);
            }
            for (int j = 0; j < Ib.size(); ++j) {
                const String &bj = beta[*Ib.begin() + j];
                double bb = diagonal(bj, h, V);
                for (int a = 0; a < Ia.size(); ++a) {
                    double q = diagonal2(alpha[*Ia.begin() + a], bj, V);
                    q = (lambda - (q + aa(a) + bb));
                    d(a,j) = (fabs(q) > 1.0e-4) ? d(a,j)/q : 0;
                }
            } 
            D(Ia,Ib) = d;
            dd += dot(d, d);
        }
        comm.sum(dd);
        D.sync();

        if (ms != 0) return;
 
        double phase = 1.0;
        double scale = 1.0/sqrt(dd);

        task.reset(new MPI::Task(comm));
#pragma omp parallel
        while (true) {

            auto next = task->next(blocks.begin(), blocks.end());
            if (next == blocks.end()) break;

            auto Ia = A.at(next->alpha);
            auto Ib = B.at(next->beta);
            if (!ci.test(Ia,Ib)) continue;

            if (next->alpha == next->beta) {
                Matrix aa = D(Ia,Ib);
                ci::symmetrize(aa, phase, scale);
                D(Ia,Ib) = aa;
            }

            if (next->alpha < next->beta) {
                Matrix ab = D(Ia,Ib);
                Matrix ba = D(Ib,Ia);
                ab *= scale;
                ba  = phase*(ab.transpose());
                D(Ia,Ib) = ab;
                D(Ib,Ia) = ba;
            }

        }
        D.sync();

    }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_PRECONDITIONER_HPP */
