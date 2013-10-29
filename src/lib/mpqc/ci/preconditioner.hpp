#ifndef MPQC_CI_PRECONDITIONER_HPP
#define MPQC_CI_PRECONDITIONER_HPP

#include "mpqc/ci/hamiltonian.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/utility/foreach.hpp"
#include "mpqc/utility/exception.hpp"

namespace mpqc {
namespace ci {
    
    template<class Type, class Index>
    void preconditioner(CI<Type, Index> &ci,
                        const mpqc::Vector &h, const mpqc::Matrix &V,
                        double lambda, ci::Vector<Type> &D) {

        const auto &comm = ci.comm;
        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;

        double dd = 0;

        std::vector< Subspace<Alpha> > A = split(ci.subspace.alpha(), ci.block);
        std::vector< Subspace<Beta> > B = split(ci.subspace.beta(), ci.block);

        struct Tuple {
            int a, b;
            Tuple(int a, int b) : a(a), b(b) {}
        };
        std::vector<Tuple> tuples;
        // (a,b) block tuples
        for (int b = 0; b < B.size(); ++b) {
            for (int a = 0; a < A.size(); ++a) {
                if (!ci.test(A.at(a), B.at(b))) continue;
                tuples.push_back(Tuple(a,b));
            }
        }

        MPI::Task task(comm);
#pragma omp parallel
        while (true) {

            auto next = task.next(tuples.begin(), tuples.end());
            if (next == tuples.end()) break;
            auto Ia = A.at(next->a);
            auto Ib = B.at(next->b);

            //std::cout << rb << " out of " << local << std::endl;
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
        D.symmetrize(1.0, 1.0/sqrt(dd));
        D.sync();

    }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_PRECONDITIONER_HPP */
