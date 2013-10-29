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

        for (const auto &B : ci.subspace.beta()) {
        for (auto rb : range::split(B & ci.local.beta, ci.block)) {

            auto A = ci.allowed(B); // allowed alpha subspaces
            for (range ra : A) {
                //std::cout << rb << " out of " << local << std::endl;
                mpqc::Matrix d = D(ra,rb);
                mpqc::Vector aa(ra.size());
                for (int a = 0; a < ra.size(); ++a) {
                    aa(a) = diagonal(alpha[ra[a]], h, V);
                }
#pragma omp parallel for schedule(dynamic,1)
                for (int j = 0; j < rb.size(); ++j) {
                    const String &bj = beta[rb[j]];
                    double bb = diagonal(bj, h, V);
                    for (int a = 0; a < ra.size(); ++a) {
                        double q = diagonal2(alpha[ra[a]], bj, V);
                        q = (lambda - (q + aa(a) + bb));
                        d(a,j) = (fabs(q) > 1.0e-4) ? d(a,j)/q : 0;
                    }
                }                
                D(ra,rb) = d;
                dd += dot(d, d);
            }

        }
        }

        comm.sum(dd);
        D.sync();
        D.symmetrize(1.0, 1.0/sqrt(dd));
        D.sync();

    }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_PRECONDITIONER_HPP */
