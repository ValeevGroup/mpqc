#ifndef MPQC_CI_PRECONDITIONER_HPP
#define MPQC_CI_PRECONDITIONER_HPP

#include "mpqc/ci/hamiltonian.hpp"
#include "mpqc/ci/array.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/utility/foreach.hpp"

namespace mpqc {
namespace ci {
                
    template<class Type>
    void preconditioner(CI<Type> &ci,
                        const Vector &h, const Matrix &V,
                        double lambda, ci::Array &D) {
        
        const auto &comm = ci.comm;
        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;

        double dd = 0;
        foreach (auto rb, range(beta).block(1)) {
            Matrix d = D.array(alpha,rb);

            Vector aa(alpha.size());
            for (int a = 0; a < alpha.size(); ++a) {
                aa(a) = diagonal(alpha[a], h, V);
            }

#pragma omp parallel for schedule(dynamic,1)
            for (int j = 0; j < rb.size(); ++j) {
                const String &bj = beta[rb[j]];
                double bb = diagonal(bj, h, V);
                for (int a = 0; a < alpha.size(); ++a) {
                    double q = diagonal2(alpha[a], bj, V);
                    q = (lambda - (q + aa(a) + bb));
                    d(a,j) = (fabs(q) > 1.0e-4) ? d(a,j)/q : 0;
                }
            }

            D.array(alpha,rb) << d;
            dd += dot(d, d);
        }
        comm.sum(dd);
        D.sync();
        if (comm.rank() == 0)
            symmetrize(D.array(), 1, 1/sqrt(dd));
        D.sync();
    }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_PRECONDITIONER_HPP */
