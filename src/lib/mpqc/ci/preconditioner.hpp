#ifndef MPQC_CI_PRECONDITIONER_HPP
#define MPQC_CI_PRECONDITIONER_HPP

#include "mpqc/ci/hamiltonian.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/utility/foreach.hpp"

namespace mpqc {
namespace ci {
                
    template<class Type>
    void preconditioner(CI<Type> &ci,
                        const mpqc::Vector &h, const mpqc::Matrix &V,
                        double lambda, mpqc::Array<double> &D) {
        
        const auto &comm = ci.comm;
        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;

        double dd = 0;
        foreach (auto rb, range(ci.local.beta).block(ci.block)) {
            //std::cout << rb << " out of " << local << std::endl;
            mpqc::Matrix d = D(alpha,rb);

            mpqc::Vector aa(alpha.size());
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

            D(alpha,rb) << d;
            dd += dot(d, d);
        }
        comm.sum(dd);
        D.sync();
        symmetrize(D, 1, 1/sqrt(dd), comm);
        D.sync();

    }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_PRECONDITIONER_HPP */
