#ifndef MPQC_CI_COLLAPSE_HPP
#define MPQC_CI_COLLAPSE_HPP

#include "mpqc/ci/string.hpp"

#include "mpqc/math.hpp"
#include "mpqc/array.hpp"
#include "mpqc/python.hpp"
#include "mpqc/file.hpp"
#include "mpqc/task.hpp"

// #define MPQC_PROFILE_ENABLE
// #include "mpqc/profile.hpp"

namespace mpqc {
namespace ci {

    /// NOT finished

    void collapse() {
        // collapse
        if (ci.collapse && (M > ci.collapse)) {
            printf("collapse %lu to %lu\n", M, ci.collapse);
            int k = 0;
            Vector c = a.col(k);

            {
                Matrix v(alpha.size(), beta.size());
                Matrix C(alpha.size(), beta.size());
                Matrix S(alpha.size(), beta.size());

                C.fill(0);
                // shift b(i+1) to b(i)
                for (int i = 0; i < M; ++i) {
                    ds.b(alpha,beta,i) >> v;
                    //ds.b(alpha,beta,j) << v;
                    C += c(i)*v;
                    printf("C* = %e*C(%i)\n", c(i), i);
                    //printf("C(%i) = C(%i)\n", j, i);
                }
                // printf("C(%lu) = C*\n", ci.collapse);
                // ds.b(alpha,beta,ci.collapse) << d;                

                S.fill(0);                
                // shift Hb(i+1) to Hb(i)
                for (int i = 0; i < M; ++i) {
                    ds.Hb(alpha,beta,i) >> v;
                    //ds.Hb(alpha,beta,j) << v;
                    S += c(i)*v;
                    printf("S* = %e*S(%i)\n", c(i), i);
                    //printf("S(%i) = S(%i)\n", j, i);
                }
                // printf("S(%lu) = S*\n", ci.collapse);
                // ds.Hb(alpha,beta,ci.collapse) << d;

                    
                for (int j = 1; j < ci.collapse; ++j) {
                    C.fill(0);
                    auto last = iters[it-j];
                    for (int i = 0; i < last.M; ++i) {
                        double a = last.a(i,k);
                        ds.b(alpha,beta,i) >> v;
                        C += a*v;
                        printf("C* = %e*C(%i)\n", a, i);
                    }
                    Matrix b(alpha.size(), beta.size());
                    ds.b(alpha,beta,i) >> b;
                    orthonormalize(alpha, beta, b, C);
                    std::cout << C << std::endl;
                }


            }

            M = ci.collapse;                

            throw;

            // orthonormalize
            for (int i = 0; i < M; ++i) {
                MPQC_PROFILE_LINE;
                Array<double> &b = C;
                b.read(ds.b[i]);
                double q = orthonormalize(alpha, beta, b, D); 
            }


            for (auto rb : range(beta).block(128)) {
                MPQC_PROFILE_LINE;
                Matrix c(alpha.size(), rb.size());
                const Matrix &s = D(alpha, rb);
                for (int j = 0; j < M; ++j) {
                    int i = it;
                    ds.b(alpha,rb,j) >> c;
                    double q = 0;
#pragma omp parallel for schedule(dynamic,1) reduction(+:q)
                    for (int b = 0; b < rb.size(); ++b) {
                        q += dot(c.col(b), s.col(b));
                    }
                    G(i,j) += q;
                    G(j,i) = G(i,j);
                }
            }
            // solve G eigenvalue
            Vector lambda = symmetric(G).eigenvalues();
            Matrix a = symmetric(G).eigenvectors();
            iters[it].lambda = lambda;
            //std::cout << "G:\n" << G << std::endl;

        }
    }

} // namespace ci
} // namespace mpqc

#endif /* MPQC_CI_COLLAPSE_HPP */
