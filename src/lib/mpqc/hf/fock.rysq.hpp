#ifndef MPQC_HF_FOCK_HPP
#define MPQC_HF_FOCK_HPP

#include "mpqc/matrix.hpp"
#include "mpqc/integrals.hpp"
#include "mpqc/omp.hpp"
#include "mpqc/utility/foreach.hpp"
#include "mpqc/utility/timer.hpp"

#include "rysq/eri.hpp"
#include <papi.h>

namespace mpqc {
namespace hf {

    typedef integrals::detail::Integral<sc::TwoBodyInt> Integral;

    inline int symmetry(int i, int j, int k, int l) {
        return 1 << ((i == j) + (k == l) + ((i == k) && (j == l)));
    }


    void fock1(Integral &int2,
               const Integral::Shell &P, const Integral::Shell &Q,
               const Integral::Shell &R, const Integral::Shell &S,
               const Matrix &D, Matrix &F, double scale) {
        const double* pqrs = int2(P,Q,R,S).data();
        foreach (int i, P) {
            foreach (int j, Q) {
                foreach (int k, R) {
                    foreach (int l, S) {
                        double v = scale*(*pqrs++);
                        F(i,j) += v*D(k,l);
                    }
                }
            }
        }
    }  

    void fock2(Integral &int2,
               const Integral::Shell &P, const Integral::Shell &Q,
               const Integral::Shell &R, const Integral::Shell &S,
               const Matrix &D, Matrix &F, double scale) {
        const double* pqrs = int2(P,Q,R,S).data();
        foreach (int i, P) {
            foreach (int j, Q) {
                foreach (int k, R) {
                    foreach (int l, S) {
                        double v = scale*(*pqrs++);
                        v *= (fabs(v) > 1e-10);
                        F(i,j) += 1.0*v*D(k,l);
                    }
                }
            }
        }
        const double* prqs = int2(P,R,Q,S).data();
        foreach (int i, P) {
            foreach (int k, R) {
                foreach (int j, Q) {
                    foreach (int l, S) {
                        double v = scale*(*prqs++);
                        v *= (fabs(v) > 1e-10);
                        F(i,j) -= 0.50*v*D(k,l);
                    }
                }
            }
        }
    }

    struct Matrix6 {
        Matrix pq, rs, pr, ps, qr, qs;
        Matrix6(size_t p, size_t q, size_t r, size_t s) {
            initialize(pq, p, q);
            initialize(rs, r, s);
            initialize(pr, p, r);
            initialize(ps, p, s);
            initialize(qr, q, r);
            initialize(qs, q, s);
        }
    private:
        void initialize(Matrix &a, size_t m, size_t n) {
            a.resize(m,n);
            a.fill(0);
        }
    };

    void fock(const Integral::Tensor4 &pqrs,
              const range &p, const range &q,
              const range &r, const range &s,
              const Matrix (&D)[6], Matrix (&F)[6],
              double scale = 1.0) {
        foreach (int i, p) {
            foreach (int j, q) {
                foreach (int k, r) {
                    foreach (int l, s) {
                        double v = scale*(pqrs(i,j,k,l));
                        //v *= (fabs(v) > 1e-15);
                        F[0](i,j) += v*D[1](k,l);
                        F[1](k,l) += v*D[0](i,j);
                        F[2](i,k) += v*D[5](j,l);
                        F[3](i,l) += v*D[4](j,k);
                        F[4](j,k) += v*D[3](i,l);
                        F[5](j,l) += v*D[2](i,k);
                    }
                }
            }
        }
    }

    void fock(Integral &int2,
              const Integral::Shell &p, const Integral::Shell &q,
              const Integral::Shell &r, const Integral::Shell &s,
              const Matrix &D, Matrix &F, double scale) {

        Matrix F6[] = {
            Matrix::Zero(p.size(), q.size()),
            Matrix::Zero(r.size(), s.size()),
            Matrix::Zero(p.size(), r.size()),
            Matrix::Zero(p.size(), s.size()),
            Matrix::Zero(q.size(), r.size()),
            Matrix::Zero(q.size(), s.size())
        };

        Matrix D6[] = {
            D(p,q), D(r,s), D(p,r), D(p,s), D(q,r), D(q,s)
        };

        fock(int2(p,q,r,s),
             range0(p), range0(q), range0(r), range0(s),
             D6, F6);

#pragma omp critical
        {
            F(p,q) += scale*1.00*F6[0];
            F(r,s) += scale*1.00*F6[1];
            F(p,r) -= scale*0.25*F6[2];
            F(p,s) -= scale*0.25*F6[3];
            F(q,r) -= scale*0.25*F6[4];
            F(q,s) -= scale*0.25*F6[5];
        }

    }

    // int max2(int d, const char *Dmax, int i, int j) {
    //     if (j < i) std::swap(i,j);
    //     return Dmax[i+(j*j+j)/2]
    // }

    void fock(const sc::GaussianBasisSet &basis,
              Integral &int2,
              const Matrix &D, Matrix &F,
              const matrix<signed char> &Dmax2,
              double tol) {

        basis.print(std::cout);

        int tol2 = (int)log2(tol);
        int N = basis.nshell();

        struct {
            std::vector< ::rysq::Shell* > basis;
            ::rysq::Eri *eri;
            double *buffer;
        } rysq;
        ::rysq::initialize();
        rysq.buffer = new double[100000];

        for (int i = 0; i < N; ++i) {
            const auto &shell = basis.shell(i);
            int nc = shell.ncontraction();
            int K = shell.nprimitive();

            shell.print(std::cout);

            rysq::type type = rysq::type(shell.max_am());
            if (nc == 2) type = rysq::type(-type);
            if (nc > 2) throw;

            double e[K];
            double C[K][2];

            for (int k = 0; k < K; ++k) {
                e[k] = shell.exponent(k);
                for (int j = 0; j < nc; ++j) {
                    C[k][j] = shell.coefficient_unnorm(j,k);
                }
                if (nc == 2) std::swap(C[k][0], C[k][1]);
            }

            int r = basis.shell_to_center(i);
            rysq::Center R = { basis.r(r,0), basis.r(r,1), basis.r(r,2) };

            rysq.basis.push_back(new rysq::Shell(type, K, e, C, R));
            std::cout << *rysq.basis.back() << std::endl;
        }

        rysq.eri = new ::rysq::Eri();

        if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) exit(1);

        int events = PAPI_NULL;
        if (PAPI_create_eventset(&events) != PAPI_OK) throw;
        if (PAPI_add_event(events, PAPI_TOT_CYC) != PAPI_OK) throw;
        if (PAPI_add_event(events, PAPI_DP_OPS) != PAPI_OK) throw;
        if (PAPI_add_event(events, PAPI_VEC_DP) != PAPI_OK) throw;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; ++j) {

                omp::task<int> task;
#pragma omp parallel 
                {
                    auto next = task++;
                    for (int k = 0, kl = 0; k <= i; ++k) {
                        for (int l = 0; l <= (k == i ? j : k); ++l, ++kl) {

                            if (next != kl) continue;
                            next = task++;

                            const auto p = Integral::Shell(i, basis.range(i));
                            const auto q = Integral::Shell(j, basis.range(j));
                            const auto r = Integral::Shell(k, basis.range(k));
                            const auto s = Integral::Shell(l, basis.range(l));

                            // {
                            //     using std::max;
                            //     int g = int2.max2(p,q,r,s);
                            //     int dmax;
                            //     dmax = max(Dmax2(i,j), Dmax2(k,l));
                            //     dmax = max(dmax, Dmax2(i,k)-1);
                            //     dmax = max(dmax, Dmax2(i,l)-1);
                            //     dmax = max(dmax, Dmax2(j,k)-1);
                            //     dmax = max(dmax, Dmax2(j,l)-1);
                            //     if (g+dmax < tol2) {
                            //         continue;
                            //     }
                            // }

                            //fock(int2, p, q, r, s, D, F, 4.0/symmetry(i,j,k,l));


                            {
                                mpqc::timer timer;
                                PAPI_start(events);

                                auto p = rysq.basis[i];
                                auto q = rysq.basis[j];
                                auto r = rysq.basis[k];
                                auto s = rysq.basis[l];

                                if (p->L < q->L) std::swap(p,q);
                                if (r->L < s->L) std::swap(r,s);

                                if ((p->L + q->L) < (r->L + s->L)) {
                                    std::swap(p,r);
                                    std::swap(q,s);
                                }

                                std::cout << "quartet "
                                          << ::rysq::shell(*p)
                                          << ::rysq::shell(*q)
                                          << ::rysq::shell(*r)
                                          << ::rysq::shell(*s)
                                          << std::endl;

                                double *eri2 = rysq.buffer;
                                (*rysq.eri)(*p, *q, *r, *s, eri2, 0);

                                long_long counter[3];
                                PAPI_stop(events, counter);

                                std::cout << " librysq:";
                                std::cout << " cycles = " << counter[0];
                                std::cout << " dp ops = " << counter[1];
                                std::cout << " vec dp = " << counter[2];
                                std::cout << " time = " << timer;
                                std::cout << std::endl;
                            }

                            {
                                mpqc::timer timer;
                                PAPI_start(events);

                                const double *eri1 = int2(p, q, r, s).data();

                                long_long counter[3];
                                PAPI_stop(events, counter);

                                std::cout << " libint:";
                                std::cout << " cycles = " << counter[0];
                                std::cout << " dp ops = " << counter[1];
                                std::cout << " vec dp = " << counter[2];
                                std::cout << " time = " << timer;
                                std::cout << std::endl;
                            }

                            std::cout << std::endl;

                            // if (p.size()*q.size()*r.size()*s.size() > 12) continue;

                            // foreach (int i, p) {
                            //     foreach (int j, q) {
                            //         foreach (int k, r) {
                            //             foreach (int l, s) {
                            //                 printf("2-e %i %i %i %i: %e %e\n",
                            //                        i, j, k, l, *eri1++, *eri2++);
                            //             }
                            //         }
                            //     }
                            // }

                        }
                    }
                }
            }
        }

        //throw;
        
        delete rysq.eri;
        delete[] rysq.buffer;

        F += Matrix(F.transpose());
        F /= 2;

    }


    void fock() {
    }


} // namespace hf
} // namespace mpqc

#endif /* MPQC_HF_FOCK_HPP */
