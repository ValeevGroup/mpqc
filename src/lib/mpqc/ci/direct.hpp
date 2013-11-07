#ifndef MPQC_CI_DIRECT_HPP
#define MPQC_CI_DIRECT_HPP

#include <util/misc/formio.h>

#include "mpqc/ci/string.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/ci/sigma.hpp"
#include "mpqc/ci/preconditioner.hpp"

#include "mpqc/math/matrix.hpp"
#include "mpqc/file.hpp"

#include "mpqc/utility/profile.hpp"

//#define MPQC_CI_VERBOSE 1

namespace mpqc {
namespace ci {

    // read local segments into V from F
    inline void read(ci::Vector &V, File::Dataspace<double> F,
                     const std::vector<mpqc::range> &local) {
        timer t;
        size_t count = 0;
        foreach (auto r, local) {
            mpqc::Vector v(r.size());
            F(r) >> v;
            V(r) << v;
            count += r.size();
        }
#if MPQC_CI_VERBOSE
        printf("read took %f s, %f mb/s\n", (double)t, count*sizeof(double)/(t*(1<<20)));
#endif
    }

    // write local segments of V to F
    inline void write(ci::Vector &V, File::Dataspace<double> F,
                      const std::vector<mpqc::range> &local) {
        timer t;
        size_t count = 0;
        foreach (auto r, local) {
            mpqc::Vector v(r.size());
            V(r) >> v;
            F(r) << v;
            count += r.size();
        }
#if MPQC_CI_VERBOSE
        printf("write took %f s, %f mb/s\n", (double)t, count*sizeof(double)/(t*(1<<20)));
#endif
    }


    template<class Type>
    std::vector<double> direct(CI<Type> &ci,
                               const mpqc::Vector &h,
                               const mpqc::Matrix &V) {

        MPQC_PROFILE_REGISTER_THREAD;

        mpqc::Matrix lambda;
        mpqc::Vector a, r;
        size_t R = ci.roots; // roots
        size_t M = 1;

        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;

        mpqc::Matrix G;
        struct Iter {
            double E, D;
            size_t M;
            mpqc::Matrix G;
            mpqc::Vector lambda;
            mpqc::Matrix a;
        };
        std::map<int, Iter> iters;
        iters[-1].E = iters[-1].D = 0;

        auto &comm = ci.comm;

        ci::Vector C("ci.C", ci.subspace, comm, (ci.incore >= 1));
        ci::Vector D("ci.D", ci.subspace, comm, (ci.incore >= 2));

        comm.barrier();

        std::vector<double> E;

        for (size_t it = 0;; ++it) {

            timer t;

            // augment G
            {
                mpqc::Matrix g = G;
                G.resize(M, M);
                G.fill(0);
                range m(0, M - 1);
                G(m, m) = g(m, m);
            }

            {
                // read local segments of b(it) to C
                read(C, ci.vector.b[it], ci.local());
                C.sync();

                sigma(ci, h, V, C, D);
                D.sync();

                // write local segments of D to Hb(it)
                write(D, ci.vector.Hb[it], ci.local());
                D.sync();
            }

            // augment G matrix
            {
                MPQC_PROFILE_LINE;
                mpqc::Vector g = mpqc::Vector::Zero(M);
                foreach (auto r, ci.local()) {
                    mpqc::Vector c(r.size());
                    mpqc::Vector s(r.size());
                    D(r) >> s;
                    for (int j = 0; j < M; ++j) {
                        ci.vector.b(r,j) >> c;
                        g(j) += c.dot(s);
                    }
                }
                comm.sum(g.data(), M);
                G.row(it) = g;
                G.col(it) = g;
                //std::cout << "G = \n" << G << std::endl;
            }

            // solve G eigenvalue
            mpqc::Vector lambda = symmetric(G).eigenvalues();
            mpqc::Matrix a = symmetric(G).eigenvectors();

            iters[it].M = M;
            iters[it].G = G;
            iters[it].lambda = lambda;
            iters[it].a = a;

            // std::cout << "lambda:\n" << lambda << std::endl;
            // std::cout << "alpha:\n" << a << std::endl;

            for (int k = 0; k < R; ++k) {

                // update d part
                for (auto r : ci.local()) {
                    MPQC_PROFILE_LINE;
                    mpqc::Vector v(r.size());
                    mpqc::Vector d(r.size());
                    d.fill(0);
                    for (int i = 0; i < M; ++i) {
                        ci.vector.b(r,i) >> v;
                        double r = -a(i,k)*lambda(k);
                        d += r*v;
                    }
                    for (int i = 0; i < M; ++i) {
                        ci.vector.Hb(r,i) >> v;
                        d += a(i,k)*v;
                    }
                    D(r) << d;
                }
                D.sync();

                iters[it].E = lambda[0];
                iters[it].D = norm(D, ci.local(), comm);

                if (comm.rank() == 0) {
                    double dc = fabs(iters[it - 1].D - iters[it].D);
                    double de = fabs(iters[it - 1].E - iters[it].E);
                    sc::ExEnv::out0()
                        << sc::indent
                        << sc::scprintf("CI iter. %3i, E=%15.12lf, "
                                        "del.E=%4.2e, del.C=%4.2e\n",
                                        (int)it, lambda[0] + ci.e_ref + ci.e_core,
                                        de, dc);
                }
                
                preconditioner(ci, h, V, lambda[0], D);

                // orthonormalize
                for (int i = 0; i < M; ++i) {
                    ci::Vector &b = C;
                    read(b, ci.vector.b[i], ci.local());
                    orthonormalize(b, D, ci.local(), ci.comm);
                }
                D.sync();

                if (it+1 == ci.max) break;

                write(D, ci.vector.b[it+1], ci.local());
                D.sync();

            }

            MPQC_PROFILE_DUMP(std::cout);

            sc::ExEnv::out0() << sc::indent << "Davidson iteration time: " << t << std::endl;

            if (fabs(iters[it - 1].E - iters[it].E) < ci.convergence) {
                E.push_back(iters[it].E);
                break;
            }

            if (it+1 == ci.max) {
                throw MPQC_EXCEPTION("CI failed to converge");
            }

            ++M;
        }

        return E;

    }

} // namespace ci
} // namespace mpqc

#endif /* MPQC_CI_DIRECT_HPP */
