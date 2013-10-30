#ifndef MPQC_CI_DIRECT_HPP
#define MPQC_CI_DIRECT_HPP

#include <util/misc/formio.h>

#include "mpqc/ci/string.hpp"
#include "mpqc/ci/preconditioner.hpp"
#include "mpqc/ci/sigma.hpp"
#include "mpqc/ci/vector.hpp"

#include "mpqc/math/matrix.hpp"
#include "mpqc/python.hpp"
#include "mpqc/file.hpp"

namespace mpqc {
namespace ci {

    template<class Type>
    std::vector<double> direct(CI<Type> &ci,
                               const mpqc::Vector &h,
                               const mpqc::Matrix &V) {

        mpqc::Matrix lambda;
        mpqc::Vector a, r;
        size_t R = ci.roots; // roots
        size_t M = 1;

        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;
        const size_t BLOCK = alpha.size()*std::max<size_t>(1, (ci.block*ci.block)/alpha.size());

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
        const auto &local = ci.local;

        ci::Vector<Type> C("ci.C", ci);
        ci::Vector<Type> D("ci.D", ci);

        comm.barrier();

        std::vector<double> E;

        for (size_t it = 0;; ++it) {

            timer t;

            if (it + 1 > ci.max)
                throw std::runtime_error("CI failed to converge");

            // augment G
            {
                mpqc::Matrix g = G;
                G.resize(M, M);
                G.fill(0);
                range m(0, M - 1);
                G(m, m) = g(m, m);
            }

            {
                C(local.determinants).read(ci.vector.b[it]);
                C.sync();

                sigma(ci, h, V, C, D);
                D.sync();

                D(local.determinants).write(ci.vector.Hb[it]);
                comm.barrier();
            }

            // augment G matrix
            {
                mpqc::Vector g = mpqc::Vector::Zero(M);
                foreach (auto r, local.determinants.block(BLOCK)) {
                    mpqc::Vector c(r.size());
                    mpqc::Vector s = D(r);
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
                foreach (auto r, local.determinants.block(BLOCK)) {
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
                    D(r).put(d);
                }
                D.sync();

                iters[it].E = lambda[0];
                iters[it].D = norm(D, comm, local.determinants, BLOCK);

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
                    ci::Vector<Type> &b = C;
                    b(local.determinants).read(ci.vector.b[i]);
                    orthonormalize(b, D, ci.comm, local.determinants, BLOCK);
                }
                D.sync();

                D(local.determinants).write(ci.vector.b[it + 1]);

            }

            std::cout << "Davidson iteration time: " << t << std::endl;

            if (fabs(iters[it - 1].E - iters[it].E) < ci.convergence) {
                E.push_back(iters[it].E);
                break;
            }

            ++M;
        }

        return E;

    }

} // namespace ci
} // namespace mpqc

#endif /* MPQC_CI_DIRECT_HPP */
