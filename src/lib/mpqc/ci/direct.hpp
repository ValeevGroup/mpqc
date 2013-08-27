#ifndef MPQC_CI_DIRECT_HPP
#define MPQC_CI_DIRECT_HPP

#include <util/misc/formio.h>

//#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"

#include "mpqc/ci/string.hpp"
#include "mpqc/ci/hamiltonian.hpp"
#include "mpqc/ci/sigma.hpp"
#include "mpqc/ci/array.hpp"

#include "mpqc/math/matrix.hpp"
#include "mpqc/python.hpp"
#include "mpqc/file.hpp"

namespace mpqc {
namespace ci {

    template<class Type>
    std::vector<double> direct(CI<Type> &ci, const Vector &h, const Matrix &V) {

        MPQC_PROFILE_REGISTER_THREAD;

        Matrix lambda;
        Vector a, r;
        size_t MAX = ci.max;
        size_t R = ci.roots; // roots
        size_t M = 1;

        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;

        Matrix G;
        struct Iter {
            double E, D;
            size_t M;
            Matrix G;
            Vector lambda;
            Matrix a;
        };
        std::map<int, Iter> iters;
        iters[-1].E = iters[-1].D = 0;

        auto &comm = ci.comm;
        auto &io = ci.io;
        range local = ci.io.local;

        ci::Array C("ci.C", alpha.size(), beta.size(), comm);
        ci::Array D("ci.D", alpha.size(), beta.size(), comm);

        // mpqc::Array<double> C("ci.C", ci.dims, ARRAY_FILE, comm);
        // mpqc::Array<double> D("ci.D", ci.dims, ARRAY_FILE, comm);

        comm.barrier();

        std::vector<double> E;

        for (size_t it = 0;; ++it) {

            timer t;

            if (it + 1 > MAX)
                throw std::runtime_error("CI failed to converge");

            // augment G
            {
                Matrix g = G;
                G.resize(M, M);
                G.fill(0);
                range m(0, M - 1);
                G(m, m) = g(m, m);
            }

            {
                C.vector(local).read(io.b[it]);
                C.sync();

                sigma(ci, h, V, C.array(), D.array());
                D.sync();

                D.vector(local).write(io.Hb[it]);
                comm.barrier();
            }

            MPQC_PROFILE_RESET;

            // augment G matrix
            {
                Vector g = Vector::Zero(M);
                foreach (auto r, local.block(alpha.size())) {
                    Vector c(r.size());
                    const Vector &s = D.vector(r);
                    for (int j = 0; j < M; ++j) {
                        io.b(r,j) >> c;
                        g(j) += c.dot(s);
                    }
                }
                //std::cout << "g: " << g << std::endl;
                comm.sum(g.data(), M);
                G.row(it) = g;
                G.col(it) = g;
            }

            // solve G eigenvalue
            Vector lambda = symmetric(G).eigenvalues();
            Matrix a = symmetric(G).eigenvectors();

            iters[it].M = M;
            iters[it].G = G;
            iters[it].lambda = lambda;
            iters[it].a = a;

            // std::cout << "lambda:\n" << lambda << std::endl;
            // std::cout << "alpha:\n" << a << std::endl;

            for (int k = 0; k < R; ++k) {

                //                 MPQC_PROFILE_LINE;

                // update d part
                foreach (auto r, local.block(alpha.size())) {
                    Vector v(r.size());
                    Vector d(r.size());
                    d.fill(0);
                    for (int i = 0; i < M; ++i) {
                        MPQC_PROFILE_LINE;
                        io.b(r,i) >> v;
                        double r = -a(i,k)*lambda(k);
                        d += r*v;
                    }
                    for (int i = 0; i < M; ++i) {
                        io.Hb(r,i) >> v;
                        d += a(i,k)*v;
                    }
                    D.vector(r) << d;
                }
                D.sync();

                iters[it].E = lambda[0];
                iters[it].D = norm(D, comm, local);

                if (comm.rank() == 0) {
                    double dc = fabs(iters[it - 1].D - iters[it].D);
                    double de = fabs(iters[it - 1].E - iters[it].E);
                    sc::ExEnv::out0()
                        << sc::indent
                        << sc::scprintf("CI iter. %3i, E=%15.12lf, del.E=%4.2e, del.C=%4.2e\n",
                                        (int) it, lambda[0] + ci.e_ref + ci.e_core, de, dc);
                }

                // preconditioner
                {
                    MPQC_PROFILE_LINE;
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
                                q = (lambda[0] - (q + aa(a) + bb));
                                d(a,j) = (fabs(q) > 1.0e-4) ? d(a,j)/q : 0;
                            }
                        }

                        D.array(alpha,rb) << d;
                        dd += dot(d, d);
                    }
                    comm.sum(dd);
                    D.sync();
                    if (comm.rank() == 0)
                        symmetrize(D.array(), 1, 1 / sqrt(dd));
                    D.sync();
                }

                // orthonormalize
                for (int i = 0; i < M; ++i) {
                    MPQC_PROFILE_LINE;
                    ci::Array &b = C;
                    b.vector(local).read(io.b[i]);
                    orthonormalize(b.vector(local), D.vector(local),
                                   ci.comm, alpha.size());
                }
                D.sync();

                D.vector(local).write(io.b[it + 1]);

            }

            MPQC_PROFILE_DUMP(std::cout);

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
