#ifndef MPQC_CI_DIRECT_HPP
#define MPQC_CI_DIRECT_HPP

#include <util/misc/formio.h>

//#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"

#include "mpqc/ci/string.hpp"
#include "mpqc/ci/hamiltonian.hpp"
#include "mpqc/ci/sigma.hpp"

#include "mpqc/matrix.hpp"
#include "mpqc/array.hpp"
#include "mpqc/python.hpp"
#include "mpqc/file.hpp"

namespace mpqc {
namespace ci {

    double norm(const mpqc::Array<double> &A, const MPI::Comm &comm) {
        MPQC_PROFILE_LINE;
        range r1 = range(0, A.dims()[0]);
        auto r2 = range(0, A.dims()[1]).block(128);
        double n = 0;
        MPI::Task task(comm);
        int j;
        while ((j = task++) < r2.size()) {
            n += norm(Matrix(A(r1, r2[j])));
        }
        comm.sum(n);
        return n;
    }

    void symmetrize(Matrix &a, double phase, double scale) {
        for (size_t j = 0; j < a.cols(); ++j) {
            a(j, j) *= scale;
            for (size_t i = 0; i < j; ++i) {
                a(i, j) = scale * a(i, j); // + a(j,i));
                a(j, i) = phase * a(i, j);
            }
        }
    }

    void symmetrize(Array<double> &A, double phase, double scale) {
        size_t N = 512;
        Matrix a;
        std::vector<range> r = range::block(range(0, A.dims()[1]), N);
        for (auto rj = r.begin(); rj < r.end(); ++rj) {
            //std::cout << *rj << std::endl;
            size_t nj = rj->size();
            for (auto ri = r.begin(); ri < rj; ++ri) {
                size_t ni = ri->size();
                a.resize(ni, nj);
                A(*ri, *rj) >> a;
                a *= scale;
                A(*ri, *rj) << a;
                a *= phase;
                A(*rj, *ri) << Matrix(a.transpose());
            }
            a.resize(nj, nj);
            A(*rj, *rj) >> a;
            symmetrize(a, phase, scale);
            A(*rj, *rj) << a;
        }
    }

    /**
     //     double db = dot(D, b);
     //     D = D - db*b;
     //     D *= 1/D.norm();
     //     @return <d,b>/||d||
     */
    double orthonormalize(range alpha, range beta, const Array<double> &b,
                          Array<double> &D, MPI::Comm &comm) {
        MPQC_PROFILE_LINE;
        // db = d*B
        double db = 0;
        //#pragma omp parallel reduction(+:db)
        foreach (auto j, beta.block(128)) {
            //MPQC_PROFILE_LINE;
            db += dot<double>(D(alpha,j), b(alpha,j));
        }
        comm.sum(db);

        // D = D - db*b;
        double dd = 0;
        // #pragma omp parallel reduction(+:dd)
        foreach (auto j, beta.block(128)) {
            //MPQC_PROFILE_LINE;
            Matrix bj = b(alpha,j);
            Matrix Dj = D(alpha,j);
            Dj -= db*bj;
            dd += dot(Dj, Dj);
            D(alpha,j) << Dj;
        }
        comm.sum(dd);

        // D = D/||D||
        dd = 1 / sqrt(dd);
        //#pragma omp parallel
        foreach (auto j, beta.block(128)) {
            //MPQC_PROFILE_LINE;
            Matrix Dj = D(alpha,j);
            Dj *= dd;
            //D(alpha,b) *= d;
            D(alpha,j) << Dj;
        }
        return db * dd;
    }

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

        mpqc::Array<double> C("ci.C", ci.dims, comm);
        mpqc::Array<double> D("ci.D", ci.dims, comm);

        // mpqc::Array<double> C("ci.C", ci.dims, ARRAY_FILE, comm);
        // mpqc::Array<double> D("ci.D", ci.dims, ARRAY_FILE, comm);

        comm.barrier();

        std::vector<double> E;

        // #ifdef HAVE_PYTHON 
        //         Python py;
        //         try {
        //             py.exec("import numpy as np\n"
        //                     "import h5py\n"
        //                     "import matplotlib.pyplot as plt\n");
        //             py.exec("mpqc.io = h5py.File(h5py.h5f.FileID(%i))\n"
        //                     "mpqc.ci = h5py.Group(h5py.h5g.GroupID(%i))\n"
        //                     "na,nb = (%i,%i)\n" %
        //                     Python::tuple(io.file(), io.id(),
        //                                   alpha.size(), beta.size()));
        //         }
        //         catch (...) {
        //             py.error();
        //         }
        // #endif // HAVE_PYTHON

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
                C(alpha, local).read(io.b[it]);
                C.sync();

                sigma(ci, h, V, C, D);
                D.sync();

                D(alpha, local).write(io.Hb[it]);
                comm.barrier();
            }

            MPQC_PROFILE_RESET;

            // augment G matrix
            {
                Vector g = Vector::Zero(M);
                foreach (auto rb, local.block(128)) {
                    MPQC_PROFILE_LINE;
                    Matrix c(alpha.size(), rb.size());
                    const Matrix &s = D(alpha, rb);
                    for (int j = 0; j < M; ++j) {
                        io.b(alpha,rb,j) >> c;
                        double q = 0; //dot(c,s);
#pragma omp parallel for schedule(dynamic,1) reduction(+:q)
                        for (int b = 0; b < rb.size(); ++b) {
                            q += dot<double>(c.col(b), s.col(b));
                        }
                        g(j) += q;
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
                foreach (auto rb, local.block(128)) {
                    MPQC_PROFILE_LINE;
                    Matrix v(alpha.size(), rb.size());
                    Matrix d(alpha.size(), rb.size());
                    d.fill(0);
                    for (int i = 0; i < M; ++i) {
                        MPQC_PROFILE_LINE;
                        io.b(alpha,rb,i) >> v;
                        double r = -a(i,k)*lambda(k);
#pragma omp parallel for schedule(dynamic,1)
                        for (int b = 0; b < rb.size(); ++b) {
                            d.col(b) += r*v.col(b);
                        }
                    }
                    for (int i = 0; i < M; ++i) {
                        MPQC_PROFILE_LINE;
                        io.Hb(alpha,rb,i) >> v;
#pragma omp parallel for schedule(dynamic,1)
                        for (int b = 0; b < rb.size(); ++b) {
                            d.col(b) += a(i,k)*v.col(b);
                        }
                    }
                    D(alpha,rb) << d;
                }
                D.sync();

                iters[it].E = lambda[0];
                iters[it].D = norm(D, comm);

                if (comm.rank() == 0) {
                    double dc = fabs(iters[it - 1].D - iters[it].D);
                    double de = fabs(iters[it - 1].E - iters[it].E);
                    sc::ExEnv::out0()
                        << sc::indent
                        << sc::scprintf("CI iter. %i, E=%12.8f, del.E=%e, del.C=%e\n",
                                        (int) it, lambda[0] + ci.e_ref, de, dc);
                }

                // preconditioner
                {
                    MPQC_PROFILE_LINE;
                    double dd = 0;
                    foreach (auto rb, local.block(128)) {
                        Matrix d = D(alpha,rb);

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

                        D(alpha,rb) << d;
                        dd += dot(d, d);
                    }
                    comm.sum(dd);
                    D.sync();
                    if (comm.rank() == 0)
                        symmetrize(D, 1, 1 / sqrt(dd));
                    D.sync();
                }

                // orthonormalize
                for (int i = 0; i < M; ++i) {
                    MPQC_PROFILE_LINE;
                    Array<double> &b = C;
                    b(alpha, local).read(io.b[i]);
                    orthonormalize(alpha, local, b, D, ci.comm);
                }
                D.sync();

                D(alpha, local).write(io.b[it + 1]);

            }

            // #ifdef HAVE_PYTHON 
            //             if (it > 0 && it%1 == 0) {
            //                 try {
            //                     py.exec("it = %i" % Python::tuple(M));
            //                     py.interactive();
            // 		    //py.exec("ci.items()");
            // 		}
            // 		catch (...) {
            // 		    py.error();
            // 		}
            // 	    }
            // #endif // HAVE_PYTHON

            //std::cout << "time: " << t << std::endl;

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
