#ifndef MPQC_CI_ARRAY_HPP
#define MPQC_CI_ARRAY_HPP

#include "mpqc/array.hpp"
#include "mpqc/ci/full.hpp"
#include "mpqc/mpi.hpp"

namespace mpqc {
namespace ci {

    template<class CI>
    struct Array;

    template<class CI>
    inline void operator<<(Array<CI>::Vector a, const mpqc::Vector &v) {
        a.put(v);
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

    void symmetrize(mpqc::Array<double> &A, double phase, double scale) {
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

    template<class CI>
    double norm(ci::Array<CI> v, const MPI::Comm &comm, range local, size_t block) {
        double n = 0;
        foreach (auto r, local.block(block)) {
            n += Vector(v.vector(r)).norm();
        }
        comm.sum(n);
        return n;
    }

    /**
     //     double db = dot(D, b);
     //     D = D - db*b;
     //     D *= 1/D.norm();
     //     @return <d,b>/||d||
     */

    template<class CI>
    double orthonormalize(ci::Array<CI> b, ci::Array<CI> D,
                          MPI::Comm &comm,
                          range local, size_t block) {
        // db = d*B
        double db = 0;
        //#pragma omp parallel reduction(+:db)
        foreach (auto rj, local.block(block)) {
            Vector Dj = D.vector(rj);
            Vector bj = b.vector(rj);
            db += Dj.dot(bj);
        }
        comm.sum(db);

        // D = D - db*b;
        double dd = 0;
        // #pragma omp parallel reduction(+:dd)
        foreach (auto rj, local.block(block)) {
            Vector bj = b.vector(rj);
            Vector Dj = D.vector(rj);
            Dj -= db*bj;
            dd += Dj.dot(Dj);
            D.vector(rj) << Dj;
        }
        comm.sum(dd);

        // D = D/||D||
        dd = 1/sqrt(dd);
        //#pragma omp parallel
        foreach (auto rj, local.block(block)) {
            Vector Dj = D.vector(rj);
            Dj *= dd;
            D.vector(rj) << Dj;
        }
        return db*dd;
    }


}
}

#endif /* MPQC_CI_ARRAY_HPP */
