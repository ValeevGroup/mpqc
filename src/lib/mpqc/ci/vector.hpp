#ifndef MPQC_CI_VECTOR_HPP
#define MPQC_CI_VECTOR_HPP

#include "mpqc/array.hpp"
#include "mpqc/ci/full/ci.hpp"
#include "mpqc/mpi.hpp"
#include "mpqc/mpi/task.hpp"

namespace mpqc {
namespace ci {
namespace detail {

    template<typename T, typename U>
    void symmetrize(Array<T> &A, U &counter, size_t B) {
        assert(A.rank() == 2);
        assert(A.dims()[0] == A.dims()[1]);

        size_t N = A.dims()[0];
        matrix<T> Aij(B,B), Aji(B,B);

        for (size_t j = 0, ij = 0, next = counter++; j < N; j += B) {
            for (size_t i = 0; i <= j; i += B, ++ij) {

                if (ij != next) continue;
                next = counter++;

                range ri(i, std::min(i+B,N));
                range rj(j, std::min(j+B,N));

                Aij.resize(ri.size(), rj.size());
                Aji.resize(rj.size(), ri.size());
                A(ri,rj) >> Aij;
                A(rj,ri) >> Aji;
                Aij = Aij + Aji.transpose();
                Aji = Aij.transpose();
                A(ri,rj) << Aij;
                A(rj,ri) << Aji;                
            }
        }
    }

}
}
}


namespace mpqc {
namespace ci {

    template<class CI>
    struct Vector;

    template<class CI>
    inline void operator<<(Vector<CI> a, const mpqc::Vector &v) {
        a.put(v);
    }

    template<typename T>
    void symmetrize(mpqc::Array<T> &A, MPI::Comm comm, size_t B = 1024) {
        MPI::Task task(comm);
        detail::symmetrize(A, task, B);
    }

    template<typename T>
    void symmetrize(mpqc::Array<T> &A, size_t B = 1024) {
        size_t counter = 0;
        detail::symmetrize(A, counter, B);
    }

    template<class A>
    void symmetrize(Eigen::MatrixBase<A> &a) {
        MPQC_CHECK(a.rows() == a.cols());
        for (size_t j = 0; j < a.cols(); ++j) {
            for (size_t i = 0; i <= j; ++i) {
                a(i,j) += a(j,i);
                a(j,i) = a(i,j);
            }
        }
    }

    template<class A>
    void symmetrize(Eigen::MatrixBase<A> &a, double phase, double scale) {
        MPQC_CHECK(a.rows() == a.cols());
        for (size_t j = 0; j < a.cols(); ++j) {
            a(j,j) *= scale;
            for (size_t i = 0; i < j; ++i) {
                a(i,j) = scale*a(i,j); // + a(j,i));
                a(j,i) = phase*a(i,j);
            }
        }
    }

    void symmetrize(mpqc::Array<double> &A,
                    double phase, double scale,
                    MPI::Comm comm,
                    size_t block = 512) {
        MPQC_CHECK(A.dims()[0] == A.dims()[1]);
        Matrix a;
        std::vector<range> r = range::split(range(0, A.dims()[1]), block);
        MPI::Task task(comm);
        int next = task++;
        int ij = 0;
        for (auto rj = r.begin(); rj < r.end(); ++rj) {
            //std::cout << *rj << std::endl;
            size_t nj = rj->size();
            for (auto ri = r.begin(); ri <= rj; ++ri) {
                if (ij++ != next) continue;
                next = task++;
                size_t ni = ri->size();
                a.resize(ni, nj);
                if (ri == rj) {
                    A(*ri,*ri) >> a;
                    symmetrize(a, phase, scale);
                    A(*ri,*ri) << a;
                }
                else {
                    A(*ri,*rj) >> a;
                    a *= scale;
                    A(*ri,*rj) << a;
                    a *= phase;
                    A(*rj,*ri) << Matrix(a.transpose());
                }
            }
        }
    }

    template<class CI>
    double norm(ci::Vector<CI> &v, const MPI::Comm &comm, range local, size_t block) {
        double n = 0;
        foreach (auto r, local.block(block)) {
            n += mpqc::Vector(v(r)).norm();
        }
        comm.sum(n);
        return n;
    }

    /// Schmidt orthogonalization
    /// d' = normalized(d - <d,b>*b)
    /// @return <d,b>*<d',d'>
    template<class CI>
    double orthonormalize(ci::Vector<CI> &b, ci::Vector<CI> &D,
                          MPI::Comm &comm,
                          range local, size_t block) {
        // N.B. foreach doesn't work with openmp
        // most likely not needed - I/O bound
        // db = d*B
        double db = 0;
        ////#pragma omp parallel reduction(+:db)
        foreach (auto rj, local.block(block)) {
            mpqc::Vector Dj = D(rj);
            mpqc::Vector bj = b(rj);
            db += Dj.dot(bj);
        }
        comm.sum(db);

        // D = D - db*b;
        double dd = 0;
        ////#pragma omp parallel reduction(+:dd)
        foreach (auto rj, local.block(block)) {
            mpqc::Vector bj = b(rj);
            mpqc::Vector Dj = D(rj);
            Dj -= db*bj;
            dd += Dj.dot(Dj);
            D(rj).put(Dj);
        }
        comm.sum(dd);

        // D = D/||D||
        dd = 1/sqrt(dd);
        ////#pragma omp parallel
        foreach (auto rj, local.block(block)) {
            mpqc::Vector Dj = D(rj);
            Dj *= dd;
            D(rj).put(Dj);
        }
        return db*dd;
    }


}
}

#endif /* MPQC_CI_VECTOR_HPP */
