#ifndef MPQC_CI_ARRAY_HPP
#define MPQC_CI_ARRAY_HPP

#include "mpqc/array.hpp"
#include "mpqc/mpi.hpp"

namespace mpqc {
namespace ci {

    struct Array  {
        Array(const std::string &name, size_t m, size_t n, MPI::Comm comm)
            : cols_(m), rows_(n),
              data_(name, extents(m,n))
        {
            std::cout << "m = " << m << ", n = " << n << std::endl;
        }

        struct Vector {
            void read(mpqc::File::Dataspace<double> ds) {
                mpqc::Vector v(this->size_);
                ds.read(v.data());
                data_.put(v.data());
            }
            void write(mpqc::File::Dataspace<double> ds) {
                const mpqc::Vector &v = *this;
                ds.write(v.data());
            }
            operator mpqc::Vector() const {
                mpqc::Vector v(this->size_);
                data_.get(v.data());
                return v;
            }
            void put(const mpqc::Vector &v) {
                data_.put(v.data());
            }
        private:
            friend class Array;
            Vector(Array *data, range r)
            {
                //std::cout << "range r = " << r << std::endl;
                MPQC_CHECK(r.front()%data->rows_ == 0);
                MPQC_CHECK(r.size()%data->rows_ == 0);
                range ri = range(r.front()%data->cols_, 1+r.back()%data->cols_);
                range rj = range(r.front()/data->cols_, 1+r.back()/data->cols_);
                this->data_ = data->array(ri, rj);
                this->size_ = r.size();
            }
            mpqc::Array<double> data_;
            size_t size_;
        };

        Vector vector(range r) {
            return Vector(this, r);
        }

        mpqc::Array<double> array(range ri, range rj) {
            return data_(ri,rj);
        }

        mpqc::Array<double>& array() {
            return data_;
        }

        void sync() {
            data_.sync();
        }

    private:
        size_t cols_, rows_;
        mpqc::Array<double> data_;
        static std::vector<size_t> extents(size_t m, size_t n) {
            std::vector<size_t> v{m,n};
            std::cout << "v=" << v.size() << std::endl;
            return v;
        }
    };

    inline void operator<<(Array::Vector a, const mpqc::Vector &v) {
        a.put(v);
    }

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

    /**
     //     double db = dot(D, b);
     //     D = D - db*b;
     //     D *= 1/D.norm();
     //     @return <d,b>/||d||
     */
    double orthonormalize(range alpha, range beta,
                          const mpqc::Array<double> &b,
                          mpqc::Array<double> &D,
                          MPI::Comm &comm) {
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


}
}

#endif /* MPQC_CI_ARRAY_HPP */
