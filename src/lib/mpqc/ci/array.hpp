#ifndef MPQC_CI_ARRAY_HPP
#define MPQC_CI_ARRAY_HPP

#include "mpqc/array.hpp"
#include "mpqc/mpi.hpp"

namespace mpqc {
namespace ci {

    struct Array  {
        Array(const std::string &name, size_t m, size_t n, MPI::Comm comm)
            : data_(name, extents(m,n), comm),
              alpha_(m), beta_(n)
              
        {
            //std::cout << "m = " << m << ", n = " << n << std::endl;
        }

        struct Vector {
            // not correct in parallel runs, fix it
            Vector operator()(range r) {
                MPQC_CHECK(r.size() < r_.size());
                auto begin = r_.front();
                return Vector(data_, alpha_, beta_,
                              range(r.front() + begin, r.back()+1+begin));
            }
            size_t size() const {
                return r_.size();
            }
            void read(mpqc::File::Dataspace<double> ds) {
                mpqc::Vector v(this->size());
                ds.read(v.data());
                this->put(v);
            }
            void write(mpqc::File::Dataspace<double> ds) {
                const mpqc::Vector &v = this->get();
                ds.write(v.data());
            }
            operator mpqc::Vector() const {
                return this->get();
            }
            void put(const mpqc::Vector &v) {
                MPQC_CHECK(v.size() == this->size());
                data_(r2()).put(v.data());
            }
            mpqc::Vector get() const {
                mpqc::Vector v(this->size());
                data_(r2()).get(v.data());
                return v;
            }
        private:
            friend class Array;
            Vector(mpqc::Array<double> data, size_t alpha, size_t beta, range r)
                : data_(data), alpha_(alpha), beta_(beta), r_(r)
            {
                //std::cout << "range r = " << r << std::endl;
            }
            std::vector<range> r2() const {
                MPQC_CHECK(r_.front()%alpha_ == 0);
                MPQC_CHECK(r_.size()%alpha_ == 0);
                std::vector<range> r2;
                r2.push_back(range(r_.front()%beta_, 1+r_.back()%beta_));
                r2.push_back(range(r_.front()/beta_, 1+r_.back()/beta_));
                return r2;
                //std::cout << "r_[] = " << r_[0] << "," << r_[1] << std::endl;
            }
        private:
            mpqc::Array<double> data_;
            size_t alpha_, beta_;
            range r_;
        };

        Vector vector(range r) {
            return Vector(this->data_, alpha_, beta_, r);
        }

        mpqc::Array<double> array(range ri, range rj) {
            //std::cout << "array " << ri << "," << rj << std::endl;
            return data_(ri,rj);
        }

        mpqc::Array<double>& array() {
            return data_;
        }

        void sync() {
            data_.sync();
        }

    private:
        size_t alpha_, beta_;
        mpqc::Array<double> data_;
        static std::vector<size_t> extents(size_t m, size_t n) {
            std::vector<size_t> v{m,n};
            //std::cout << "v=" << v.size() << std::endl;
            return v;
        }
    };

    inline void operator<<(Array::Vector a, const mpqc::Vector &v) {
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

    double norm(ci::Array v, const MPI::Comm &comm, range local, size_t block) {
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
    double orthonormalize(ci::Array b, ci::Array D,
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
