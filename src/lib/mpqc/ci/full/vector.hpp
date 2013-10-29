#ifndef MPQC_CI_FULL_VECTOR_HPP
#define MPQC_CI_FULL_VECTOR_HPP

#include "mpqc/ci/full/ci.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/array.hpp"

namespace mpqc {
namespace ci {

    template<>
    struct Vector<Full> {

        typedef mpqc::Array<double> Array;

        template<class Index>
        Vector(std::string name, const CI<Full, Index> &ci)
            : array_(name, dims(ci), ci.comm),
              alpha_(ci.alpha.size()),
              beta_(ci.beta.size())
        {
        }

        struct Slice1d {
            size_t size() const {
                return range1d_.size();
            }
            operator mpqc::Vector() const {
                mpqc::Vector v(this->size());
                array_.get(v.data());
                return v;
            }
            void read(mpqc::File::Dataspace<double> ds) {
                io<READ>(ds);
            }

            void write(mpqc::File::Dataspace<double> ds) {
                io<WRITE>(ds);
            }

            void put(const mpqc::Vector &v) {
                MPQC_CHECK(v.size() == this->size());
                array_.put(v.data());
            }

        private:
            enum { READ = -1, WRITE = 1 };
            template<int OP>
            void io(mpqc::File::Dataspace<double> ds) {
                range alpha = this->range2d_[0];
                range beta = this->range2d_[1];
                // N.B. using buffer of 8MB, make adjustable
                size_t B = std::max<size_t>((1<<20)/alpha.size(), 1);
                double *buffer = new double[alpha.size()*B];
                auto begin = this->range1d_.front();
                foreach (auto b, beta.block(B)) {
                    size_t size = alpha.size()*b.size();
                    range r = range(begin, begin+size);
                    if (OP == READ) {
                        ds(r).read(buffer);
                        array_(alpha, b).put(buffer);
                    }
                    if (OP == WRITE) {
                        array_(alpha, b).get(buffer);
                        ds(r).write(buffer);
                    }
                    begin += size;
                }
                delete [] buffer;
            }

        private:
            friend class Vector;
            range range1d_;
            std::vector<range> range2d_;
            Array array_;
            Slice1d(Vector &v, range r)
                : range1d_(r), range2d_(range2d(v, r)),
                  array_(v.array_(range2d_)) {}
            static std::vector<range> range2d(const Vector &v, range r) {
                size_t alpha = v.alpha_;
                size_t beta = v.beta_;
                MPQC_CHECK(r.front()%alpha == 0);
                MPQC_CHECK(r.size()%alpha == 0);
                std::vector<range> r2;
                r2.push_back(range(r.front()%beta, 1+r.back()%beta));
                r2.push_back(range(r.front()/beta, 1+r.back()/beta));
                return r2;
            }
        };

        Slice1d operator()(range r) {
            return Slice1d(*this, r);
        }

        struct Slice2d : mpqc::Matrix::Assignable {
        public:
            void operator=(const mpqc::Matrix &m) {
                this->array_ << m;
            }
            void assign_to(mpqc::Matrix &m) const {
                m.resize(this->array_.dims()[0], this->array_.dims()[1]);
                this->array_ >> m;
            }
        private:
            friend class Vector;
            Array array_;
            Slice2d(const Array &a) : array_(a) {}
        };

        Slice2d operator()(range ri, range rj) {
            return this->array_(ri,rj);
        }

        Slice2d operator()(range ri, range rj) const {
            return this->array_(ri,rj);
        }

        // operator Array&() {
        //     return this->array_;
        // }

        void sync() {
            array_.sync();
        }

        void symmetrize(size_t block = 512) {
            ci::symmetrize(this->array_, this->array_.comm(), block);
        }

        void symmetrize(double phase, double scale, size_t block = 512) {
            const MPI::Comm &comm = this->array_.comm();
            ci::symmetrize(this->array_, phase, scale, comm, block);
        }

    private:
        Array array_;
        size_t alpha_, beta_;
        template<class Index>
        static std::vector<size_t> dims(const CI<Full, Index> &ci) {
            std::vector<size_t> v{ ci.alpha.size(), ci.beta.size() };
            //std::cout << "v=" << v.size() << std::endl;
            return v;
        }
    };

}
}

#endif /* MPQC_CI_FULL_VECTOR_HPP */
