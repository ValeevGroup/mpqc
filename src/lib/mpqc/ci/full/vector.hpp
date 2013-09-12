#ifndef MPQC_CI_FULL_VECTOR_HPP
#define MPQC_CI_FULL_VECTOR_HPP

#include "mpqc/ci/full/ci.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/array.hpp"

namespace mpqc {
namespace ci {

    template<class Index>
    struct Vector< Full<Index> > {

        typedef mpqc::Array<double> Array;

        Vector(std::string name, const CI< Full<Index> > &ci)
            : array_(name, dims(ci), ci.comm),
              alpha_(ci.alpha.size()),
              beta_(ci.beta.size())
        {
        }

        struct Slice {
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
            Slice(Vector &v, range r)
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

        Slice operator()(range r) {
            return Slice(*this, r);
        }

        operator Array&() {
            return this->array_;
        }

        void sync() {
            array_.sync();
        }

    private:
        Array array_;
        size_t alpha_, beta_;
        static std::vector<size_t> dims(const CI< Full<Index> > &ci) {
            std::vector<size_t> v{ ci.alpha.size(), ci.beta.size() };
            //std::cout << "v=" << v.size() << std::endl;
            return v;
        }
    };

}
}

#endif /* MPQC_CI_FULL_VECTOR_HPP */
