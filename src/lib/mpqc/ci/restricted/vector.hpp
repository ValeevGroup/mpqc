#ifndef MPQC_CI_BLOCK_ARRAY_HPP
#define MPQC_CI_BLOCK_ARRAY_HPP

#include "mpqc/ci/restricted/ci.hpp"
#include "mpqc/ci/vector.hpp"
#include "mpqc/array.hpp"

namespace mpqc {
namespace ci {

    template<class Index>
    struct Vector< Restricted<Index> > {

        Vector(std::string name, const CI< Restricted<Index> > &ci)
            // : array_(name, dims(ci), ci.comm),
            //   alpha_(ci.alpha.size()),
            //   beta_(ci.beta.size())
        {
            data_.resize(ci.dets);
            data_.fill(0);
        }

        struct Slice {
        private:
            friend class Vector;
            Slice(mpqc::Vector &v, range r)
                : data_(v.data() + r.front(), r.size()), r_(r) {}
        public:
            size_t size() const {
                return r_.size();
            }
            operator mpqc::Vector() const {
                return this->data_;
            }
            void read(mpqc::File::Dataspace<double> ds) {
                ds.read(this->data_.data());
            }
            void write(mpqc::File::Dataspace<double> ds) {
                ds.write(this->data_.data());
            }
            void put(const mpqc::Vector &v) {
                data_ = v;
            }
        private:
            range r_;
            Eigen::Map<Eigen::VectorXd> data_;
        };

        Slice operator()(range r) {
            return Slice(this->data_, r);
        }

        void sync() {}

    private:
        mpqc::Vector data_;

    };


}
}

#endif /* MPQC_CI_BLOCK_ARRAY_HPP */
