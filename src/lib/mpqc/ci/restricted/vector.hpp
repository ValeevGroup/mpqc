#ifndef MPQC_CI_BLOCK_ARRAY_HPP
#define MPQC_CI_BLOCK_ARRAY_HPP

#include "mpqc/array.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/mpi.hpp"

namespace mpqc {
namespace ci {

    struct BlockArray  {
        BlockArray(const std::string &name, mpqc::matrix<size_t> blocks)
        {
            //std::cout << "m = " << m << ", n = " << n << std::endl;
        }

        struct Vector {
            size_t size() const {
                return data_.size();
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
                //data_(r2()).put(v.data());
            }
            mpqc::Vector get() const {
                return this->data_;
            }
        private:
            friend class BlockArray;
            Vector(Eigen::Map<Eigen::VectorXd> data)
                : data_(data)
            {
                //std::cout << "range r = " << r << std::endl;
            }
        private:
            Eigen::Map<Eigen::VectorXd> data_;
        };

        void sync() {
        }

    private:
        mpqc::Vector data_;
    };


}
}

#endif /* MPQC_CI_BLOCK_ARRAY_HPP */
