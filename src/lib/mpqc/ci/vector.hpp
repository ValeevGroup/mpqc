#ifndef MPQC_CI_VECTOR_HPP
#define MPQC_CI_VECTOR_HPP

#include "mpqc/ci/subspace.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/array.hpp"
#include "mpqc/mpi.hpp"

#include "mpqc/utility/profile.hpp"

namespace mpqc {
namespace ci {

    /// @addtogroup CI
    /// @{

    /// Block CI Vector, with 1-d (vector) and 2-d (matrix) access.
    struct Vector : boost::noncopyable {

        /// Construct vector
        /// @param name vector name
        /// @param G vector grid (blocking and sparsity)
        /// @param comm MPI communicator
        /// @param incore store vector in core (incore=true) or in file
        Vector(std::string name, const SubspaceGrid &G, MPI::Comm comm, bool incore) {
            blocks_.resize(G.alpha().size(), G.beta().size());
            size_t dets = 0;
            for (size_t j = 0; j < G.beta().size(); ++j) {
                for (size_t i = 0; i < G.alpha().size(); ++i) {
                    Block b;
                    b.rows = G.alpha().at(i).size();
                    b.cols = G.beta().at(j).size();
                    b.begin = dets;
                    b.allowed = G.allowed(i,j);
                    if (b.allowed) {
                        dets += b.rows*b.cols;
                    }
                    blocks_(i,j) = b;
                }
            }
            MPQC_CHECK(dets == G.dets());
            BOOST_FOREACH (const auto &s, G.alpha())
                this->alpha_.push_back(s);
            BOOST_FOREACH (const auto &s, G.beta())
                this->beta_.push_back(s);
            std::vector<size_t> extents;
            extents.push_back(dets);
            this->array_ = (incore)
                ? Array<double>(name, extents, ARRAY_CORE, comm)
                : Array<double>(name, extents, ARRAY_FILE, comm);
        }

        
        /// 1-d vector sub-block
        struct Block1d {
        private:
            Array<double> array_;
            friend void operator<<(Vector::Block1d block, const mpqc::Vector &v);
            friend void operator>>(Vector::Block1d block, mpqc::Vector &v);
        private:
            friend class Vector;
            Block1d(Array<double> array) : array_(array) {}
        };

        /// Returns 1-d sub-block of vector
        /// @param r sub-block range
        Block1d operator()(range r) {
            return Block1d(this->array_(r));
        }


        /// 2-d vector sub-block
        struct Block2d : mpqc::Matrix::Assignable {
            /// Assign matrix to this sub-block
            void operator=(const mpqc::Matrix &m) const {
                MPQC_CHECK(m.rows() == this->rows_);
                MPQC_CHECK(m.cols() == this->cols_);
                this->array_ << m;
            }
            /// mpqc::Matrix::Assignable::assign_to implementation
            void assign_to(mpqc::Matrix &m) const {
                m.resize(this->rows_, this->cols_);
                this->array_ >> m;
            }
        private:
            friend class Vector;
            Block2d(Array<double> array, size_t rows, size_t cols)
                : array_(array), rows_(rows), cols_(cols) {}
        private:
            Array<double> array_;
            size_t rows_, cols_;
        };

        /// Returns 2-d sub-block of vector
        Block2d operator()(mpqc::range A, mpqc::range B) {
            Block b = this->block(A,B);
            return Block2d(this->array_(b.range()), b.rows, b.cols);
        }

        void sync() {
            this->array_.sync();
        }

    private:

        struct Block {
            size_t rows, cols;
            size_t begin;
            bool allowed;
            mpqc::range range() const {
                return mpqc::range(begin, begin+rows*cols);
            }
        };

        std::vector<mpqc::range> alpha_;
        std::vector<mpqc::range> beta_;
        mpqc::Array<double> array_;
        mpqc::matrix<Block> blocks_;

    private:

        /// resolve (A,B) into block
        Block block(mpqc::range A, mpqc::range B) const {
            size_t i = range_to_block(A, this->alpha_);
            size_t j = range_to_block(B, this->beta_);
            Block b = this->blocks_(i,j);
            if (!b.allowed) {
                throw MPQC_EXCEPTION("ci::Vector: block (%s,%s) is not allowed\n",
                                     string_cast(A).c_str(), string_cast(B).c_str());
            }
            return b;
        }

        /// resolve subspace to position in subspace vector S.
        static size_t range_to_block(mpqc::range r, const std::vector<mpqc::range> &R) {
            auto it = std::find(R.begin(), R.end(), r);
            if (it == R.end()) {
                throw MPQC_EXCEPTION("ci::Vector: range r=(%s) does not map an exact block",
                                     string_cast(mpqc::range(r)).c_str());
            }
            return it-R.begin();
        }

    };

    /// Set vector block to v
    void operator<<(Vector::Block1d block, const mpqc::Vector &v) {
        block.array_ << v;
    }

    /// Get vector block into v
    void operator>>(Vector::Block1d block, mpqc::Vector &v) {
        block.array_ >> v;
    }

    /// Compute CI vector norm
    /// @param[in] V CI Vector
    /// @param local vector of sub-blocks the node will process (implies parallelization)
    /// @param comm MPI comm
    double norm(ci::Vector &V,
                const std::vector<mpqc::range> &local,
                const MPI::Comm &comm) {
        double norm = 0;
        BOOST_FOREACH (auto r, local) {
            mpqc::Vector v(r.size());
            V(r) >> v;
            norm += v.norm();
        }
        comm.sum(norm);
        return norm;
    }

    /// Schmidt orthogonalization
    /// d' = normalized(d - <d,b>*b)
    /// @param[in] b orthonormal CI Vector (notice ortho*normal* - MUST be norm-1)
    /// @param[in,out] d CI Vector to orthonormalize
    /// @param local vector of sub-blocks the node will process (implies parallelization)
    /// @param comm MPI comm
    /// @return <d,b>*<d',d'>
    double orthonormalize(ci::Vector &b, ci::Vector &D,
                          const std::vector<mpqc::range> &local,
                          const MPI::Comm &comm)
    {
        MPQC_PROFILE_LINE;
        double db = 0;
#pragma omp parallel for reduction(+:db)
        for (int j = 0; j < local.size(); ++j) {
            mpqc::range rj = local.at(j);
            mpqc::Vector Dj(rj.size()), bj(rj.size());
            D(rj) >> Dj;
            b(rj) >> bj;
            db += Dj.dot(bj);
        }
        comm.sum(db);

        // D = D - db*b;
        double dd = 0;
#pragma omp parallel for reduction(+:dd)
        for (int j = 0; j < local.size(); ++j) {
            mpqc::range rj = local.at(j);
            mpqc::Vector Dj(rj.size()), bj(rj.size());
            D(rj) >> Dj;
            b(rj) >> bj;
            Dj -= db*bj;
            dd += Dj.dot(Dj);
            D(rj) << Dj;
        }
        comm.sum(dd);

        // D = D/||D||
        dd = 1/sqrt(dd);
#pragma omp parallel for
        for (int j = 0; j < local.size(); ++j) {
            mpqc::range rj = local.at(j);
            mpqc::Vector Dj(rj.size());
            D(rj) >> Dj;
            Dj *= dd;
            D(rj) << Dj;
        }
        return db*dd;
    }

    /// @}

}
}

#endif /* MPQC_CI_VECTOR_HPP */
