#ifndef MPQC_CI_VECTOR_HPP
#define MPQC_CI_VECTOR_HPP
#include "mpqc/ci/subspace.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/array.hpp"
#include "mpqc/mpi.hpp"

namespace mpqc {
namespace ci {

    /// Block CI Vector.
    /// This Vector should be renamed Block Sparse Vector as it is applicable beyond R-CI.
    /// Matrix access spanning more than one block is not allowed.
    struct BlockVector : boost::noncopyable {

        BlockVector(std::string name, const SubspaceGrid &G, MPI::Comm comm) {
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
            foreach (const auto &s, G.alpha())
                this->alpha_.push_back(s);
            foreach (const auto &s, G.beta())
                this->beta_.push_back(s);
            std::vector<size_t> extents;
            extents.push_back(dets);
            this->array_ = Array<double>(name, extents, comm);
        }


        struct Block1d {
            Block1d(Array<double> array)
                : array_(array) {}
            friend void operator<<(BlockVector::Block1d block, const mpqc::Vector &v);
            friend void operator>>(BlockVector::Block1d block, mpqc::Vector &v);
        private:
            Array<double> array_;
        };

        Block1d operator()(range r) {
            return Block1d(this->array_(r));
        }


        struct Block2d : mpqc::Matrix::Assignable {
            Block2d(Array<double> array, size_t rows, size_t cols)
                : array_(array), rows_(rows), cols_(cols) {}
            void operator=(const mpqc::Matrix &m) const {
                MPQC_CHECK(m.rows() == this->rows_);
                MPQC_CHECK(m.cols() == this->cols_);
                this->array_ << m;
            }
            void assign_to(mpqc::Matrix &m) const {
                m.resize(this->rows_, this->cols_);
                this->array_ >> m;
            }
        private:
            Array<double> array_;
            size_t rows_, cols_;
        };

        Block2d operator()(mpqc::range A, mpqc::range B) {
            Block b = this->block(A,B);
            return Block2d(this->array_(b.range()), b.rows, b.cols);
        }


        void sync() {}

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


    void operator<<(BlockVector::Block1d block, const mpqc::Vector &v) {
        block.array_ << v;
    }

    void operator>>(BlockVector::Block1d block, mpqc::Vector &v) {
        block.array_ >> v;
    }


    double norm(ci::BlockVector &V,
                const std::vector<mpqc::range> &local,
                const MPI::Comm &comm) {
        double norm = 0;
        foreach (auto r, local) {
            mpqc::Vector v(r.size());
            V(r) >> v;
            norm += v.norm();
        }
        comm.sum(norm);
        return norm;
    }

    /// Schmidt orthogonalization
    /// d' = normalized(d - <d,b>*b)
    /// @return <d,b>*<d',d'>
    double orthonormalize(ci::BlockVector &b, ci::BlockVector &D,
                          const std::vector<mpqc::range> &local,
                          const MPI::Comm &comm)
    {
        // N.B. foreach doesn't work with openmp
        // most likely not needed - I/O bound
        // db = d*B
        double db = 0;
        ////#pragma omp parallel reduction(+:db)
        foreach (auto rj, local) {
            mpqc::Vector Dj(rj.size()), bj(rj.size());
            D(rj) >> Dj;
            b(rj) >> bj;
            db += Dj.dot(bj);
        }
        comm.sum(db);

        // D = D - db*b;
        double dd = 0;
        ////#pragma omp parallel reduction(+:dd)
        foreach (auto rj, local) {
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
        ////#pragma omp parallel
        foreach (auto rj, local) {
            mpqc::Vector Dj(rj.size());
            D(rj) >> Dj;
            Dj *= dd;
            D(rj) << Dj;
        }
        return db*dd;
    }


}
}

#endif /* MPQC_CI_VECTOR_HPP */
