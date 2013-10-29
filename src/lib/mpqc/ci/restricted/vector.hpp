#ifndef MPQC_CI_RESTRICTED_VECTOR_HPP
#define MPQC_CI_RESTRICTED_VECTOR_HPP

#include "mpqc/ci/subspace.hpp"
#include "mpqc/ci/restricted/ci.hpp"
#include "mpqc/ci/vector.hpp"

#include <boost/lambda/lambda.hpp>

namespace mpqc {
namespace ci {

    /// Restricted CI Vector.
    /// This Vector should be renamed Block Sparse Vector as it is applicable beyond R-CI.
    /// Matrix access spanning more than one block is not allowed.
    template<>
    struct Vector<Restricted> : boost::noncopyable {

        Vector(std::string name, const SubspaceGrid &grid) {
            initialize(name, grid);
        }

        template<class Index>
        Vector(std::string name, const CI<Restricted, Index> &ci)
            : grid_(ci.subspace)
        {
            initialize(name, ci.subspace);
            MPQC_CHECK(vector_.size() == ci.dets());
        }

        struct Block1d {
        private:
            friend class Vector;
            Block1d(mpqc::Vector &v, range r)
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
                assert(data_.size() == v.size());
                for (int i = 0; i < v.size(); ++i) {
                    data_[i] = v[i];
                }
            }
        private:
            range r_;
            Eigen::Map<Eigen::VectorXd> data_;
        };

        Block1d operator()(range r) {
            return Block1d(this->vector_, r);
        }

        template<class Matrix>
        struct Block2d : mpqc::Matrix::Assignable {
        public:
            typedef Eigen::Map<Matrix> Map;
            typedef typename Map::PointerType PointerType;

            template<class B>
            void operator=(const Eigen::DenseBase<B> &b) {
                namespace lambda = boost::lambda;
                apply<Matrix>(lambda::_1 = lambda::_2, *this, b);
            }
            
            void assign_to(mpqc::Matrix &m) const {
                namespace lambda = boost::lambda;
                m.resize(this->ri_.size(), this->rj_.size());
                apply<const Matrix>(lambda::_2 = lambda::_1, *this, m);
            }

        private:

            friend class Vector;
            mpqc::range ri_, rj_;
            mpqc::matrix<PointerType> tiles_;
            std::vector< Subspace<Alpha> > alpha_;
            std::vector< Subspace<Beta> > beta_;

            explicit Block2d(mpqc::range ri, mpqc::range rj,
                             const mpqc::matrix<PointerType> &tiles,
                             const std::vector< Subspace<Alpha> > &alpha,
                             const std::vector< Subspace<Beta> > &beta)
                : ri_(ri), rj_(rj), alpha_(alpha), beta_(beta), tiles_(tiles) {}

            template<class M, class F, class A_, class B_>
            static void apply(F f, A_ &A, B_ &B) { 
                assert(B.size());
                range ri = A.ri_, rj = A.rj_;
                for (int j = 0, Y = 0; j < A.beta_.size(); ++j) {
                    auto beta = A.beta_.at(j);
                    auto y = (rj & beta);
                    for (int i = 0, X = 0; i < A.alpha_.size(); ++i) {
                        auto alpha = A.alpha_.at(i);
                        auto x = (ri & alpha);

                        struct { size_t X, Y; } offset;
                        offset.X = x.begin() - alpha.begin();
                        offset.Y = y.begin() - beta.begin();

                        Eigen::Map<M> Aij(A.tiles_(i,j), alpha.size(), beta.size());

                        // printf("Aij(%lu,%lu)\n", Aij.rows(), Aij.cols());
                        // printf("offset.X=%i, offset.Y=%i x.size()=%i, y.size()=%i\n",
                        //        offset.X, offset.Y, x.size(), y.size());
                        // printf("B(%lu,%lu)\n", B.rows(), B.cols());
                        // printf("X=%i, Y=%i, x.size()=%i, y.size()=%i\n",
                        //        X, Y, x.size(), y.size());

                        auto aij = Aij.block(offset.X, offset.Y, x.size(), y.size());
                        auto bij = B.block(X, Y, x.size(), y.size());
                        f(aij, bij);

                        X += x.size();
                    }
                    Y += y.size();
                }
            }
            
        };

        typedef Eigen::Block< Eigen::Map<Eigen::MatrixXd> > EigenBlock;
        typedef Eigen::Block< Eigen::Map<const Eigen::MatrixXd> > ConstEigenBlock;

        Block2d<const Eigen::MatrixXd> operator()(range ri, range rj) const {
            return block<const Eigen::MatrixXd>(ri, rj, *this);
        }

        Block2d<Eigen::MatrixXd> operator()(range ri, range rj) {
            return block<Eigen::MatrixXd>(ri, rj, *this);
        }

        void sync() {}

        void symmetrize(size_t block = 512);

        void symmetrize(double phase, double scale, size_t block = 512);

    private:

        mpqc::matrix<size_t> blocks_;
        SubspaceGrid grid_;
        mpqc::Vector vector_;

    private:

        struct Block {
            size_t offset;
            mpqc::range alpha, cols;
        };

        /// Initialize Vector according to grid
        void initialize(std::string name, const SubspaceGrid &G) {
            this->grid_ = G;
            blocks_.resize(grid_.alpha().size(), grid_.beta().size());
            for (size_t j = 0, offset = 0; j < blocks_.cols(); ++j) {
                for (size_t i = 0; i < blocks_.rows(); ++i) {
                    blocks_(i,j) = size_t(-1);
                    if (!grid_.allowed(i,j)) continue;
                    blocks_(i,j) = offset;
                    offset += grid_.alpha(i).size()*grid_.beta(j).size();
                }
            }
            //printf("dets = %lu, ci.dets = %lu\n", dets, ci.dets);
            vector_.resize(grid_.dets());
            vector_.fill(0);
        }

        Eigen::Map<Eigen::MatrixXd> block(int i, int j) {
            typedef Eigen::MatrixXd Matrix;
            const SubspaceGrid &grid = this->grid_;
            if (!grid.allowed(i,j))
                throw MPQC_EXCEPTION("ci::Vector<Restricted> forbidden block (%i,%i)", i, j);
            return Eigen::Map<Matrix>(this->vector_.data() + this->blocks_(i,j),
                                      grid.alpha(i).size(), grid.beta(j).size());
        }        

        /// Resolve row/col range to containing 2-d subspace block
        template<class Matrix, class Vector>
        static Block2d<Matrix> block(range ri, range rj, Vector &V) {
            auto x = intersections(ri, V.grid_.alpha());
            auto y = intersections(rj, V.grid_.beta());
            mpqc::matrix<typename Block2d<Matrix>::PointerType> tiles(x.size(), y.size());
            for (int j = 0; j < y.size(); ++j) {
                for (int i = 0; i < x.size(); ++i) {
                    auto xi = x.at(i);
                    auto xj = y.at(j);
                    if (!V.grid_.allowed(xi,xj)) {
                        throw MPQC_EXCEPTION("ci::Vector<Restricted> range (%s,%s) "
                                             "spans forbidden subspace (%i,%i)",
                                             string_cast(ri).c_str(),
                                             string_cast(rj).c_str(),
                                             xi, xj);
                    }
                    tiles(i,j) = V.vector_.data() + V.blocks_(xi,xj);
                }
            }
            return Block2d<Matrix>(ri, rj, tiles,
                                   slice(V.grid_.alpha(), x),
                                   slice(V.grid_.beta(), y));
        }

        /// resolve range to indices in subspace vector R.
        template <class Spin>
        static std::vector<int> intersections(range r, const std::vector< Subspace<Spin> > &R) {
            std::vector<int> intersections;
            for (int i = 0; i < R.size(); ++i) {
                range x = (r & R.at(i));
                // no intersection
                if (!x.empty()) intersections.push_back(i);
            }
            return intersections;
        }

        template<class Spin>
        static std::vector< Subspace<Spin> > slice(const std::vector< Subspace<Spin> > &V,
                                                const std::vector<int> &indices) {
            std::vector< Subspace<Spin> > slice;
            foreach (int i, indices) {
                slice.push_back(V.at(i));
            }
            return slice;
        }

    };

    void Vector<Restricted>::symmetrize(size_t block) {
        MPQC_CHECK(this->grid_.symmetric());
        for (int j = 0; j < this->grid_.beta().size(); ++j) {
            for (int i = 0; i < j; ++i) {
                if (!this->grid_.allowed(i,j)) continue;
                this->block(i,j) += this->block(j,i).transpose();
                this->block(j,i)  = this->block(i,j).transpose();
            }
            if (!this->grid_.allowed(j,j)) continue;
            auto bjj = this->block(j,j);
            ci::symmetrize(bjj);
        }
    }

    void Vector<Restricted>::symmetrize(double phase, double scale, size_t block) {
        MPQC_CHECK(this->grid_.symmetric());
        for (int j = 0; j < this->grid_.beta().size(); ++j) {
            for (int i = 0; i < j; ++i) {
                if (!this->grid_.allowed(i,j)) continue;
                auto bij = this->block(i,j);
                auto bji = this->block(j,i);
                bij *= scale;
                bji = phase*(bij.transpose());
            }
            if (!this->grid_.allowed(j,j)) continue;
            auto bjj = this->block(j,j);
            ci::symmetrize(bjj, phase, scale);
        }
    }


}
}

#endif /* MPQC_CI_RESTRICTED_VECTOR_HPP */
