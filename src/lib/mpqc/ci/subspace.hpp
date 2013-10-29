#ifndef MPQC_CI_SUBSPACE_HPP
#define MPQC_CI_SUBSPACE_HPP

#include <util/misc/formio.h>
#include "mpqc/ci/string.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/mpi.hpp"
#include "mpqc/file.hpp"

#include "mpqc/utility/foreach.hpp"
#include "mpqc/utility/check.hpp"
#include "mpqc/utility/exception.hpp"

#include <boost/mpl/int.hpp>

namespace mpqc {
namespace ci {

    /// @addtogroup CI
    /// @{

    /// Electron spin
    template<int S>
    struct Spin : boost::mpl::int_<S> {};

    typedef Spin<1> Alpha;
    typedef Spin<-1> Beta;
    typedef Spin<0> Any;

    /// A CI space, marked by Spin S and rank.
    template<class Spin>
    struct Space {
        explicit Space(int rank) : rank_(rank) {}
        int rank() const { return rank_; }
        bool operator==(const Space &s) const {
            return this->rank_ == s.rank_;
        }
    private:
        int rank_;
    };

    /// Compare two spaces by their rank
    template<class Spin>
    inline bool operator<(const Space<Spin> &a, const Space<Spin> &b) {
        return (a.rank() < b.rank());
    }

    /// A range of a space where all objects in the subspace
    /// range are assumed to have the same space rank.
    template<class Spin>
    struct Subspace : ci::Space<Spin>, mpqc::range {
        typedef ci::Space<Spin> Space;
        /// Constructs subspace
        /// @arg s Subspace space
        /// @arg Subspace range
        Subspace(Space s, mpqc::range r)
            : Space(s), mpqc::range(r) {}
        /// Compares two subspaces by their range and space
        bool operator==(const Subspace &a) const {
            return ((mpqc::range(a) == mpqc::range(*this)) && (Space(a) == Space(*this)));
        }
        /// test if index is in subspace range
        bool test(int idx) const {
            return ((*mpqc::range::begin() <= idx) && (idx < *mpqc::range::end()));
        }
        /// Split the subspace into a vector of subspaces
        std::vector<Subspace> split(size_t block) const {
            std::vector<Subspace> s;
            foreach (auto r, mpqc::range::split(*this, block)) {
                s.push_back(Subspace(*this, r));
            }
            return s;
        }
    };

    /// A subspace specialization that "discard" the spin parameter.
    /// The specialization intent is to represent subspace where spin is not applicable
    template<>
    struct Subspace<Any> : ci::Space<Any>, mpqc::range {
        /// Copies subspace while discarding the Spin parameter
        template<class Spin>
        Subspace(const Subspace<Spin> &s) 
            : ci::Space<Any>(s.rank()), mpqc::range(s) {}
        /// Compares two subspaces by their range and space
        bool operator==(const Subspace &b) const {
            return ((mpqc::range(*this) == mpqc::range(b)) && (Space<Any>(*this) == Space<Any>(b)));
        }
    };

    template<class Spin>
    std::vector< Subspace<Spin> > split(const std::vector< Subspace<Spin> > &V, size_t block) {
        std::vector< Subspace<Spin> > S;
        foreach (auto v, V) {
            foreach (auto s, v.split(block)) {
                S.push_back(s);
            }
        }
        return S;
    }

    /// Extent of the vector of subspace, ie the [0:last) range
    template<class Spin>
    mpqc::range extents(const std::vector< Subspace<Spin> > &R) {
        size_t last = 0;
        foreach (auto r, R) {
            last = std::max<size_t>(last, *r.end());
        }
        return mpqc::range(0, last);
    }

    /// Grid of subspaces, represented as blocks of determinants defined by alpha/beta pair,
    /// along with sparsity information.
    struct SubspaceGrid {

        /// Construct empty subspace grid
        SubspaceGrid() {}

        /// Construct subspace grid from alpha and beta subspace vectors and sparsity matrix
        SubspaceGrid(const std::vector< Subspace<Alpha> > &A,
                     const std::vector< Subspace<Beta> > &B,
                     const mpqc::matrix<bool> &mask)
            : alpha_(A), beta_(B), mask_(mask)
        {
            verify(alpha_, this->mask_.rows());
            verify(beta_, this->mask_.cols());
            // count total size
            this->dets_ = 0;
            for (int i = 0; i < mask_.rows(); ++i) {
                for (int j = 0; j < mask_.cols(); ++j) {
                    if (this->allowed(i,j))
                        this->dets_ += this->alpha(i).size()*this->beta(j).size();
                }
            }
            std::vector< Subspace<Any> > a(alpha_.begin(), alpha_.end());
            std::vector< Subspace<Any> > b(beta_.begin(), beta_.end());
            this->symm_ = ((a == b) && (mask_ == mask_.transpose()));
        }

        /// Resolve string index to its alpha space
        Space<Alpha> alpha_string_space(int idx) const {
            return space(idx, this->alpha_);
        }

        /// Resolve string index to its beta space
        Space<Beta> beta_string_space(int idx) const {
            return space(idx, this->beta_);
        }

        /// Returns whenever a subspace alpha/beta block is allowed
        bool allowed(int a, int b) const {
            return mask_(a,b);
        }

        /// Returns all alpha subspaces
        const std::vector< Subspace<Alpha> >& alpha() const {
            return alpha_;
        }

        /// Returns all beta subspaces
        const std::vector< Subspace<Beta> >& beta() const {
            return beta_;
        }

        /// Returns i-th alpha subspace
        Subspace<Alpha> alpha(int i) const {
            return alpha_.at(i);
        }

        /// Returns i-th beta subspace
        Subspace<Beta> beta(int i) const {
            return beta_.at(i);
        }

        /// Returns true if the grid is symmetric
        bool symmetric() const {
            return this->symm_;
        }

        /// Returns number of determinants in the grid
        size_t dets() const {
            return this->dets_;
        }

    private:
        std::vector< Subspace<Alpha> > alpha_;
        std::vector< Subspace<Beta> > beta_;
        mpqc::matrix<bool> mask_;
        size_t dets_;
        bool symm_;

    private:

        /// Verify that the subspaces are contigous and their number is N
        /// @throws MPQC_EXCEPTION
        template<class Spin>
        static void verify(const std::vector< Subspace<Spin> > &V, int N) {
            int begin = 0;
            if (V.size() != N)
                throw MPQC_EXCEPTION("Subspace vector size must be equal to %i", N);
            foreach (auto r, V) {
                if (begin != *r.begin())
                    throw MPQC_EXCEPTION("Subspace vector is not contigous");
                begin = *r.end();
            }
        }

    private:

        /// Resolves string index to a subspace in V
        /// @throws MPQC_EXCEPTION if string index does not belong to any subspace
        template<class Spin>
        Space<Spin> space(int idx, const std::vector< Subspace<Spin> > &V) const {
            foreach (const auto &s, V) {
                //std::cout << "idx=" << idx << " s=" << s.rank() << " " << mpqc::range(s) << std::endl;
                if (*s.begin() <= idx && idx < *s.end()) return s;
            }
            throw MPQC_EXCEPTION("String index %i does not belong to any subspace\n", idx);
            return Space<Spin>(-1);
        } 

    };

    /// @}

} // namespace ci
} // namespace mpqc

#endif // MPQC_CI_SUBSPACE_HPP
