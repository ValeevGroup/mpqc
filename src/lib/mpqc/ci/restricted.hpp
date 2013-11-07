#ifndef MPQC_CI_RESTRICTED_HPP
#define MPQC_CI_RESTRICTED_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/subspace.hpp"
#include "mpqc/math/matrix.hpp"

namespace mpqc {
namespace ci {
    
    /// @addtogroup CI
    /// @{

    /// Restricted CI Functor
    struct Restricted {
        /// Create Restricted CI subspace grid
        template<class Index>
        static SubspaceGrid grid(const CI<Restricted, Index> &ci,
                                 const std::vector< Subspace<Alpha> > &A,
                                 const std::vector< Subspace<Beta> > &B)
        {
            mpqc::matrix<bool> mask(A.size(), B.size());
            for (int j = 0; j < B.size(); ++j) {
                for (int i = 0; i < A.size(); ++i) {
                    mask(i,j) = test(ci, A.at(i), B.at(j));
                }
            }
            return SubspaceGrid(A, B, mask);
        }
        /// tests if the excitation to a is allowed
        template<class Index>
        static bool test(const CI<Restricted, Index> &ci,
                         const String &a)
        {
            return (ci.excitation(a) <= ci.config.rank);
        }
        /// tests if simultaneous excitation to space a and space b is allowed
        template<class Index>
        static bool test(const CI<Restricted, Index> &ci,
                         const Space<Alpha> &a,
                         const Space<Beta> &b)
        {
            return ((a.rank() + b.rank()) <= ci.config.rank);
        }
    };
    
    /// @}

} // namespace ci
} // namespace mpqc

#endif // MPQC_CI_RESTRICTED_HPP
