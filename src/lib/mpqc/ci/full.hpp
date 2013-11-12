#ifndef MPQC_CI_FULL_HPP
#define MPQC_CI_FULL_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/subspace.hpp"
#include "mpqc/math/matrix.hpp"

namespace mpqc {
namespace ci {    

    struct Full {
        template<class Index>
        static SubspaceGrid grid(const CI<Full, Index> &ci,
                                 const std::vector< Subspace<Alpha> > &A,
                                 const std::vector< Subspace<Beta> > &B)
        {
            mpqc::matrix<bool> mask = mpqc::matrix<bool>::Constant(A.size(), B.size(), true);
            return SubspaceGrid(A, B, mask);
        }
        /// tests if the excitation to a is allowed
        template<class Index>
        static bool test(const CI<Full, Index> &ci,
                         const String &a)
        {
            return true;
        }
        /// tests if simultaneous excitation to space a and space b is allowed
        template<class Index>
        static bool test(const CI<Full, Index> &ci,
                         const Space<Alpha> &a,
                         const Space<Beta> &b)
        {
            return true;
        }
    };

} // namespace ci
} // namespace mpqc

#endif // MPQC_CI_FULL_HPP
