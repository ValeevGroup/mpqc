#ifndef MPQC_CI_RESTRICTED_SIGMA_HPP
#define MPQC_CI_RESTRICTED_SIGMA_HPP

#include "mpqc/ci/restricted/ci.hpp"
#include "mpqc/utility/exception.hpp"

namespace mpqc {
namespace ci {

    template<class Index>
    void sigma(const CI< Restricted<Index> > &ci,
               const mpqc::Vector &h, const mpqc::Matrix &V,
               const Vector< Restricted<Index> > &C,
               Vector< Restricted<Index> > &S) {
        throw MPQC_EXCEPTION("not implemented");
    }

}
}

#endif /* MPQC_CI_RESTRICTED_SIGMA_HPP */
