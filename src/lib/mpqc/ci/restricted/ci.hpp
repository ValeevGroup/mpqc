#ifndef MPQC_CI_RESTRICTED_CI_HPP
#define MPQC_CI_RESTRICTED_CI_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"

#include "mpqc/math/matrix.hpp"
#include "mpqc/file.hpp"

namespace mpqc {
namespace ci {

    struct Restricted;

    template<class Index>
    struct CI<Restricted, Index> : CI<void, Index> {

        typedef CI<void, Index> Base;

        struct {
            range determinants;
            range alpha, beta;
        } local;

        CI(const Config &config, MPI::Comm comm, File::Group io)
            : CI<void, Index>(config, comm)
        {
            // sort/space strings according to rank (excitation)
            const auto &sa = Base::template sort<Alpha>(Base::alpha);
            const auto &sb = Base::template sort<Beta>(Base::beta);
            mpqc::matrix<bool> allowed(sa.size(), sb.size());
            for (int j = 0; j < sa.size(); ++j) {
                for (int i = 0; i < sb.size(); ++i) {
                    allowed(i,j) = test(sa.at(i), sb.at(j));
                }
            }
            Base::subspace = SubspaceGrid(sa, sb, allowed);
            local.determinants = range(0, Base::dets());
            local.alpha = range(0, Base::alpha.size());
            local.beta = range(0, Base::beta.size());
            Base::initialize(io, local.determinants);
            Base::summary();
        }

        using Base::test;

        /// tests if the excitation to a is allowed
        bool test(const String &a) const {
            return (Base::excitation(a) <= Base::rank);
        }

        /// tests if simultaneous excitation to a and b is allowed
        bool test(const String &a, const String &b) const {
            // std::cout << "rank = " << this->rank_
            //           << " x(a) " << excitation(a)
            //           << " x(b) " << excitation(b) << std::endl;
            return ((Base::excitation(a) + Base::excitation(b)) <= Base::rank);
        }

        /// tests if simultaneous excitation to space a and string b is allowed
        bool test(const Space<Alpha> &a, const Space<Beta> &b) const {
            return ((a.rank() + b.rank()) <= Base::rank);
        }

    };

}
}

#include "mpqc/ci/restricted/vector.hpp"

#endif // MPQC_CI_RESTRICTED_CI_HPP
