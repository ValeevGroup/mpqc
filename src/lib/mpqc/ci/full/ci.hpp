#ifndef MPQC_CI_FULL_HPP
#define MPQC_CI_FULL_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/utility/exception.hpp"

namespace mpqc {
namespace ci {    

    struct Full;


    template<class Index>
    struct CI<Full, Index> : CI<void, Index> {

        typedef CI<void, Index> Base;

        struct {
            range determinants;
            range alpha, beta;
        } local;

        CI(const Config &config, MPI::Comm comm, File::Group io)
            : CI<void, Index>(config, comm)
        {
            auto &alpha = Base::alpha;
            auto &beta = Base::beta;

            // // sort/space strings according to rank (excitation)
            // Base::template sort<Alpha>(alpha);
            // Base::template sort<Beta>(beta);
            // mpqc::matrix<bool> allowed = mpqc::matrix<bool>::Constant(1, 1, true);
            // std::vector< Subspace<Alpha> >
            //     sa{ Subspace<Alpha>(Space<Alpha>(0), mpqc::range(0, alpha.size())) };
            // std::vector< Subspace<Beta> >
            //     sb{ Subspace<Beta>(Space<Beta>(0), mpqc::range(0, beta.size())) };

            // sort/space strings according to rank (excitation)
            const auto &sa = Base::template sort<Alpha>(alpha);
            const auto &sb = Base::template sort<Beta>(beta);
            mpqc::matrix<bool> allowed = mpqc::matrix<bool>::Constant(sa.size(), sb.size(), true);

            Base::subspace = SubspaceGrid(sa, sb, allowed);
            this->local.beta =
                range(beta.size()).split2(Base::comm.size()).at(Base::comm.rank());
            this->local.determinants =
                range(*local.beta.begin()*alpha.size(),
                      *local.beta.end()*alpha.size());
            Base::initialize(io, local.determinants);
            Base::summary();
        }

        using Base::test;

        /// tests if the excitation to a is allowed
        bool test(const String &a) const {
            return true;
        }

        /// tests if simultaneous excitation to a and b is allowed
        bool test(const String &a, const String &b) const {
            return true;
        }

        /// tests if simultaneous excitation to space a and space b is allowed
        bool test(const Space<Alpha> &a, const Space<Beta> &b) const {
            return true;
        }

    };

} // namespace ci
} // namespace mpqc

#include "mpqc/ci/full/vector.hpp"

#endif // MPQC_CI_FULL_HPP
