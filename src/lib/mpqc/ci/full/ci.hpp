#ifndef MPQC_CI_FULL_HPP
#define MPQC_CI_FULL_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"

namespace mpqc {
namespace ci {    

    template<class Index = ci::String::Index>
    struct Full {
        ci::String::List<Index> alpha, beta;
        size_t dets;
        struct {
            range beta;
            range determinants;
        } local;
        explicit Full(const Config &config, MPI::Comm comm)
            : alpha(ci::strings(config.orbitals, config.electrons.alpha, config.rank)),
              beta(ci::strings(config.orbitals, config.electrons.beta, config.rank))
        {
            this->dets = alpha.size()*beta.size();
            this->local.beta =
                range(beta.size()).split2(comm.size()).at(comm.rank());
            this->local.determinants =
                range(*local.beta.begin()*alpha.size(), *local.beta.end()*alpha.size());
        }
        bool test(const String &a, const String &b) const {
            return true;
        }
    };

    typedef CI< Full<> > FullCI;

} // namespace ci
} // namespace mpqc

#include "mpqc/ci/full/vector.hpp"

#endif // MPQC_CI_FULL_HPP
