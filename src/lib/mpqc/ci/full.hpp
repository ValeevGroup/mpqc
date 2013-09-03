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
        explicit Full(const Config &config)
            : alpha(ci::strings(config.orbitals, config.electrons.alpha, config.rank)),
              beta(ci::strings(config.orbitals, config.electrons.beta, config.rank))
        {
            this->dets = alpha.size()*beta.size();
        }
        bool test(const String &a, const String &b) const {
            return true;
        }
    protected:
        range local(MPI::Comm comm) const {
            range r = range(beta.size()).split2(comm.size()).at(comm.rank());
            return range(*r.begin()*alpha.size(), *r.end()*alpha.size());
        }
    };

    typedef CI< Full<> > FullCI;

} // namespace ci
} // namespace mpqc

#endif // MPQC_CI_FULL_HPP
