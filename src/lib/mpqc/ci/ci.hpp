#ifndef MPQC_CI_CONFIG_HPP
#define MPQC_CI_CONFIG_HPP

#include <util/misc/formio.h>
#include "mpqc/ci/string.hpp"
#include "mpqc/matrix.hpp"
#include "mpqc/file.hpp"
#include "mpqc/array.hpp"

namespace mpqc {
namespace ci {
    
    struct Config {
        size_t core, orbitals;
        size_t alpha, beta;
        size_t level;
        size_t roots;
        size_t max;
        size_t collapse;
        double e_ref;
        double convergence;
        double cutoff;
        size_t block, block2;
        Config() {
            core = 0;
            orbitals = 0;
            alpha = 0;
            beta = 0;
            level = 0;
            roots = 1;
            max = 10;
            collapse = 0;
            e_ref = 0;
            convergence = 1e-10;
            cutoff = convergence;
            block = 128;
            block2 = 128;
        }
    };

    template<class Type = void, class Index = ci::String::Index>
    struct CI;

    template<>
    struct CI<void> : Config {

        typedef ci::String::Index Index;

        struct IO : boost::noncopyable {
            range local;
            File::Dataset<double> b, Hb;
            IO(MPI::Comm comm, File::Group io,
               size_t alpha, size_t beta, size_t N) {
                std::vector<range> extents(3);

                extents[0] = range(0, alpha);
                extents[1] = range(0, beta);
                extents[2] = range(0, N);

                extents[1] = extents[1].split2(comm.size()).at(comm.rank());
                local = extents[1];

                b = File::Dataset<double>(io, "b", extents);
                Hb = File::Dataset<double>(io, "Hb", extents);
            }
        };

        CI(const Config &config,
           MPI::Comm comm, File::Group io,
           const ci::String::List<Index> &alpha,
           const ci::String::List<Index> &beta)
            : Config(config),
              comm(comm),
              io(comm, io, alpha.size(), beta.size(), config.max+1),
              alpha(alpha),
              beta(beta)
        {
            dims.push_back(alpha.size());
            dims.push_back(beta.size());

            sc::ExEnv::out0() << sc::indent
                << sc::scprintf("CI space = (%lu %lu):  alpha = %lu, beta = %lu, dets = %lu\n",
                                config.orbitals, config.alpha+config.beta,
                                alpha.size(), beta.size(),
                                alpha.size()*beta.size());

            initialize();

        }

        //void guess(Array<double> C) const {}

        bool test(const String &ex) const {
            return true;
        }

        void initialize() {
            Vector one(1);
            one[0] = 1;
            if (*io.local.begin() <= 0 && 0 < *io.local.end())
                io.b(0,0,0) << one;
            comm.barrier();
        }

    public:
        ci::String::List<ci::String::Index> alpha, beta;
        MPI::Comm comm;
        IO io;
        std::vector<size_t> dims;
    };

    struct Full;
    struct Truncated;

    template<>
    struct CI<Full> : CI<> {
        CI(const Config &config, MPI::Comm comm, File::Group io)
            : CI<>(config, comm, io,
                   ci::strings(config.orbitals, config.alpha),
                   ci::strings(config.orbitals, config.beta))
        {}
    };

    template<>
    struct CI<Truncated> : CI<> {
        CI(const Config &config, MPI::Comm comm, File::Group io)
            : CI<>(config, comm, io,
                   ci::strings(config.orbitals, config.alpha, config.level),
                   ci::strings(config.orbitals, config.alpha, config.level)),
              ref_(CI<>::alpha[0])
        {}
        bool test(const String &ex) const {
            //return 0;
            if (String::difference(ex, ref_) > this->level) return false;
            return true;
        }
    private:
        String ref_;
    };

} // namespace ci
} // namespace mpqc


namespace mpqc {
namespace ci {

}
}

#endif // MPQC_CI_CONFIG_HPP
