#ifndef MPQC_CI_CI_HPP
#define MPQC_CI_CI_HPP

#include <util/misc/formio.h>
#include "mpqc/ci/string.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/file.hpp"
#include "mpqc/array.hpp"

namespace mpqc {
namespace ci {
    
    struct Config {
        size_t core, orbitals;
        struct {
            size_t alpha, beta;
        } electrons; //!< number of electrons of each spin in CI
        size_t rank;
        size_t roots;
        size_t max;
        size_t collapse;
        double e_ref;
        mutable double e_core;
        double convergence;
        double cutoff;
        size_t block, block2;
        Config() {
            core = 0;
            orbitals = 0;
            electrons.alpha = 0;
            electrons.beta = 0;
            rank = 0;
            roots = 1;
            max = 10;
            collapse = 0;
            e_ref = 0.0;
            e_core = 0.0;
            convergence = 1e-10;
            cutoff = convergence;
            block = 128;
            block2 = 128;
        }
        void print(std::ostream& o = sc::ExEnv::out0()) const {
            o << sc::indent << "rank       = " << rank << std::endl;
            o << sc::indent << "# of roots = " << roots << std::endl;
        }
    };

    template<class Base>
    struct CI : Base, Config {

        using Base::alpha;
        using Base::beta;
        using Base::dets;
        using Base::local;

        struct IO : boost::noncopyable {
            File::Dataset<double> b, Hb;
        };

        CI(const Config &config, MPI::Comm comm, File::Group io) 
            : Config(config), Base(config, comm), comm(comm)
        {
            sc::ExEnv::out0()
                << sc::indent
                << sc::scprintf("CI space = (%lu %lu): ",
                                config.orbitals, (config.electrons.alpha+
                                                  config.electrons.beta))
                << sc::scprintf("alpha = %lu, beta = %lu, dets = %lu\n",
                                alpha.size(), beta.size(), dets);

            // can't have size-0 datasets (size-0 has special meaning)
            MPQC_CHECK(this->local.determinants.size() > 1);
            std::vector<range> extents{this->local.determinants, range(0,config.max)};
            this->vector.b = File::Dataset<double>(io, "b", extents);
            this->vector.Hb = File::Dataset<double>(io, "Hb", extents);

            initialize();
            
        }

        bool test(const String &a) const {
            return Base::test(a, String(a.size(), a.count()));
        }

        void initialize() {
            Vector one(1);
            one[0] = 1;
            if (*local.determinants.begin() <= 0 && 0 < *local.determinants.end())
                vector.b(0,0) << one;
            comm.barrier();
        }

    public:
        MPI::Comm comm;
        IO vector;
    };

} // namespace ci
} // namespace mpqc


namespace sc {

    /// writes Config to sc::StateOut
    inline void ToStateOut(const mpqc::ci::Config &a, StateOut &so, int &count) {
        const char* a_cast = reinterpret_cast<const char*>(&a);
        count += so.put(a_cast, sizeof(mpqc::ci::Config));
    }

    /// reads Config from sc::StateIn
    inline void FromStateIn(mpqc::ci::Config &a, StateIn &si, int &count) {
        char* a_cast = reinterpret_cast<char*>(&a);
        count += si.get(a_cast);
    }

}

#endif // MPQC_CI_CI_HPP
