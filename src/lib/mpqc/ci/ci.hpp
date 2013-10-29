#ifndef MPQC_CI_CI_HPP
#define MPQC_CI_CI_HPP

#include <util/misc/formio.h>
#include "mpqc/ci/subspace.hpp"
#include "mpqc/ci/string.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/mpi.hpp"
#include "mpqc/file.hpp"
#include "mpqc/utility/check.hpp"
#include "mpqc/utility/exception.hpp"


/// @defgroup CI mpqc.CI
/// Configuration Interaction (CI) implementation 

namespace mpqc {
namespace ci {

    /// @addtogroup CI
    /// @{
    
    /// CI configuration
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

    /// CI class template.
    /// Specific CI cases should specialize this template.
    template<class Type, class Index = String::Index>
    struct CI;

    /// Base CI class, specific CIs (Full, Restricted, etc) should be derived from this base.
    template<class Index>
    struct CI<void, Index> : Config, boost::noncopyable {

    public:

        ci::String::List<Index> alpha; /// Alpha string list
        ci::String::List<Index> beta; /// Beta string list

        // /// alpha/beta subspaces
        // struct {
        //     std::vector< Subspace<Alpha> > alpha;
        //     std::vector< Subspace<Beta> > beta;
        // } space;
        SubspaceGrid subspace;

        MPI::Comm comm;
        
        /// CI vectors b(C), Hb(sigma) file datasets.
        struct IO : boost::noncopyable {
            File::Dataset<double> b, Hb;
        } vector;

    public:

        template<class Spin>
        const ci::String::List<Index>& strings() const {
            if (Spin::value == Alpha::value)
                return alpha;
            if (Spin::value == Beta::value)
                return beta;
            throw MPQC_EXCEPTION("");
        }
        
    public:

        size_t dets() const {
            return subspace.dets();
        }
        
        /// rank/excitation of the string relative to its ground state, ie [11..00]
        int excitation(const String &s) const {
            return String::difference(s, String(s.size(), s.count()));
        }

        /// Allowed Alpha spaces given Beta space b.
        std::vector< Subspace<Alpha> > allowed(const Space<Beta> &b) const {
            std::vector< Subspace<Alpha> > A;
            for (auto a : subspace.alpha()) {
                if (!test(a,b)) continue;
                A.push_back(a);
            }
            return A;
        }

    public:

        virtual ~CI() {}

        /// test if space configuration is allowed
        template<class Spin>
        static int diff(const Space<Spin> &a, const Space<Spin> &b) {
            return abs(a.rank() - b.rank());
        }

        /// test if space configuration is allowed
        virtual bool test(const Space<Alpha> &a, const Space<Beta> &b) const = 0;

    protected:
        
        /// Construct CI base given configuration and communicator.
        /// Strings will be sorted lexicographically AND by excitation
        CI(const Config &config, MPI::Comm comm) 
            : Config(config),
              alpha(ci::strings(config.orbitals, config.electrons.alpha, config.rank)),
              beta(ci::strings(config.orbitals, config.electrons.beta, config.rank)),
              comm(comm)
        {
        }

        /// initialize IO
        void initialize(File::Group io, range local) {
            // can't have size-0 datasets (size-0 has special meaning)
            MPQC_CHECK(local.size() > 1);
            std::vector<range> extents{local, range(0,Config::max)};
            this->vector.b = File::Dataset<double>(io, "b", extents);
            this->vector.Hb = File::Dataset<double>(io, "Hb", extents);
            guess(local);
        }

        /// initialize guess vector
        void guess(range local) {
            Vector one(1);
            one[0] = 1;
            if (*local.begin() <= 0 && *local.end() > 0)
                vector.b(0,0) << one;
            comm.barrier();
        }

        /// prints CI configuration summary
        void summary() const {
            auto &out = sc::ExEnv::out0();
            // general
            out << sc::indent
                << sc::scprintf("CI space = (%lu %lu): ",
                                Config::orbitals, (Config::electrons.alpha+
                                                   Config::electrons.beta))
                << sc::scprintf("alpha = %lu, beta = %lu, dets = %lu\n",
                                alpha.size(), beta.size(), dets())
                << std::endl;
            // alpha
            out << "alpha excitation population" << std::endl;
            for (int i = 0; i < subspace.alpha().size(); ++i) {
                range r = subspace.alpha(i);
                out << sc::scprintf("%i: size = %lu, \t range = ", i, r.size())
                    << r
                    << std::endl;
            }
            out << std::endl;
            // beta
            out << "beta excitation population" << std::endl;
            for (int i = 0; i < subspace.beta().size(); ++i) {
                range r = subspace.beta(i);
                out << sc::scprintf("%i: size = %lu, \t range = ", i, r.size())
                    << r
                    << std::endl;
            }
            out << std::endl;
            // allowed excitation
            out << "allowed excitations" << std::endl;
            for (auto a : subspace.alpha()) {
                out << "| ";
                for (auto b : subspace.beta()) {
                    out << test(a,b) << " ";
                }
                out << "|" << std::endl;
            }
            out << std::endl;
        }

    protected:

        /// sort string list by rank/excitation and return vector of spaces,
        /// where each space corresponds to range of strings of the same rank
        template<class Spin>
        std::vector< Subspace<Spin> > sort(ci::String::List<Index> &S) const {
            S.sort(Sort(this));
            int r = 0;
            std::vector< Subspace<Spin> > R;
            int begin = 0, end = 0;
            for (const auto &s : S) {
                int x = this->excitation(s);
                if (x == r+1) {
                    R.push_back(Subspace<Spin>(Space<Spin>(r), mpqc::range(begin, end)));
                    begin = end;
                    r = x;
                }
                else if (x == r) {
                    /* do nothing */
                }
                else {
                    throw MPQC_EXCEPTION("internal error");
                }
                ++end;
            }
            if (begin != end)
                R.push_back(Subspace<Spin>(Space<Spin>(r), mpqc::range(begin, end)));
            assert(S.size() == *R.back().end());
            return R;
        }

        /// compare-by-rank functor 
        struct Sort {
            explicit Sort(const CI *ci) : ci_(ci) {}
            bool operator()(const String &a, const String &b) const {
                return ci_->excitation(a) < ci_->excitation(b);
            }
        private:
            const CI *ci_;
        };

    protected:

        // /// Returns allowed excitation range in space 2 given excitation of space 1.
        // template<Spin S1, Spin S2>
        // mpqc::range allowed(const Space<S1> &s1, const std::vector< Subspace<S2> > &s2) const {
        //     int last = 0;
        //     for (const Subspace<S2> &s : s2) {
        //         if (!test(s,s1)) continue;
        //         // must be contigous
        //         MPQC_CHECK(last == *s.begin());
        //         last += s.size();
        //     }
        //     return mpqc::range(0, last);
        // }

    };

    /// @}

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
