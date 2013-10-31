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
        size_t ms;
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
            ms = 0;
            rank = 0;
            roots = 1;
            max = 10;
            collapse = 0;
            e_ref = 0.0;
            e_core = 0.0;
            convergence = 1e-10;
            cutoff = convergence;
            block = 1024*4;
        }
        void print(std::ostream& o = sc::ExEnv::out0()) const {
            o << sc::indent << "rank       = " << rank << std::endl;
            o << sc::indent << "# of roots = " << roots << std::endl;
        }
    };

    /// CI class template.
    /// Specific CI cases should specialize this template.
    template<class CIFunctor, class Index = String::Index>
    struct CI;

    /// Base CI class, specific CIs (Full, Restricted, etc) should be derived from this base.
    template<class CIFunctor, class Index>
    struct CI : Config, boost::noncopyable {

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

    protected:
        mpqc::range local_; // range of local determinants

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

        std::vector<mpqc::range> local() const {
            return mpqc::range::split(this->local_, Config::block*Config::block);
        }
        
        /// rank/excitation of the string relative to its ground state, ie [11..00]
        int excitation(const String &s) const {
            return String::difference(s, String(s.size(), s.count()));
        }

    public:
        
        /// Construct CI base given configuration and communicator.
        /// Strings will be sorted lexicographically AND by excitation
        CI(const Config &config, MPI::Comm comm, File::Group io) 
            : Config(config),
              alpha(ci::strings(config.orbitals, config.electrons.alpha, config.rank)),
              beta(ci::strings(config.orbitals, config.electrons.beta, config.rank)),
              comm(comm)
        {
            // sort/space strings according to rank (excitation)
            auto sa = sort<Alpha>(this->alpha);
            auto sb = sort<Beta>(this->beta);

            this->subspace = CIFunctor::grid(*this,
                                             split(sa, Config::block),
                                             split(sb, Config::block));
            
            {
                auto &out = sc::ExEnv::out0();
                // general
                out << sc::indent
                    << sc::scprintf("CI space = (%lu %lu): ",
                                    Config::orbitals, (Config::electrons.alpha+
                                                       Config::electrons.beta))
                    << sc::scprintf("alpha = %lu, beta = %lu, dets = %lu\n",
                                    alpha.size(), beta.size(), this->subspace.dets())
                    << std::endl;
                // alpha
                print(out, "Alpha excitations:\n", sb);
                print(out, "Beta excitations:\n", sa);
            }
            
            // I/O
            {
                size_t dets = this->subspace.dets();
                this->local_ = mpqc::range::split2(mpqc::range(0, dets), this->comm.size())
                    .at(this->comm.rank()); // segment belonging to this node
                MPQC_CHECK(this->local_.size() > 0);
                std::vector<range> extents(2);
                extents[0] = this->local_;
                extents[1] = range(0,Config::max);
                this->vector.b = File::Dataset<double>(io, "b", extents);
                this->vector.Hb = File::Dataset<double>(io, "Hb", extents);
            }

            // initial guess
            {
                guess();
            }

        }

        /// test if space configuration is allowed
        template<class Spin>
        static int diff(const Space<Spin> &a, const Space<Spin> &b) {
            return abs(a.rank() - b.rank());
        }

        /// test if string is allowed
        bool test(const String &a) const {
            return CIFunctor::test(*this, a);
        }

        /// test if space configuration is allowed
        bool test(const Space<Alpha> &a, const Space<Beta> &b) const {
            return CIFunctor::test(*this, a, b);
        }


    protected:

        // /// initialize IO
        // void initialize(File::Group io, range local) {
        //     std::cout << "local = " << local << std::endl;
        //     // can't have size-0 datasets (size-0 has special meaning)
        // }

        /// initialize guess vector
        void guess() {
            Vector one(1);
            one[0] = 1;
            if (this->local_.test(0))
                vector.b(0,0) << one;
            comm.barrier();
        }

        /// prints CI configuration summary
        template<class Spin>
        static void print(std::ostream &out, const std::string &header,
                          const std::vector< Subspace<Spin> > &S)
        {
            out << header;
            for (int i = 0; i < S.size(); ++i) {
                range r = S.at(i);
                out << sc::scprintf("%i: size = %lu, \t range = %s\n",
                                    i, r.size(), string_cast(r).c_str());
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
            foreach (const auto &s, S) {
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
