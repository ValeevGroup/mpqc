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

    template<class Index = ci::String::Index>
    struct Truncated {
        ci::String::List<Index> alpha, beta;
        size_t dets;
        explicit Truncated(const Config &config)
            : alpha(ci::strings(config.orbitals, config.electrons.alpha, config.rank)),
              beta(ci::strings(config.orbitals, config.electrons.beta, config.rank)),
              rank_(config.rank)
        {
            alpha.sort(SortByExcitation(alpha[0]));
            beta.sort(SortByExcitation(beta[0]));
            // int N = Config::rank+1; // number of excitation blocks
            // std::vector<int> A(N, 0);
            // std::vector<int> B(N, 0);
            // // count # of strings in each alpha/beta excitation
            // // and do some sanity checks
            // {
            //     int r = 0;
            //     foreach (const String &s, this->alpha) {
            //         int x = String::difference(s, alpha[0]); // excitation
            //         std::cout << s << " rank = " << x << std::endl;
            //         assert(x <= Config::rank);
            //         assert(r <= x);
            //         r = x;
            //         ++A.at(x);
            //     }
            // }
            // {
            //     int r = 0;
            //     foreach (const String &s, this->beta) {
            //         int x = String::difference(s, beta[0]);
            //         assert(x <= Config::rank);
            //         assert(r <= x);
            //         r = x;
            //         ++B.at(x);
            //     }
            // }
            // mpqc::matrix<int> AB(N,N);
            // AB.fill(0);
            // for (int j = 0; j < N; ++j) {
            //     for (int i = 0; i < N-j; ++i) {
            //         AB(i,j) = A.at(i)*B.at(j);
            //         this->dets += AB(i,j);
            //     }
            // }
            // std::cout << "AB = \n" << AB << std::endl;            
        }
        bool test(const String &a, const String &b) const {
            return (String::difference(a,b) <= this->rank_);
        }
    protected:
        range local(MPI::Comm) const {
            throw;
        }
    private:
        int rank_;
        struct SortByExcitation {
            explicit SortByExcitation(const String &ref) : ref_(ref) {}
            bool operator()(const String &a, const String &b) const {
                return (String::difference(a, ref_) < String::difference(b, ref_));
            }
        private:
            String ref_;
        };
    };

    template<class Base>
    struct CI : Base, Config {

        using Base::alpha;
        using Base::beta;
        using Base::dets;

        struct IO : boost::noncopyable {
            range local;
            File::Dataset<double> b, Hb;
        };

        CI(const Config &config, MPI::Comm comm, File::Group io) 
            : Config(config), Base(config), comm(comm)
        {
            sc::ExEnv::out0()
                << sc::indent
                << sc::scprintf("CI space = (%lu %lu): ",
                                config.orbitals, (config.electrons.alpha+
                                                  config.electrons.beta))
                << sc::scprintf("alpha = %lu, beta = %lu, dets = %lu\n",
                                alpha.size(), beta.size(), dets);

            this->vector.local = Base::local(comm);
            std::vector<range> extents{this->vector.local, range(0,config.max)};
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
            if (*vector.local.begin() <= 0 && 0 < *vector.local.end())
                vector.b(0,0) << one;
            comm.barrier();
        }

    public:
        MPI::Comm comm;
        IO vector;
    };

    typedef CI< Full<> > FullCI;
    typedef CI< Truncated<> > TruncatedCI;

    // template<>
    // struct CI<Truncated> : CI<> {

    //     CI(const Config &config, MPI::Comm comm, File::Group io)
    //         : CI<>(config, comm, io,
    //                ci::strings(config.orbitals, config.alpha, config.rank),
    //                ci::strings(config.orbitals, config.alpha, config.rank)),
    //           ref_(CI<>::alpha[0])
    //     {
    //         alpha.sort(Sort(this->ref_));
    //         beta.sort(Sort(this->ref_));
    //         int N = config.rank+1; // number of excitation blocks
    //         std::vector<int> A(N, 0);
    //         std::vector<int> B(N, 0);
    //         // count # of strings in each alpha/beta excitation
    //         // and do some sanity checks
    //         {
    //             int r;
    //             r = 0;
    //             foreach (const String &s, this->alpha) {
    //                 int x = String::difference(s, ref_); // excitation
    //                 ++A.at(x);
    //                 //std::cout << s << " rank = " << x << std::endl;
    //                 assert(r <= x);
    //                 r = x;
    //             }
    //             r = 0;
    //             foreach (const String &s, this->beta) {
    //                 int x = String::difference(s, ref_);
    //                 ++B.at(x);                
    //                 assert(r <= x);
    //                 r = x;
    //             }
    //         }
    //         size_t determinants = 0;
    //         mpqc::matrix<int> AB(N,N);
    //         AB.fill(0);
    //         for (int j = 0; j < N; ++j) {
    //             for (int i = 0; i < N-j; ++i) {
    //                 AB(i,j) = A.at(i)*B.at(j);
    //                 determinants += AB(i,j);
    //             }
    //         }
    //         std::cout << "AB = \n" << AB << std::endl;
    //         std::cout << "# Dets: " << determinants << std::endl;
    //     }

    //     bool test(const String &ex) const {
    //         //return 0;
    //         if (String::difference(ex, ref_) > this->rank) return false;
    //         return true;
    //     }
    // private:
    //     String ref_;
    //     struct Sort {
    //         explicit Sort(const String &ref) : ref_(ref) {}
    //         bool operator()(const String &a, const String &b) const {
    //             return (String::difference(a, ref_) < String::difference(b, ref_));
    //         }
    //     private:
    //         String ref_;
    //     };
    // };

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

#endif // MPQC_CI_CONFIG_HPP
