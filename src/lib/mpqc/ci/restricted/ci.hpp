#ifndef MPQC_CI_RESTRICTED_CI_HPP
#define MPQC_CI_RESTRICTED_CI_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"

#include "mpqc/math/matrix.hpp"
#include "mpqc/file.hpp"

namespace mpqc {
namespace ci {

    template<class Index = ci::String::Index>
    struct Restricted {
        ci::String::List<Index> alpha, beta;
        size_t dets;
        struct {
            range determinants;
            range alpha, beta;
        } local;
        explicit Restricted(const Config &config, MPI::Comm comm)
            : alpha(ci::strings(config.orbitals, config.electrons.alpha, config.rank)),
              beta(ci::strings(config.orbitals, config.electrons.beta, config.rank)),
              dets(0), local(),
              rank_(config.rank)
        {
            int N = rank_+1; // number of excitation blocks
            std::vector<int> A = reorder(alpha);
            std::vector<int> B = reorder(beta);
            assert(N == A.size());
            assert(N == B.size());
            mpqc::matrix<int> AB(N,N);
            AB.fill(0);
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < N-j; ++i) {
                    AB(i,j) = A.at(i)*B.at(j);
                    this->dets += AB(i,j);
                }
            }
            std::cout << "AB = \n" << AB << std::endl;
            local.determinants = range(0,this->dets);
        }
        bool test(const String &a, const String &b) const {
            return (String::difference(a,b) <= this->rank_);
        }

    protected:

        int excitation(const String &a) const {
            String r(a.size(), a.count());
            int diff = String::difference(r,a);
            //std::cout << a << "-" << r << "=" << diff << std::endl;
            return diff;
        }

        std::vector<int> reorder(ci::String::List<Index> &S) const {
            // std::cout << "unsorted" << std::endl;
            // foreach (const String &s, S.data()) {
            //     std::cout << s << " . " << excitation(s) << std::endl;
            // }
            S.sort(SortByExcitation(this));
            // std::cout << "sorted" << std::endl;
            // foreach (const String &s, S.data()) {
            //     std::cout << s << std::endl;
            // }
            std::vector<int> X(this->rank_+1);
            int r = 0;
            foreach (const String &s, S) {
                int x = this->excitation(s); // excitation
                //std::cout << s << " rank = " << x << std::endl;
                assert(x <= this->rank_);
                assert(r <= x);
                r = x;
                ++X.at(x);
            }
            return X;
        }

    private:
        int rank_;

    private:
        struct SortByExcitation {
            explicit SortByExcitation(const Restricted *ci) : ci(ci) {}
            bool operator()(const String &a, const String &b) const {
                return (ci->excitation(a) < ci->excitation(b));
            }
        private:
            const Restricted *ci;
        };
    };

    typedef CI< Restricted<> > RestrictedCI;

}
}

#include "mpqc/ci/restricted/sigma.hpp"
#include "mpqc/ci/restricted/vector.hpp"

#endif // MPQC_CI_RESTRICTED_CI_HPP
