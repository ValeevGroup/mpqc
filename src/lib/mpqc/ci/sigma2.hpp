#ifndef MPQC_CI_SIGMA2_HPP
#define MPQC_CI_SIGMA2_HPP

#include "mpqc/ci/string.hpp"
#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/full/ci.hpp"

#include "mpqc/utility/timer.hpp"
#include "mpqc/range.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/omp.hpp"

#include "mpqc/array.hpp"
#include "mpqc/array/functions.hpp"

namespace mpqc {
namespace ci {

    template<class CI, class Spin>
    void sigma12(const CI &ci,
                 const String &I, const Subspace<Spin> &S,
                 const mpqc::Vector &h, const mpqc::Matrix &V,
                 mpqc::Vector &F)
    {
        size_t count = I.count();
        const auto &list = ci.template strings<Spin>();
        int idx = list[I];

        std::vector<int> O, E;
        for (size_t l = 0; l < I.size(); ++l) {
            if (I[l])
                O.push_back(l); // occ. orbs
            if (!I[l]) {
                E.push_back(l); // exc. orbs
            }
        }

        //asm("#andrey");
        for (auto k = O.begin(); k < O.end(); ++k) {

            E.push_back(*k); // k->k

            for (auto l = E.begin(); l < E.end(); ++l) {
                String J = I.swap(*k,*l);

                // out-of-subspace
                if (abs(ci.excitation(J) - S.rank()) > 1) continue;

                double sgn_kl = sgn(I,*k,*l);
                int kl = index(*k,*l);
                int jdx = (ci.test(J) ? list[J] : -1);

                std::swap(*k,*l); // k->l

                bool singles = false;
                if (S.test(jdx)) {
                    singles = true;
                    F(jdx-*S.begin()) += sgn_kl*h(kl);
                    // l->k, i->i
                    foreach (int i, O) {
                        F(jdx-*S.begin()) += 0.5*sgn_kl*V(index(i,i),kl);
                    }
                }

                // l->k, i->j
                for (auto j = E.begin(); j < E.end()-1; ++j) {
                    //if (!ci.test(J) && *j > *l) continue;
                    for (auto i = O.begin(); i < O.end(); ++i) {
                        String K = J.swap(*i,*j);
                        if (ci.excitation(K) != S.rank()) continue;
                        int kdx = list[K];
                        if (S.test(kdx))
                            F(kdx-*S.begin()) += 0.5*sgn_kl*sgn(J,*i,*j)*V(index(*i,*j),kl);
                    }
                }

                // restore original vectors
                std::swap(*k,*l);

            }

            E.pop_back();

        } // k
    }

    template<class CI, class Spin>
    void sigma12(const CI &ci, Subspace<Spin> I, Subspace<Spin> J,
                 const mpqc::Vector &H, const mpqc::Matrix &V,
                 const mpqc::Matrix &C,
                 mpqc::Matrix &S)
    {
        mpqc::Vector F = mpqc::Vector::Zero(J.size());
        for (int i = 0; i < (int)I.size(); ++i) {
            sigma12(ci, ci.template strings<Spin>()[i+*I.begin()], J, H, V, F);
            for (size_t j = 0; j < J.size(); ++j) {
                double f = F(j);
                F(j) = 0;
                if (fabs(f) < 1e-14) continue;
                S.col(i) += f*C.col(j);
            }
        }
    }

}
}

#endif /* MPQC_CI_SIGMA2_HPP */

