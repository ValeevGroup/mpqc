#ifndef MPQC_CI_HAMILTONIAN_HPP
#define MPQC_CI_HAMILTONIAN_HPP

#include "mpqc/ci/string.hpp"
#include "mpqc/math/matrix.hpp"
#include <boost/foreach.hpp>

// #define MPQC_PROFILE_ENABLE
// #include "mpqc/profile.hpp"


namespace mpqc {
namespace ci {

    double diagonal(const String &alpha,
                    const mpqc::Vector &h,
                    const mpqc::Matrix &V) {
        //MPQC_PROFILE_LINE;
        double q = 0;
        const auto &o = alpha.occ();
	for (auto i = o.begin(); i < o.end(); ++i) {
	    size_t ii = index(*i,*i);
	    q += h(ii);
	    for (auto j = o.begin(); j < i; ++j) {
		size_t jj = index(*j,*j);
		size_t ij = index(*i,*j);
		q += V(ii,jj) - V(ij,ij);
	    }
	}
	//printf("Hd %e\n", q);
	return q;
    }

    template<typename Index>
    void diagonal2(const String::List<Index> &alpha, const String &beta,
                   const mpqc::Matrix &V, mpqc::Vector &d) {
        //MPQC_PROFILE_LINE;
	const auto &b = beta.occ();
        BOOST_FOREACH (auto j, b) {
            size_t jj = index(j,j);
            auto const &Vj = V.col(jj);
            for (int k = 0; k < alpha.size(); ++k) {
                double q = 0;
                const auto &a = alpha[k].occ();
                BOOST_FOREACH (auto i, a) {
                    size_t ii = index(i,i);
                    q += Vj(i);
                }
                d(k) += q;
            }
        }
    }

    double diagonal2(const String &alpha, const String &beta,
                     const mpqc::Matrix &V) {
        //MPQC_PROFILE_LINE;
        double q = 0;
        const auto &a = alpha.occ();
	const auto &b = beta.occ();
        BOOST_FOREACH (auto j, b) {
            size_t jj = index(j,j);
            auto const &Vj = V.col(jj);
            BOOST_FOREACH (auto i, a) {
                size_t ii = index(i,i);
                q += Vj(ii);
            }
        }
        return q;
    }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_HAMILTONIAN_HPP */
