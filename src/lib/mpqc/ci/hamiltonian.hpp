#ifndef MPQC_CI_HAMILTONIAN_HPP
#define MPQC_CI_HAMILTONIAN_HPP

#include "mpqc/ci/string.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/utility/foreach.hpp"

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
        foreach (auto j, b) {
            size_t jj = index(j,j);
            auto const &Vj = V.col(jj);
            for (int k = 0; k < alpha.size(); ++k) {
                double q = 0;
                const auto &a = alpha[k].occ();
                foreach (auto i, a) {
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
        foreach (auto j, b) {
            size_t jj = index(j,j);
            auto const &Vj = V.col(jj);
            foreach (auto i, a) {
                size_t ii = index(i,i);
                q += Vj(ii);
            }
        }
        return q;
    }

    // template<class Index>
    // void diagonal(const String::List<Index> &alpha,
    // 		  const String &beta,
    // 		  const Vector &h, const Matrix &V,
    // 		  Vector &H) {
    // 	auto &b = beta;
    // 	size_t k = 0;
    // 	for (auto &a : alpha) {
    // 	    double q = 0;
    // 	    for (auto &i : a.occ()) {
    // 	    	size_t ii = index(i,i);
    // 	    	q += h(ii);
    // 	    	for (auto &j : a.occ()) {
    // 	    	    size_t jj = index(j,j);
    // 	    	    size_t ij = index(i,j);
    // 	    	    //q += V(ii,jj) - V(ij,ij);
    // 	    	}
    // 	    	for (auto &j : b.occ()) {
    // 	    	    size_t jj = index(j,j);
    // 	    	    size_t ij = index(i,j);
    // 	    	    //q += V(ii,jj);
    // 	    	}
    // 	    }
    // 	    for (int i : b.occ()) {
    // 	    	size_t ii = index(i,i);
    // 	    	q += h(ii);
    // 	    	for (auto &j : b.occ()) {
    // 	    	    size_t jj = index(j,j);
    // 	    	    size_t ij = index(i,j);
    // 	    	    //q += V(ii,jj) - V(ij,ij);
    // 	    	}
    // 	    }
    // 	    H[k++] = q;
    // 	}
    // }

    // template<class Index>
    // Vector diagonal(const String::List<Index> &alpha,
    // 		    const String &beta,
    // 		    const Vector &h, const Matrix &V) {
    // 	Vector H(alpha.size());
    // 	diagonal(alpha, beta, h, V, H);
    // 	return H;
    // }

} // namespace ci
} // namespace mpqc


#endif /* MPQC_CI_HAMILTONIAN_HPP */
