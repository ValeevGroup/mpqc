#ifndef MPQC_BLAS_HPP
#define MPQC_BLAS_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/storage_adaptors.hpp>

// #include <boost/numeric/bindings/ublas/vector.hpp>
// #include <boost/numeric/bindings/ublas/matrix.hpp>
// #include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
// #include <boost/numeric/bindings/ublas/matrix_expression.hpp>

#ifdef HAVE_MKL
#define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
#define BOOST_NUMERIC_BINDINGS_BLAS_MKL
#endif

#ifdef HAVE_CBLAS
#define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
#endif

#ifdef HAVE_BLAS
#define BIND_FORTRAN_INTEGER_8 // or not
#include <boost/numeric/bindings/blas.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mpqc {

namespace blas {
namespace detail {

    template<typename Alpha, class A, class B, typename Beta, class C>
    void prod(const Alpha &alpha, const A &a, const B &b, const Beta &beta, C &c) {
	BOOST_AUTO(ab, alpha*boost::numeric::ublas::prod(a,b));
	if (beta == Beta(0)) {
	    c.assign(ab);
	}
	else {
	    c *= beta;
	    c.plus_assign(ab);
	}
    }

    template<typename Alpha, class A, class B, class C>
    void outer_prod(const Alpha &alpha, const A &a, const B &b, C &c) {
	c.plus_assign(alpha*boost::numeric::ublas::outer_prod(a,b));
    }

}
}

namespace blas {

    struct Context {};

    template<typename T>
    void clear(size_t n, T *x) {
	std::fill_n(x, n, 0);
    }

    template<class A, class B>
    typename A::value_type dot(const A &a, const B &b) {
	return inner_prod(a,b);
// #ifdef HAVE_BLAS
// 	boost::numeric::bindings::blas::axpy(alpha, a, b);
// #else
// 	b.plus_assign(alpha*a);
// #endif // HAVE_BLAS
    }

    template<typename Alpha, class A, class B>
    void axpy(const Alpha &alpha, const A &a, B &b,
	     blas::Context context = blas::Context()) {
#ifdef HAVE_BLAS
	boost::numeric::bindings::blas::axpy(alpha, a, b);
#else
	b.plus_assign(alpha*a);
#endif // HAVE_BLAS
    }


    template<typename Alpha, class A, class B, class C>
    void ger(const Alpha &alpha, const A &a, const B &b, C &c) {
#ifdef HAVE_BLAS
	boost::numeric::bindings::blas::ger(alpha, a, b, c);
#else
	detail::outer_prod(alpha, a, b, c);
#endif // HAVE_BLAS
    }


    template<typename Alpha, class A, class B, typename Beta, class C>
    void gemv(const Alpha &alpha, const A &a, const B &b, const Beta &beta, C &c) {
#ifdef HAVE_BLAS
	boost::numeric::bindings::blas::gemv(alpha, a, b, beta, c);
#else
	detail::prod(alpha, a, b, beta, c);
#endif // HAVE_BLAS
    }


    template<typename Alpha, class A, class B, typename Beta, class C>
    void gemm(const Alpha &alpha, const A &a, const B &b, const Beta &beta, C &c) {
#ifdef HAVE_BLAS
	boost::numeric::bindings::blas::gemm(alpha, a, b, beta, c);
#else
	detail::prod(alpha, a, b, beta, c);
#endif // HAVE_BLAS
    }

    inline void set_num_threads(size_t size) {
// #ifdef HAVE_MKL
// 	mkl_set_num_threads(size);
#if defined(_OPENMP)
	omp_set_num_threads(size);
#endif
    }

}
}

#endif // MPQC_BLAS_HPP
