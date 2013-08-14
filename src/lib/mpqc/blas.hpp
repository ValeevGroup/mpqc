#ifndef MPQC_BLAS_HPP
#define MPQC_BLAS_HPP

#ifdef HAVE_CONFIG_H
#include "mpqc_config.h"
#endif

#ifdef HAVE_BLAS

#if (F77_INTEGER_WIDTH == 8)
#define BIND_FORTRAN_INTEGER_8 // or not
#endif
#include <boost/numeric/bindings/blas.hpp>

#include "mpqc/matrix.hpp"
#include <Eigen/Core>
#include <boost/numeric/bindings/eigen/matrix.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

// specialize numeric bindings traits for Eigen::Map
template< class Matrix, typename Id, typename Enable>
struct adaptor< Eigen::Map< Matrix >, Id, Enable >
    : adaptor<Matrix, Id, Enable> {};

// specialize numeric bindings traits for mpqc::matrix
template< typename T, int Options, typename Id, typename Enable>
struct adaptor< mpqc::matrix<T, Options>, Id, Enable >
    : adaptor< typename mpqc::matrix<T, Options>::EigenType, Id, Enable> {};

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif // HAVE_BLAS

namespace mpqc {
namespace blas {

    template<typename T>
    void clear(size_t n, T *x) {
	std::fill_n(x, n, 0);
    }

    template<class A, class B>
    typename A::Scalar dot(const Eigen::MatrixBase<A> &a,
                           const Eigen::MatrixBase<B> &b) {
#ifdef HAVE_BLAS
 	return boost::numeric::bindings::blas::dot(a.derived(), b.derived());
#else
 	return a.dot(b);
#endif // HAVE_BLAS
    }

    template<typename Alpha, class A, class B, typename Beta, class C>
    void gemm(const Alpha &alpha,
              const Eigen::MatrixBase<A> &a,
              const Eigen::MatrixBase<B> &b,
              const Beta &beta,
              Eigen::MatrixBase<C> &c) {
#ifdef HAVE_BLAS
	boost::numeric::bindings::blas::gemm(alpha, a.derived(), b.derived(),
                                             beta, c.derived());
#else
        // N.B. fill (rather than scale) by zero to clear NaN/Inf
	if (beta == Beta(0)) c.fill(Beta(0));
        c = alpha*a*b + beta*c;
#endif // HAVE_BLAS
    }

}
}

#endif /* MPQC_BLAS_HPP */
