#define BOOST_TEST_MODULE BLAS bindings

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "mpqc/math/blas.hpp"
#include "mpqc/math/matrix.hpp"

BOOST_AUTO_TEST_CASE(dot) {

    using namespace mpqc;

    size_t m = 1<<30;

    double tol = 1e-10;
    Vector a = Vector::Random(m);
    Vector b = Vector::Random(m);

    BOOST_TEST_MESSAGE("Testing dot(a,b)");
    {
        double c = blas::dot(a,b);
        double r = a.dot(b);
        BOOST_REQUIRE_SMALL(c-r, tol);
    }

}


BOOST_AUTO_TEST_CASE(gemm) {

    using namespace mpqc;

    size_t m = 3;
    size_t n = 2;
    size_t k = 3;

    double tol = 1e-10;
    double alpha = 2.5;
    double beta = 1.4;
    Matrix a = Matrix::Random(m,k);
    Matrix b = Matrix::Random(k,n);

    BOOST_TEST_MESSAGE("Testing alpha*a*b + beta*c");
    {
        Matrix c = Matrix::Random(m,n);
        Matrix r = alpha*a*b + beta*c;
        blas::gemm(alpha, a, b, beta, c);
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                BOOST_REQUIRE_SMALL(c(i,j)-r(i,j), tol);
            }
        }
    }

    // BOOST_TEST_MESSAGE("Testing alpha*b'*a'' + beta*c");
    // {
    //     Matrix c = Matrix::Random(n,m);
    //     Matrix r = alpha*b.transpose()*a.transpose() + beta*c;
    //     blas::gemm(alpha, b.transpose(), a.transpose(), beta, c);
    //     for (int j = 0; j < c.cols(); ++j) {
    //         for (int i = 0; i < c.rows(); ++i) {
    //             BOOST_CHECK_SMALL(c(i,j)-r(i,j), 0.0);
    //         }
    //     }
    // }

    

}
