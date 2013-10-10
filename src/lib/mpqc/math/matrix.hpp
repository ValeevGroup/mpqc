#ifndef MPQC_MATRIX_HPP
#define MPQC_MATRIX_HPP

#include "mpqc/range.hpp"
#include "math/scmat/matrix.h"

#define EIGEN_DONT_PARALLELIZE

#include <Eigen/Eigen>
// #define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
// #include <Eigen/Sparse>

namespace mpqc {

    /// @addtogroup MathMatrix
    /// @{

    /// Matrix class derived from Eigen::Matrix with additional MPQC integration
    template<typename T, int Order = Eigen::ColMajor>
    struct matrix :
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Order>
    {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Order> EigenType;

	matrix() : EigenType() {}

	/// Construct <i>unititialized</i> matrix
	/// @param m number of rows
	/// @param n number of columns
	/// @warning NOT initialized to zeroes.
	matrix(size_t m, size_t n) : EigenType(m,n) {}

	/// Construct matrix from Eigen type
        template<class A>
	matrix(const Eigen::EigenBase<A> &a) : EigenType(a) {}

	/// Construct matrix from sc::RefSCMatrix matrix
        matrix(sc::RefSCMatrix a) {
            this->resize(a.nrow(), a.ncol());
            apply(assign(), this->rows(), this->cols(), *this, a);
        }

	/// Construct full matrix from sc::RefSCMatrix matrix
        matrix(const sc::RefSymmSCMatrix &a) {
            this->resize(a.n(), a.n());
            apply(assign(), this->rows(), this->cols(), *this, a);
        }

#ifdef DOXYGEN

	/// Operators to access matrix element/block, see @ref MatrixOperators "here".
	Type operator()(i, j);

#else // DOXYGEN

        using EigenType::operator();

	/// Test if both template parameters are integral types
        template<typename T_, typename U_>
        struct is_index : boost::mpl::and_<
            boost::is_integral<T_>,
            boost::is_integral<U_> > {};

	/// Access a matrix block, see @ref MatrixDetail "details".
        template<class Ri, class Rj>
        typename boost::disable_if<
            is_index<Ri, Rj>,
            Eigen::Block<EigenType>
            >::type 
        operator()(const Ri &i, const Rj &j) {
            range ri = range_cast(i);
            range rj = range_cast(j);
            // printf("i=(%i,%i),j=(%i:%i)\n",
            //        *ri.begin(), *ri.end(),
            //        *rj.begin(), *rj.end());
	    return this->block(*ri.begin(), *rj.begin(), ri.size(), rj.size());
        }

        template<class Ri, class Rj>
        typename boost::disable_if<
            is_index<Ri, Rj>,
            Eigen::Block<const EigenType>
            >::type 
        operator()(const Ri &i, const Rj &j) const {
            range ri = range_cast(i);
            range rj = range_cast(j);
            // printf("i=(%i,%i),j=(%i:%i)\n",
            //        *i.begin(), *i.end(),
            //        *j.begin(), *j.end());
	    return this->block(*ri.begin(), *rj.begin(), ri.size(), rj.size());
        }

#endif // DOXYGEN

        void reshape(int m, int n);

    private:
        template<class F, class A, class B>
        static void apply(const F &f, size_t m, size_t n, A &a, B &b) {
	    for (int j = 0; j < n; ++j) {
		for (int i = 0; i < m; ++i) {
		    f(a(i, j), b(i, j));
		}
	    }
        }
        struct assign {
            template<class A, class B>
            void operator()(A &a, const B &b) const {
		a = b;
            }
        };
    };

    /// Vector class derived from Eigen::Matrix with additional MPQC integration
    /// @tparam T vector type
    template<typename T>
    struct vector: Eigen::Matrix<T, Eigen::Dynamic, 1> {

	/// Eigen base type.
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> EigenType;

	/// Construct <i>unititialized</i> vector
	/// @param m vector size
	/// @warning NOT initialized to zeroes
	explicit vector(size_t m = 0) : EigenType(m) {}

	/// Construct vector from Eigen type
        template<class A>
	vector(const Eigen::EigenBase<A> &a) : EigenType(a) {}

	/// Construct vector from iterator range
        vector(const T *begin, const T *end) {
	    EigenType::resize(end - begin, 1);
	    std::copy(begin, end, EigenType::data());
        }

	/// Element access operators, inherited from Eigen base.
        using EigenType::operator();

	/// range operator.
        Eigen::Block<EigenType> operator()(range i) {
	    return this->block(*i.begin(), 0, i.size(), 1);
        }

	/// const range operator.
        Eigen::Block<const EigenType> operator()(range i) const {
	    return this->block(*i.begin(), 0, i.size(), 1);
        }

    };


    /// Convience double matrix type
    typedef matrix<double> Matrix;
    /// Convience double vector type
    typedef vector<double> Vector;

    // typedef Eigen::SparseMatrix<double> Sparse;

    /// absolute max of an Eigen type
    /// @todo Refactor to accept EigenBase types and return proper type
    template<class E>
    double absmax(const E &e) {
        return std::max(fabs(e.maxCoeff()), fabs(e.minCoeff()));
    }

    /// element-wise dot product of two matrices
    /// @todo Refactor to work with EigenBase types
    template<class T>
    T dot(const matrix<T> &a, const matrix<T> &b) {
	T q = 0;
	assert(a.cols() == b.cols());
	for (size_t j = 0; j < a.cols(); ++j) {
	    q += a.col(j).dot(b.col(j));
	}
	return q;
    }

    /// Computes (Eigen::SelfAdjointEigenSolver) eigensystem of a matrix.
    /// Matrix must be symmetric.
    /// @todo Find a better name
    template<class T>
    Eigen::SelfAdjointEigenSolver<Matrix::EigenType> symmetric(const matrix<T> &a) {
	Eigen::SelfAdjointEigenSolver<Matrix::EigenType> es(a);
	if (es.info() != Eigen::Success)
	    throw std::runtime_error("Eigen solver failed");
	return es;
    }

    /// Matrix norm
    /// @todo Refactor to work with EigenBase types
    template<class T>
    T norm(const matrix<T> &a) {
	return a.norm();
    }

    /// Normalize matrix
    /// @todo Refactor to work with EigenBase types
    template<class T>
    void normalize(matrix<T> &a) {
	a *= 1/a.norm();
    }

    /// orthormalize matrix d wrt to *normalized* matrix b
    /// d = normalize(d - (<d|b>*b))
    /// @todo Refactor to work with EigenBase types
    template<class T>
    void orthonormalize(matrix<T> &d, const matrix<T> &b) {
	T db = dot(d, b);
	T bb = 1; //dot(b,b); // should be normalized already
	d = d - (db/bb)*b;
	normalize(d);
        //d /= sqrt(d.norm());
    }

    /// @} Matrix

} // namespace mpqc

namespace sc {
    using mpqc::Vector;
    using mpqc::Matrix;
    using mpqc::dot;
}

#endif /* MPQC_MATRIX_HPP */
