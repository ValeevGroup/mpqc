#ifndef MPQC_MATH_HPP
#define MPQC_MATH_HPP

#include "mpqc/range.hpp"
#include "math/scmat/matrix.h"

#define EIGEN_DONT_PARALLELIZE

#include <Eigen/Eigen>
// #define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
// #include <Eigen/Sparse>

namespace mpqc {

    template<typename T>
    struct matrix :
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigen_base;

	matrix() : eigen_base() {}

	matrix(size_t m, size_t n) : eigen_base(m,n) {}

        template<class A>
	matrix(const Eigen::EigenBase<A> &a) : eigen_base(a) {}

        matrix(const sc::SCMatrix &a) {
            this->resize(a.nrow(), a.ncol());
            apply(assign(), this->rows(), this->cols(), *this, a);
        }

        matrix(const sc::SymmSCMatrix &a) {
            this->resize(a.n(), a.n());
            apply(assign(), this->rows(), this->cols(), *this, a);
        }

        using eigen_base::operator();

        template<typename T_, typename U_>
        struct is_index : boost::mpl::and_<
            boost::is_integral<T_>,
            boost::is_integral<U_> > {};

        template<class Ri, class Rj>
        typename boost::disable_if<
            is_index<Ri, Rj>,
            Eigen::Block<eigen_base>
            >::type 
        operator()(const Ri &i, const Rj &j) {
            range ri = rangify1(i);
            range rj = rangify1(j);
            // printf("i=(%i,%i),j=(%i:%i)\n",
            //        *ri.begin(), *ri.end(),
            //        *rj.begin(), *rj.end());
	    return this->block(*ri.begin(), *rj.begin(), ri.size(), rj.size());
        }

        template<class Ri, class Rj>
        typename boost::disable_if<
            is_index<Ri, Rj>,
            Eigen::Block<const eigen_base>
            >::type 
        operator()(const Ri &i, const Rj &j) const {
            range ri = rangify1(i);
            range rj = rangify1(j);
            // printf("i=(%i,%i),j=(%i:%i)\n",
            //        *i.begin(), *i.end(),
            //        *j.begin(), *j.end());
	    return this->block(*ri.begin(), *rj.begin(), ri.size(), rj.size());
        }

        void reshape(int m, int n);

    private:
        template<class F, class A, class B>
        static void apply(const F &f, size_t m, size_t n, A &a, const B &b) {
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

    template<typename T>
    struct vector: Eigen::Matrix<T, Eigen::Dynamic, 1> {
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> eigen_base;

	explicit vector(size_t m = 0) : eigen_base(m) {}

        template<class A>
	vector(const Eigen::EigenBase<A> &a) : eigen_base(a) {}

        vector(const T *begin, const T *end) {
	    eigen_base::resize(end - begin, 1);
	    std::copy(begin, end, eigen_base::data());
        }

        using eigen_base::operator();

        Eigen::Block<eigen_base> operator()(range i) {
	    return this->block(*i.begin(), 0, i.size(), 1);
        }

        Eigen::Block<const eigen_base> operator()(range i) const {
	    return this->block(*i.begin(), 0, i.size(), 1);
        }

    };


    //typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef matrix<double> Matrix;
    typedef vector<double> Vector;

    // typedef Eigen::SparseMatrix<double> Sparse;

    template<class E>
    double absmax(const E &e) {
        return std::max(fabs(e.maxCoeff()), fabs(e.minCoeff()));
    }

    template<class T>
    T dot(const matrix<T> &a, const matrix<T> &b) {
	T q = 0;
	assert(a.cols() == b.cols());
	for (size_t j = 0; j < a.cols(); ++j) {
	    q += a.col(j).dot(b.col(j));
	}
	return q;
    }

    template<class T>
    Eigen::SelfAdjointEigenSolver<Matrix::eigen_base> symmetric(const matrix<T> &a) {
	Eigen::SelfAdjointEigenSolver<Matrix::eigen_base> es(a);
	if (es.info() != Eigen::Success)
	    throw std::runtime_error("Eigen solver failed");
	return es;
    }

    template<class T>
    T norm(const matrix<T> &a) {
	return a.norm();
    }

    template<class T>
    void normalize(matrix<T> &a) {
	a *= 1/a.norm();
    }

    template<class T>
    void orthonormalize(matrix<T> &d, const matrix<T> &b) {
	T db = dot(d, b);
	T bb = 1; //dot(b,b); // should be normalized already
	d = d - (db/bb)*b;
	normalize(d);
        //d /= sqrt(d.norm());
    }

} // namespace mpqc

namespace sc {
    using mpqc::Vector;
    using mpqc::Matrix;
    using mpqc::dot;
}

#endif /* MPQC_MATH_HPP */
