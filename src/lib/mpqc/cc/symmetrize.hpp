#ifndef CC_SYMMETRIZE_HPP
#define CC_SYMMETRIZE_HPP

#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/mpl/vector_c.hpp>

#include "cc/tensor.hpp"

namespace cchem {

namespace cc {
namespace detail {

    template<typename S>
    struct symmetrize {
	symmetrize(const S &a, const S &b) : a_(a), b_(b) {}
	template<typename T, typename U>
	void operator()(T &a, U &b) const {
	    T t = a;
	    T u = b;
	    a = a_*a + b_*u;
	    b = a_*b + b_*t;
	}
    private:
	S a_, b_;
    };

    template<size_t N>
    struct symmetric;    

    template<>
    struct symmetric<2> {

    private:
	template<typename T, size_t N>
	struct tile {
	    struct range {
		size_t start, finish;
	    };
	    template<class A>
	    void load(range ri, range rj, const A &a) {
		size_t ni = ri.finish - ri.start;
		size_t nj = rj.finish - rj.start;
		for (size_t j = 0; j < nj; ++j) {
		    BOOST_AUTO(aj, a[j+rj.start]);
		    for (size_t i = 0; i < ni; ++i) {
			data[j][i] = aj[i+ri.start];
		    }
		}
	    }
	    template<class A>
	    void store(range ri, range rj, A &a) {
		size_t ni = ri.finish - ri.start;
		size_t nj = rj.finish - rj.start;
		for (size_t j = 0; j < nj; ++j) {
		    BOOST_AUTO(aj, a[j+rj.start]);
		    for (size_t i = 0; i < ni; ++i) {
			//assert(j < N && i < N);
			//assert(i+ri.start < aj.shape()[0]);
			aj[i+ri.start] = data[j][i];
		    }
		}
	    }
	    T data[N][N];
	};

	template<class F, typename T, size_t N>
	static void apply(const F &f, T (&a)[N][N], T (&b)[N][N]) {
	    for (size_t j = 0; j < N; ++j) {
		for (size_t i = 0; i < N; ++i) {
		    f(a[j][i], b[i][j]);
		}
	    }
	}
    public:
	template<class F, class A>
	static void apply(const F &f, A &a) {
	    assert(a.shape()[0] == a.shape()[1]);
	    size_t N = a.shape()[0];
	    // for (size_t j = 0; j < N; ++j) {
	    // 	for (size_t i = 0; i <= j; ++i) {
	    // 	    f(a[j][i], a[i][j]);
	    // 	}
	    // }
	    // return;
	    static const size_t B = 16;
	    typedef tile<typename A::element,B> T;
	    for (size_t j = 0; j < N; j += B) {
		for (size_t i = 0; i <= j; i += B) {
		    typename T::range ri = { i, std::min(i+B, N) };
		    typename T::range rj = { j, std::min(j+B, N) };
		    // std::cout << ri.start << "-" << ri.finish << std::endl;
		    // std::cout << rj.start << "-" << rj.finish << std::endl;
		    T aij;
		    T aji;
		    aij.load(ri, rj, a);
		    aji.load(rj, ri, a);
		    apply(f, aij.data, aji.data);
		    aij.store(ri, rj, a);
		    aji.store(rj, ri, a);			
		}
	    }
	}
    };


    template<>
    struct symmetric<3> {
	template<class F, class A>
	static void apply(const F &f, A &a, boost::mpl::vector_c<int,1,2>) {
	    assert(a.shape()[0] == a.shape()[1]);
	    size_t N = a.shape()[1];
	    size_t K = a.shape()[2];
	    for (size_t j = 0; j < N; ++j) {
		for (size_t i = 0; i <= j; ++i) {
		    BOOST_AUTO(aij, a[i][j]);
		    BOOST_AUTO(aji, a[j][i]);
		    for (size_t k = 0; k < K; ++k) {
			f(aji[k], aij[k]);
		    }
		}
	    }
	}
	template<class F, class A>
	static void apply(const F &f, A &a, boost::mpl::vector_c<int,0,1>) {
	    assert(a.shape()[1] == a.shape()[2]);
	    size_t K = a.shape()[0];
	    for (size_t k = 0; k < K; ++k) {
		BOOST_AUTO(ak, a[k]);
		symmetric<2>::apply(f, ak);
	    }
	}
    };


}
}

namespace cc {

    template<size_t I, size_t J, typename T, size_t N, typename U>
    void symmetrize2(Tensor<N,T> &A, const U &a, const U &b) {
	boost::mpl::vector_c<int,I,J> ij;
	detail::symmetric<N>::apply(detail::symmetrize<U>(a,b), A, ij);
    }

    template<size_t I, size_t J, typename T, size_t N, typename U>
    void symmetrize2(tensor_reference<N,T> &A, const U &a, const U &b) {
	boost::mpl::vector_c<int,I,J> ij;
	detail::symmetric<N>::apply(detail::symmetrize<U>(a,b), A, ij);
    }

    template<size_t I, size_t J, typename T>
    void symmetrize2(T *p, size_t n0, size_t n1, size_t n2) {
	boost::mpl::vector_c<int,I,J> ij;
	boost::multi_array_ref<T,3> a(p, boost::extents[n2][n1][n0]);
	detail::symmetric<3>::apply(detail::symmetrize<T>(2,-1), a, ij);
    }

}

} // namespace cchem

#endif /* CC_SYMMETRIZE_HPP */
