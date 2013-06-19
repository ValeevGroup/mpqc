#ifndef MPQC_TENSOR_HPP
#define MPQC_TENSOR_HPP

#include "array/permute.hpp"

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/adaptors.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/tuple/tuple.hpp>

namespace mpqc {
namespace tensor {
namespace detail {


    template<class A>
    bool check_bounds(const A &a, int i) {
	return (*a.index_bases() <= i && i < int(*a.index_bases()+*a.shape()));
    }


    template<class A, long I, long J, long K>
    void permute(A &a, boost::mpl::vector_c<int,I,J,K>) {
	BOOST_STATIC_ASSERT(A::dimensionality == 3);
	array::permute<I,J,K>(a);
    }

    template<class A, long I, long J, long K>
    void permute(A &a, boost::mpl::vector_c<int,I,J,K,3>) {
	BOOST_STATIC_ASSERT(A::dimensionality == 4);
	for (size_t i = 0; i < a.shape()[0]; ++i) {
	    BOOST_AUTO(ai, a[i]);
	    permute(ai, boost::mpl::vector_c<int,I,J,K>());
	}
    }

    struct functor {
	struct copy {
	    template<typename T, typename U>
	    void operator()(const T &t, U &u) const {
		u = t;
	    }
	};
    };

    template<int J, int K, int L, class F, class A, class B>
    void apply(const F &f, const A &a, B &b) {
	BOOST_STATIC_ASSERT(A::dimensionality == 3);
	BOOST_STATIC_ASSERT(B::dimensionality == 3);
	size_t shape[] = { a.shape()[2-L], a.shape()[2-K], a.shape()[2-J] };
	// for (int i = 0; i < 4; ++i) {
	//     std::cout << i << std::endl;
	//     std::cout << a.shape()[i] << std::endl;
	//     std::cout << b.shape()[i] << std::endl;
	//     std::cout <<  shape[i] << std::endl;
	//     //BOOST_ASSERT(a.shape()[i] == shape[i]);
	// }
	for (size_t l = 0; l < b.shape()[0]; ++l) {
	    for (size_t k = 0; k < b.shape()[1]; ++k) {
		for (size_t j = 0; j < b.shape()[2]; ++j) {
		    BOOST_AUTO(t, boost::tie(l, k, j));
		    using boost::get;
		    f(a[get<2-L>(t)][get<2-K>(t)][get<2-J>(t)], b[l][k][j]);
		}
	    }
	}
    }

    template<int I, int J, int K, int L,
	     class F, class A, class B>
    void apply(const F &f, const A &a, B &b) {
	BOOST_STATIC_ASSERT(A::dimensionality == 4);
	BOOST_STATIC_ASSERT(B::dimensionality == 4);
	size_t shape[] = { a.shape()[3-L], a.shape()[3-K],
			   a.shape()[3-J], a.shape()[3-I] };
	// for (int i = 0; i < 4; ++i) {
	//     std::cout << i << std::endl;
	//     std::cout << a.shape()[i] << std::endl;
	//     std::cout << b.shape()[i] << std::endl;
	//     std::cout <<  shape[i] << std::endl;
	//     //BOOST_ASSERT(a.shape()[i] == shape[i]);
	// }
	for (size_t l = 0; l < b.shape()[0]; ++l) {
	    for (size_t k = 0; k < b.shape()[1]; ++k) {
		for (size_t j = 0; j < b.shape()[2]; ++j) {
		    for (size_t i = 0; i < b.shape()[3]; ++i) {
			BOOST_AUTO(t, boost::tie(l, k, j, i));
			using boost::get;
			f(a[get<3-L>(t)][get<3-K>(t)][get<3-J>(t)][get<3-I>(t)],
			  b[l][k][j][i]);
		    }
		}
	    }
	}
    }


    template<class A, class B>
    void plus_assign(A &a, const B &b, boost::mpl::vector_c<int,0,2,1>) {
	assert(a.shape()[0] == b.shape()[1]);
	assert(a.shape()[1] == b.shape()[0]);
	assert(a.shape()[2] == b.shape()[2]);
	size_t L = a.shape()[0];
	size_t K = a.shape()[1];
	size_t N = a.shape()[2];
	for (size_t l = 0; l < L; ++l) {
	    for (size_t k = 0; k < K; ++k) {
		BOOST_AUTO(akl, a[l][k]);
		BOOST_AUTO(blk, b[k][l]);
		for (size_t i = 0; i < N; ++i) {
		    akl[i] += blk[i];
		}
	    }
	}
    }

    template<class A, class B>
    void plus_assign(A &a, const B &b, boost::mpl::vector_c<int,1,0,2>) {
	assert(a.shape()[0] == b.shape()[0]);
	assert(a.shape()[1] == b.shape()[2]);
	assert(a.shape()[2] == b.shape()[1]);
	size_t L = a.shape()[0];
	size_t K = a.shape()[1];
	size_t N = a.shape()[2];
	for (size_t l = 0; l < L; ++l) {
	    BOOST_AUTO(al, a[l]);
	    BOOST_AUTO(bl, b[l]);
	    for (size_t k = 0; k < K; ++k) {
		for (size_t i = 0; i < N; ++i) {
		    al[k][i] += bl[i][k];
		}
	    }
	}
    }

    struct functional {
	struct plus_assign {
	    template<class A, class B>
	    void operator()(A &a, const B &b) const {
		a += b;
	    }
	};

	template<class F, class A, class B>
	static void apply(const F &f, A &a, const B &b,
			  boost::mpl::vector_c<int,0,1,2>) {
	    assert(a.shape()[0] == b.shape()[0]);
	    assert(a.shape()[1] == b.shape()[1]);
	    assert(a.shape()[2] == b.shape()[2]);
	    size_t N = (a.shape()[0]*a.shape()[1]*a.shape()[2]);
	    for (size_t i = 0; i < N; ++i) {
		f(a.data()[i], b.data()[i]);
	    }
	}

	template<class F, class A, class B>
	static void apply(const F &f, A &a, const B &b,
			  boost::mpl::vector_c<int,1,0,2>) {
	    assert(a.shape()[0] == b.shape()[0]);
	    assert(a.shape()[1] == b.shape()[2]);
	    assert(a.shape()[2] == b.shape()[1]);
	    size_t L = a.shape()[0];
	    size_t K = a.shape()[1];
	    size_t N = a.shape()[2];
	    for (size_t l = 0; l < L; ++l) {
		BOOST_AUTO(al, a[l]);
		BOOST_AUTO(bl, b[l]);
		for (size_t k = 0; k < K; ++k) {
		    for (size_t i = 0; i < N; ++i) {
			f(al[k][i], bl[i][k]);
		    }
		}
	    }
	}

	template<class F, class A, class B>
	static void apply(const F &f, A &a, const B &b,
			  boost::mpl::vector_c<int,0,2,1>) {
	    assert(a.shape()[0] == b.shape()[1]);
	    assert(a.shape()[1] == b.shape()[0]);
	    assert(a.shape()[2] == b.shape()[2]);
	    size_t L = a.shape()[0];
	    size_t K = a.shape()[1];
	    size_t N = a.shape()[2];
	    for (size_t l = 0; l < L; ++l) {
		for (size_t k = 0; k < K; ++k) {
		    BOOST_AUTO(akl, a[l][k]);
		    BOOST_AUTO(blk, b[k][l]);
		    for (size_t i = 0; i < N; ++i) {
			f(akl[i], blk[i]);
		    }
		}
	    }
	}
    };


    template<class T, size_t M, size_t N, class R>
    T as_matrix(R &a) {
	size_t size1 = 1;
	size_t size2 = 1;
	BOOST_AUTO(shape, a.shape()+M+N);
	for (size_t i = 0; i < M; ++i) {
	    size1 *= *(--shape);
	}
	for (size_t i = M; i < M+N; ++i) {
	    size2 *= *(--shape);
	}
	//std::cout <<  size1 << " "  << size2 << " " << a.origin() <<std::endl;
	boost::numeric::ublas::column_major O;
	return boost::numeric::ublas::make_matrix(size1, size2, a.data(), O);
    }


}
}
}


namespace cc {
namespace tensor {

    template<class T, typename U>
    void fill(T &t, const U &value) {
	std::fill(t.data(), t.data()+t.num_elements(), value);
    }

    template<class A, class B>
    void copy(const A &a, B &b) {
	BOOST_STATIC_ASSERT((A::dimensionality == B::dimensionality));
	for (size_t i = 0; i < A::dimensionality; ++i) {
	    assert(a.shape()[i] == b.shape()[i]);
	}
	std::copy(a.data(), a.data()+a.num_elements(), b.data());
    }

    template<int I, int J, int K, class A, class B>
    void copy(const A &a, B &b) {
	detail::apply<I,J,K>(detail::functor::copy(), a, b);
    }

    template<int I, int J, int K, int L, class A, class B>
    void copy(const A &a, B &b) {
    	detail::apply<I,J,K,L>(detail::functor::copy(), a, b);
    }

    template<int I, int J, int K, class A>
    void permute(A &a) {
	BOOST_STATIC_ASSERT((A::dimensionality == 3));
	cc::tensor::detail::permute(a, boost::mpl::vector_c<int,I,J,K>());
    }

    template<int I, int J, int K, int L, class A>
    void permute(A &a) {
	BOOST_STATIC_ASSERT((A::dimensionality == 4));
	detail::permute(a, boost::mpl::vector_c<int,I,J,K,L>());
    }

    template<size_t I, size_t J, size_t K, class A, class B>
    void plus_assign(A &a, const B &b) {
	using detail::functional;
	functional::apply(functional::plus_assign(), a, b,
			  boost::mpl::vector_c<int,I,J,K>());
    }

    template<class A, class B>
    void plus_assign(A &a, const B &b) {
	plus_assign<0,1,2>(a,b);
    }


} // namespace tensor
} // namespace cc


namespace cc {

    template<size_t N, typename T = double>
    struct tensor_const_reference;

    template<size_t N, typename T = double>
    struct tensor_reference;

    template<class T>
    struct tensor_traits;


    template<size_t N, typename T>
    struct tensor_traits< tensor_const_reference<N,T> > {
	typedef typename boost::mpl::if_c<
	    (N > 1), tensor_const_reference<N-1,T>, const T&
		>::type const_reference;
    };

    template<size_t N, typename T>
    struct tensor_traits< tensor_reference<N,T> > {
	typedef typename tensor_traits<
	    tensor_const_reference<N,T> >::const_reference const_reference;
	typedef typename boost::mpl::if_c<
	    (N > 1), tensor_reference<N-1,T>, T&
		>::type reference;
    };

    template<class R, class A>
    typename boost::enable_if_c<(A::dimensionality == 1), R>::type
    make_tensor(A &a, int i) {
	BOOST_ASSERT(tensor::detail::check_bounds(a,i));
	return *(a.data() + (i-*a.index_bases())*a.strides()[0]);
    }

    template<class R, class A>
    typename boost::enable_if_c<(A::dimensionality > 1), R>::type
    make_tensor(A &a, int i) {
	BOOST_ASSERT(tensor::detail::check_bounds(a,i));
	const size_t N = A::dimensionality;
	boost::array<size_t,N-1> e;
	boost::array<int,N-1> indices;
	for (size_t j = 1; j < N; ++j) {
	    e[j-1] = a.shape()[j];
	    indices[j-1] = a.index_bases()[j];
	}
	R r(a.data() + (i-*a.index_bases())*a.strides()[0], e);
	r.reindex(indices);
	return r;
    }


    template<size_t N, typename T>
    struct tensor_const_reference : boost::const_multi_array_ref<T,N> {
	typedef typename boost::const_multi_array_ref<T,N> array_ref;
	typedef typename tensor_traits<
	    tensor_const_reference>::const_reference const_reference;
	typedef const T* pointer;
	template<class E>
	tensor_const_reference(pointer data, const E &extents) :
	    array_ref(data, extents) {}
	template<class A>
	tensor_const_reference(const A &a) : array_ref(a) {}
	const_reference operator[](int i) const {
	    return make_tensor<const_reference>(*this, i);
	}
	array_ref& array() const { return *this; }
    };

    template<size_t N, typename T>
    struct tensor_reference : boost::multi_array_ref<T,N> {
	typedef typename boost::multi_array_ref<T,N> array_ref;
	typedef typename tensor_traits<
	    tensor_reference>::reference reference;
	typedef typename tensor_traits<
	    tensor_reference>::const_reference const_reference;
	typedef T* pointer;
	template<class E>
	tensor_reference(pointer data, const E &extents) :
	    array_ref(data, extents) {}
	template<class A>
	tensor_reference(const A &a) : array_ref(a) {}
	reference operator[](int i) {
	    return make_tensor<reference>(*this, i);
	}
	const_reference operator[](int i) const {
	    return make_tensor<const_reference>(*this, i);
	}
	operator tensor_const_reference<N,T>() const {
	    return tensor_const_reference<N,T>(static_cast<const array_ref&>(*this));
	}
	array_ref& array() { return *this; }
	template<int I, int J, int K>
	void permute() {
	    tensor::permute<I,J,K>(*this);
	}
	template<int I, int J, int K, class A>
	void assign(const A &a) {
	    tensor::copy<I,J,K>(a, *this);
	}
	template<class A>
	void plus_assign(const A &a) {
	    tensor::plus_assign(*this, a);
	}
	template<int I, int J, int K, class A>
	void plus_assign(const A &a) {
	    tensor::plus_assign<I,J,K>(*this, a);
	}
    };

    template<size_t N, typename T = double>
    struct Tensor : boost::multi_array<T,N> {
	typedef boost::multi_array<T,N> Array;
	typedef typename tensor_traits<
	    tensor_reference<N,T> >::reference reference;
	typedef typename tensor_traits<
	    tensor_reference<N,T> >::const_reference const_reference;
	Tensor() {}
	template<class E>
	Tensor(const E &extents) : Array(extents) {}
	Array& array() { return *this; }
	template<class E>
	void resize(const E &extents) {
	    Array::resize(extents);
	}
	operator tensor_reference<N,T>() {
	    return tensor_reference<N,T>(*this);
	}
	operator tensor_const_reference<N,T>() const {
	    return tensor_const_reference<N,T>(*this);
	}
	reference operator[](int i) {
	    return make_tensor<reference>(*this, i);
	}
	const_reference operator[](int i) const {
	    return make_tensor<const_reference>(*this, i);
	}
    };
    
}

namespace cc {
namespace tensor {

    template<typename T, class L, class A>
    typename boost::numeric::ublas::vector_adaptor<T>::type
    as_vector(boost::numeric::ublas::matrix<T,L,A> &m) {
	namespace ublas = boost::numeric::ublas;
	return ublas::make_vector(m.size1()*m.size2(), m.data().begin());
    }

    template<typename T, class L, class A>
    typename boost::numeric::ublas::vector_adaptor<const T>::type
    as_vector(const boost::numeric::ublas::matrix<T,L,A> &m) {
	namespace ublas = boost::numeric::ublas;
	return ublas::make_vector(m.size1()*m.size2(), m.data().begin());
    }

    template<size_t M, size_t N>
    typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type
    as_matrix(Tensor<M+N> &a) {
	typedef typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type T;
	return detail::as_matrix<T,M,N>(a);
    }

    template<size_t M, size_t N>
    typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type
    as_matrix(tensor_reference<M+N> a) {
	typedef typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type T;
	return detail::as_matrix<T,M,N>(a);
    }

    template<size_t M, size_t N>
    typename boost::numeric::ublas::matrix_adaptor<
	const double, boost::numeric::ublas::column_major>::type
    as_matrix(tensor_const_reference<M+N> a) {
	typedef typename boost::numeric::ublas::matrix_adaptor<
	const double, boost::numeric::ublas::column_major>::type T;
	return detail::as_matrix<T,M,N>(a);
    }


    template<size_t N, typename T>
    typename boost::numeric::ublas::vector_adaptor<T>::type
    as_vector(Tensor<N,T> &a) {
	namespace ublas = boost::numeric::ublas;
	BOOST_AUTO(m, (as_matrix<N-1,1>(a)));
	return as_vector(m);
    }

    template<size_t N, typename T>
    typename boost::numeric::ublas::vector_adaptor<T>::type
    as_vector(tensor_reference<N,T> a) {
	namespace ublas = boost::numeric::ublas;
	BOOST_AUTO(m, (as_matrix<N-1,1>(a)));
	return as_vector(m);
    }



}
} // namespace cc
} // namespace cchem


#endif /* MPQC_TENSOR_HPP */
