#ifndef MPQC_TENSOR_HPP
#define MPQC_TENSOR_HPP

#include <assert.h>
#include <algorithm>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/and.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/ref.hpp>
#include <boost/array.hpp>

#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/array.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/io.hpp>

#include "mpqc/range.hpp"

namespace mpqc {

    template<typename T, size_t N>
    class tensor;

}


namespace mpqc {
namespace detail {
namespace tensor {

    template<class T>
    struct is_tensor : boost::mpl::false_ {};

    template<typename T, size_t N>
    struct is_tensor< mpqc::tensor<T, N> > : boost::mpl::true_ {};

    template<typename T, size_t N>
    struct is_tensor< const mpqc::tensor<T, N> > : boost::mpl::true_ {};


    struct assign {
        template<class T, class U, class S>
        void operator()(T &t, U &u, const S &s) const {
            //std::cout << "t" << s << " = " << "u" << s << std::endl;
            t(s) = u(s);
        }
    };


    template<class R, size_t K>
    struct Tie {
        struct Ref : boost::reference_wrapper<R> {
            Ref(R &r) : boost::reference_wrapper<R>(r) {}
        };
        Ref data[K];
        R& operator[](int i) const {
            return data[i];
        }
    };


    template<class F, class T, class U, class S>
    void apply(F f, T &t, U &u, const S &s, boost::mpl::int_<1>) {
        assert(t.dims()[0] == u.dims()[0]);
        int n = t.dims()[0];
        for (int i = 0; i < n; ++i) {
            f(t, u, boost::fusion::push_front(s, boost::cref(i)));
        }
    }

    template<class F, class T, class U, class S, int K>
    void apply(F f, T &t, U &u, const S &s, boost::mpl::int_<K>) {
        const static int J = K-1; // col major traversal
        assert(t.dims()[J] == u.dims()[J]);
        int n = t.dims()[J];
        for (int i = 0; i < n; ++i) {
            apply(f, t, u,
                  boost::fusion::push_front(s, boost::cref(i)),
                  boost::mpl::int_<K-1>());
        }
    }

} // namespace tensor
} // namespace detail
} // namespace mpqc



namespace mpqc {


    template<class F, class T, class U, size_t N>
    void apply(F f, tensor<T,N> t, tensor<U,N> u) {
        detail::tensor::apply(f, t, u,
                              boost::fusion::vector<>(),
                              boost::mpl::int_<N>());
    }

    enum TENSOR_STORAGE_ORDER { TENSOR_STORAGE_UNDEFINED,
                                TENSOR_COLUMN_MAJOR,
                                TENSOR_ROW_MAJOR };

    template<typename T, size_t N>
    class tensor {

    public:

        static const size_t RANK = N;

        tensor(T *data, const size_t (&dims)[N],
               TENSOR_STORAGE_ORDER order = TENSOR_ROW_MAJOR) {
            this->order_ = order;
            this->data_ = data;
            for (int i = 0; i < N; ++i) {
                this->base_[i] = 0;
                this->dims_[i] = dims[i];
            }
            if (order == TENSOR_COLUMN_MAJOR) {
                size_t stride = 1;
                for (int i = 0; i < N; ++i) {
                    this->strides_[i] = stride;
                    stride *= dims[i];
                }
            }
            else {
                size_t stride = 1;
                for (int i = N-1; i >= 0; --i) {
                    this->strides_[i] = stride;
                    stride *= dims[i];
                }
            }
        }

        tensor(T *data,
               const size_t *dims,
               const size_t *strides) {
            this->order_ = TENSOR_STORAGE_UNDEFINED;
            this->data_ = data;
            std::fill(this->base_, this->base_+N, 0);
            std::copy(dims, dims+N, this->dims_);
            std::copy(strides, strides+N, this->strides_);
        }

        template <typename U>
        void operator=(const tensor<U,N> &u) {
            apply(detail::tensor::assign(), *this, u); 
        }

        const T* data() const { return data_; }

        const size_t* dims() const {
            return dims_;
        }


#define MPQC_TENSOR_RANGE_OPERATOR(Z, N, DATA)                                  \
        tensor<T,N> operator()(BOOST_PP_ENUM_PARAMS(N, const range &r)) {       \
            detail::tensor::Tie<const range, N> t =                             \
                {{ BOOST_PP_ENUM_PARAMS(N,r) }};                                \
            return generate< tensor<T,N> >(*this, t);                           \
        }                                                                       \
                                                                                \
        tensor<const T,N> operator()(BOOST_PP_ENUM_PARAMS(N, const range &r)) const { \
            detail::tensor::Tie<const range, N> t =                             \
                {{ BOOST_PP_ENUM_PARAMS(N,r) }};                                \
            return generate< tensor<const T,N> >(*this, t);                     \
        }                                                                       \

        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_RANGE_OPERATOR, nil)



#define MPQC_TENSOR_INDEX_OPERATOR(Z, N, DATA)                                  \
        T& operator()(BOOST_PP_ENUM_PARAMS(N, int i)) {                         \
            detail::tensor::Tie<const int, N> t =                               \
                {{ BOOST_PP_ENUM_PARAMS(N,i) }};                                \
            return this->operator()(t.data);                                    \
        }                                                                       \
                                                                                \
        const T& operator()(BOOST_PP_ENUM_PARAMS(N, int i)) const {             \
            detail::tensor::Tie<const int, N> t =                               \
                {{ BOOST_PP_ENUM_PARAMS(N,i) }};                                \
            return this->operator()(t.data);                                    \
        }                                                                       \

        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_INDEX_OPERATOR, nil)

        template<class Index>
        T& operator()(const Index &idx) {
            return this->data_[this->index(idx)];
        }

        template<class Index>
        const T& operator()(const Index &idx) const {
            return this->data_[this->index(idx)];
        }

    private:

        TENSOR_STORAGE_ORDER order_;
        T *data_;
        size_t base_[N];
        size_t dims_[N];
        size_t strides_[N];

        tensor() {}

        template<class Index>
        size_t index(const Index &idx, boost::mpl::int_<0>) const {
            return 0;
        }

        template<class Index, int K>
        size_t index(const Index &idx, boost::mpl::int_<K>) const {
            static const int J = K-1;
            return ((this->base_[J] + boost::fusion::at_c<J>(idx))*this->strides_[J] +
                    (index(idx, boost::mpl::int_<J>())));
        }

        template<class Index>
        size_t index(const Index &idx) const {
            static_assert(boost::fusion::result_of::size<Index>::value == N, "");
            size_t i = index(boost::fusion::as_vector(idx), boost::mpl::int_<N>());
            //std::cout << idx << "->" << i << std::endl;
            return i;
        }

        template<class U, size_t K>
        static U generate(U &u, detail::tensor::Tie<const range, K> &rt) {
            static_assert(K == RANK, "Tensor rank mismatch in range operator()");
            U t;
            t.data_ = u.data_;
            for (int i = 0; i < N; ++i) {
                const range &r = rt[i];
                t.dims_[i] = r.size();
                t.base_[i] = u.base_[i] + *r.begin();
                t.strides_[i] = u.strides_[i];
                //printf("generate: %lu %lu %lu\n", t.dims_[i], t.base_[i], t.strides_[i]);
            }
            return t;
        }

    };

}


#endif /* MPQC_TENSOR_HPP */
