#ifndef MPQC_TENSOR_INDEX_HPP
#define MPQC_TENSOR_INDEX_HPP

namespace mpqc {

    template<int Label>
    struct index {
        static const label = Label;
    };


    template<class A, class B>
    struct contraction {
        typedef typename
        fusion::as_vector<
            fusion::result_of::find<
 labels
    };

    template<class A, class B, class C>
    void contract(const A &a, const B &b, C &c) {
        typedef typename contraction<A,B> AB;

    }

}

#endif /* MPQC_TENSOR_INDEX_HPP */
