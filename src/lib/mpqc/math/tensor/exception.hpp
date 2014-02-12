#ifndef MPQC_MATH_TENSOR_EXCEPTION_HPP
#define MPQC_MATH_TENSOR_EXCEPTION_HPP

#include "mpqc/utility/string.hpp"
#include <util/misc/scexception.h>

namespace mpqc {

    /// @addtogroup MathTensor
    /// @{

    struct TensorDimensionsException : sc::SCException {
        template<typename Rank, typename Dim1, typename Dim2>
        TensorDimensionsException(Rank rank, Dim1 dim1, Dim2 dim2)
            : sc::SCException(("rank<" + string_cast(rank) + "> " +
                                 "dimensions mismatch (" +
                                 string_cast(dim1) + " != " + string_cast(dim2) +
                                 ")").c_str())
        {}
        template<typename Rank1, typename Rank2, typename Dim1, typename Dim2>
        TensorDimensionsException(Rank1 rank1, Rank2 rank2, Dim1 dim1, Dim2 dim2)
            : sc::SCException(("rank<" + string_cast(rank1) + "> and" +
                                 "rank<" + string_cast(rank2) + "> " +
                                 "dimensions mismatch (" +
                                 string_cast(dim1) + " != " + string_cast(dim2) +
                                 ")").c_str())
        {}
    };

    struct TensorIndexException : sc::SCException {
        template<typename Rank, typename Index, typename Begin, typename End>
        TensorIndexException(Rank rank, Index index, Begin begin, End end)
            : sc::SCException(("rank<" + string_cast(rank) + "> " +
                                 "index=" + string_cast(index) + " " +
                                 "outside the dimensions [" +
                                 string_cast(begin) + ":" + string_cast(end) + ")").c_str())
        {}
    };

    struct TensorRangeException : sc::SCException {
        template<typename Rank, typename Range, typename Begin, typename End>
        TensorRangeException(Rank rank, Range range, Begin begin, End end)
            : sc::SCException(("rank<" + string_cast(rank) + "> " +
                                 "range=" + string_cast(range) + " " +
                                 "outside the dimensions [" +
                                 string_cast(begin) + ":" + string_cast(end) + ")").c_str())
        {}
    };

    /// @}

}


#endif /* MPQC_MATH_TENSOR_EXCEPTION_HPP */
