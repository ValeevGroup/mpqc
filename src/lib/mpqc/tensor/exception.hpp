#ifndef MPQC_TENSOR_EXCEPTION_HPP
#define MPQC_TENSOR_EXCEPTION_HPP

#include <stdexcept>
#include "mpqc/utility/string.hpp"

namespace mpqc {

    /// @addtogroup Tensor
    /// @{

    struct TensorDimensionsException : std::runtime_error {
        template<typename Rank, typename Dim1, typename Dim2>
        TensorDimensionsException(Rank rank, Dim1 dim1, Dim2 dim2)
            : std::runtime_error("rank<" + string_cast(rank) + "> " +
                                 "dimensions mismatch (" +
                                 string_cast(dim1) + " != " + string_cast(dim2) +
                                 ")")
        {}
        template<typename Rank1, typename Rank2, typename Dim1, typename Dim2>
        TensorDimensionsException(Rank1 rank1, Rank2 rank2, Dim1 dim1, Dim2 dim2)
            : std::runtime_error("rank<" + string_cast(rank1) + "> and" +
                                 "rank<" + string_cast(rank2) + "> " +
                                 "dimensions mismatch (" +
                                 string_cast(dim1) + " != " + string_cast(dim2) +
                                 ")")
        {}
    };

    struct TensorIndexException : std::runtime_error {
        template<typename Rank, typename Index, typename Begin, typename End>
        TensorIndexException(Rank rank, Index index, Begin begin, End end)
            : std::runtime_error("rank<" + string_cast(rank) + "> " +
                                 "index=" + string_cast(index) + " " +
                                 "outside the dimensions [" +
                                 string_cast(begin) + ":" + string_cast(end) + ")")
        {}
    };

    struct TensorRangeException : std::runtime_error {
        template<typename Rank, typename Range, typename Begin, typename End>
        TensorRangeException(Rank rank, Range range, Begin begin, End end)
            : std::runtime_error("rank<" + string_cast(rank) + "> " +
                                 "range=" + string_cast(range) + " " +
                                 "outside the dimensions [" +
                                 string_cast(begin) + ":" + string_cast(end) + ")")
        {}
    };

    /// @}

}


#endif /* MPQC_TENSOR_EXCEPTION_HPP */
