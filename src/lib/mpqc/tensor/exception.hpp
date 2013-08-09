#ifndef MPQC_TENSOR_EXCEPTION_HPP
#define MPQC_TENSOR_EXCEPTION_HPP

#include <stdexcept>
#include "mpqc/utility/string.hpp"

namespace mpqc {

    /// @addtogroup Tensor
    /// @{

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
