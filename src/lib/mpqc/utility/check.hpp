#ifndef MPQC_UTILITY_CHECK_HPP
#define MPQC_UTILITY_CHECK_HPP

#include "mpqc/utility/exception.hpp"

#include <string>
#include <boost/current_function.hpp>

namespace mpqc {

    inline void check_failed(char const *file, long line,
                             char const *function,
                             char const *expr) {
        throw mpqc::Exception(file, line, (std::string(function) + expr).c_str());
    }

}

#define MPQC_CHECK(expr)						\
    ((expr) ? ((void)0) : ::mpqc::check_failed                  \
     (__FILE__, __LINE__, BOOST_CURRENT_FUNCTION, #expr))

#endif /* MPQC_UTILITY_CHECK_HPP */
