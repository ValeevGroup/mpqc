#ifndef MPQC_ASSERT_HPP
#define MPQC_ASSERT_HPP

#include <stdexcept>
#include "mpqc/utility/string.hpp"

namespace mpqc {

    /// throw exception "file:line: 'expr' failed"
    inline void assert_failed(const char *file, int line, const char *expr) {
        throw std::runtime_error(string_cast(file) + ":" +
                                 string_cast(line) + ": " +
                                 "\'" + expr + "\'" + " failed");
    }

}

/// Check if expression evaluates to true or throw exception
#define MPQC_CHECK(expr)                                                        \
    ((expr) ? ((void)0) : ::mpqc::assert_failed(__FILE__, __LINE__, #expr))

#endif /* MPQC_ASSERT_HPP */
