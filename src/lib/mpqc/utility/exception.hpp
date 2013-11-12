#ifndef MPQC_UTILITY_EXCEPTION_HPP
#define MPQC_UTILITY_EXCEPTION_HPP

#include <stdio.h>
#include <stdarg.h>

#include <string>
#include <exception>
#include <iostream>
#include <sstream>

#include "mpqc/utility/string.hpp"

namespace mpqc {
    
    /// @addtogroup CoreUtility
    /// @{

    /// MPQC exception class
    struct Exception : std::exception {
        /// Constructs exception
        explicit Exception(const std::string &msg = "")
            : what_(msg) {}
        /// Constructs exceptin with an optional printf-style format and arguments
        explicit Exception(char const *file, long line, const char *fmt = NULL, ...) {
            what_ = what_ + file + ":" + string_cast(line);
            if (fmt) {
                char buffer[1024] = { };
                va_list args;
                va_start(args, fmt);
                vsnprintf(buffer, 1024-1, fmt, args);
                va_end(args);
                what_ = what_ + ": " + buffer;
            }
        }
        ~Exception() throw() {}
        const char* what() const throw() {
            return what_.c_str();
        }
    private:
        std::string what_;
    };

    /// @}

}

/// @ingroup CoreUtility
/// Constructs mpqc::Exception with file, line information
/// and an optional printf-style format and arguments
#define MPQC_EXCEPTION(...) mpqc::Exception(__FILE__, __LINE__, __VA_ARGS__)

#endif /* MPQC_UTILITY_EXCEPTION_HPP */
