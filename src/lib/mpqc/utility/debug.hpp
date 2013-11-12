#ifndef MPQC_UTILITY_DEBUG_HPP
#define MPQC_UTILITY_DEBUG_HPP

#include "mpqc/utility/string.hpp"

#include <stdio.h>
#include <stdarg.h>
#include <string>

namespace mpqc {
    
#ifndef DOXYGEN

    void debug(char const *file, long line, const char *fmt = NULL, ...) {
        std::string s = "DEBUG: " + string_cast(file) + ":" + string_cast(line);
        if (!fmt) {
            printf("%s", s.c_str());
        }
        else {
            s = s + ": " + fmt;
            va_list args;
            va_start(args, fmt);
            vprintf(s.c_str(), args);
            va_end(args);
        }
        fflush(stdout);
    }

#endif // DOXYGEN

}

/// @ingroup CoreUtility
/// Debug/print a printf-formatted statement with file:line prepended.
/// @code
/// MPQC_DEBUG(); // print only file:line information
/// MPQC_DEBUG("(%i,%i) = %f\n", i, j, v(i,j));
/// @endcode
#define MPQC_DEBUG(...) mpqc::debug(__FILE__, __LINE__, __VA_ARGS__)

#endif /* MPQC_UTILITY_DEBUG_HPP */
