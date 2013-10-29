#ifndef MPQC_UTILITY_DEBUG_HPP
#define MPQC_UTILITY_DEBUG_HPP

#include "mpqc/utility/string.hpp"

#include <stdio.h>
#include <stdarg.h>
#include <string>

namespace mpqc {
    
    /// @addtogroup Utility
    /// @{

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

    /// @}

}

/// @ingroup Utility
#define MPQC_DEBUG(...) mpqc::debug(__FILE__, __LINE__, __VA_ARGS__)

#endif /* MPQC_UTILITY_DEBUG_HPP */
