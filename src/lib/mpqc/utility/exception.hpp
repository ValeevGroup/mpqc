#ifndef MPQC_UTILITY_EXCEPTION_HPP
#define MPQC_UTILITY_EXCEPTION_HPP

#include <stdio.h>
#include <stdarg.h>

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
//#include <boost/current_function.hpp>

namespace mpqc {

    struct Exception : std::runtime_error {
        Exception(const char *file, long line, const char *fmt = NULL, ...)
            : std::runtime_error(str(file, line, fmt)) {}
    private:
        static std::string str(char const *file, long line, char const *fmt, ...) {
            std::string str;
            str = str + file + ":" + string_cast(line);
            if (fmt) {
                char buffer[1024] = { };
                va_list args;
                va_start(args, fmt);
                vsnprintf(buffer, 1024-1, fmt, args);
                va_end(args);
                str = str + ": " + buffer;
            }
            return str;
        }
    };

}

#define MPQC_EXCEPTION(...) mpqc::Exception(__FILE__, __LINE__, __VA_ARGS__)

#endif /* MPQC_UTILITY_EXCEPTION_HPP */
