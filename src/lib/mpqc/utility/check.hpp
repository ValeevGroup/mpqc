#ifndef MPQC_CHECK_HPP
#define MPQC_CHECK_HPP

#include <stdio.h>
#include <stdarg.h>

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <boost/current_function.hpp>

namespace mpqc {

    namespace {
        inline std::string str(char const *file, long line,
			       char const *arg0 = NULL,
			       char const *arg1 = NULL) {
	    std::stringstream ss;
	    ss << file << ":" << line;
	    if (arg0) ss << ": " << arg0;
	    if (arg1) ss << ": " << arg1;
	    return ss.str();
	}
    }

    void check_failed(char const *file, long line,
                      char const *function,
                      char const *expr) {
        throw std::runtime_error(str(file, line, function, expr));
    };

}

#define CHECK(expr)						\
    ((expr) ? ((void)0) : ::mpqc::check_failed                  \
     (__FILE__, __LINE__, BOOST_CURRENT_FUNCTION, #expr))

#endif /* MPQC_CHECK_HPP */
