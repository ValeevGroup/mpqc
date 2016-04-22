/*
 * exception.h
 *
 *  Created on: Apr 21, 2016
 *      Author: evaleev
 */

#ifndef SRC_MPQC_UTIL_MISC_EXCEPTION_HPP_
#define SRC_MPQC_UTIL_MISC_EXCEPTION_HPP_

#include <stdexcept>
#include <string>

namespace mpqc {
  namespace exception {
    struct bad_input : public std::runtime_error {
        bad_input(const char* _what);
        bad_input(const std::string& _what);
        virtual ~bad_input();
    };
  } // namespace exception
} // namespace mpqc


#endif /* SRC_MPQC_UTIL_MISC_EXCEPTION_HPP_ */
