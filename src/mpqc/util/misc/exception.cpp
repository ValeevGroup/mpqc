/*
 * exception.cpp
 *
 *  Created on: Apr 21, 2016
 *      Author: evaleev
 */

#include <mpqc/util/misc/exception.hpp>

using mpqc::exception::bad_input;

bad_input::bad_input(const char* _what) : std::runtime_error(_what) {}

bad_input::bad_input(const std::string& _what) : std::runtime_error(_what) {}

bad_input::~bad_input() {}
