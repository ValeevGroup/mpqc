#ifndef MPQC_STRING_HPP
#define MPQC_STRING_HPP

#include <string>
#include <boost/lexical_cast.hpp>

namespace mpqc {
    
    template <typename T>
    std::string string_cast(const T& value) {
	return boost::lexical_cast<std::string>(value);
    }

}

#endif /* MPQC_STRING_HPP */
