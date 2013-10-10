#ifndef MPQC_UTILITY_STRING_HPP
#define MPQC_UTILITY_STRING_HPP

#include <string>
#include <boost/lexical_cast.hpp>

namespace mpqc {

    /// @addtogroup CoreUtility
    /// @{
    
    /// cast type T to string
    template <typename T>
    std::string string_cast(const T& value) {
	return boost::lexical_cast<std::string>(value);
    }

    /// @}

}

#endif /* MPQC_UTILITY_STRING_HPP */
