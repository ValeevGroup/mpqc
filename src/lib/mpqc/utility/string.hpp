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

    /// cstring-like object to cast a value to <c>const char*<\c> string.
    struct cstring {
        template <typename T>
        explicit cstring(const T& value)
            : str_(string_cast(value)) {}
        operator const char*() const {
            return this->str();
        }
        const char* str() const {
            return str_.c_str();
        }
    private:
        std::string str_;
    };

    /// @}

}

#endif /* MPQC_UTILITY_STRING_HPP */
