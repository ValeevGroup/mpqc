
#ifndef SRC_MPQC_UTIL_MISC_STRING_H_
#define SRC_MPQC_UTIL_MISC_STRING_H_

#include <string>
#include <sstream>
#include <locale>
//#include <codecvt> // no gcc support
#include <boost/locale/encoding_utf.hpp>

namespace mpqc {
namespace utility {

namespace detail {

template <typename Ch, typename Arg>
void
__concat(std::basic_ostringstream<Ch>& oss, const Ch* separator, Arg&& arg) {
  oss << arg;
}

template <typename Ch, typename Arg, typename ... Args>
void
__concat(std::basic_ostringstream<Ch>& oss, const Ch* separator, Arg&& arg, Args&&... args) {
  oss << arg;
  if (sizeof...(args) != 0)
    oss << separator;
  __concat(oss, separator, std::forward<Args>(args)...);
}

template <typename Ch, typename ... Args>
std::basic_string<Ch>
_concat(const Ch* separator, Args&&... args) {
  std::basic_ostringstream<Ch> oss;
  detail::__concat(oss, separator, std::forward<Args>(args)...);
  return oss.str();
}

static constexpr const char* default_char_separator = "";
static constexpr const wchar_t* default_wchar_separator = L"";

}

/// @brief converts an utf-8 encoded std::string to an utf-8 encoded std::wstring
inline std::wstring to_wstring(const std::string& str_utf8) {
//  std::wstring_convert<std::codecvt_utf8<wchar_t>> myconv;
//  return myconv.from_bytes(str_utf8);
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<wchar_t>(str_utf8.c_str(), str_utf8.c_str() + str_utf8.size() );
}

/// @brief converts an utf-8 encoded char[] to a utf-8 encoded std::wstring
template <size_t N>
std::wstring to_wstring(const char str_utf8[N]) {
//  std::wstring_convert<std::codecvt_utf8<wchar_t>> myconv;
//  return myconv.from_bytes(&str_utf8[0]);
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<wchar_t>(&str_utf8[0], &str_utf8[0] + N);
}

/// @brief dummy converter
inline std::wstring to_wstring(const std::wstring& wstr_utf8) {
  return wstr_utf8;
}

/// @brief converts UTF-8 encoded wstring to UTF-8 encoded string
inline std::string to_string(const std::wstring& wstr_utf8) {
//  using convert_type = std::codecvt_utf8<wchar_t>;
//  std::wstring_convert<convert_type, wchar_t> converter;
//  return converter.to_bytes(wstr_utf8);
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<char>(wstr_utf8.c_str(), wstr_utf8.c_str() + wstr_utf8.size());
}

/// @return a std::string obtained by streaming \c args to a
/// std::ostringstream
template <typename... Args>
std::string concat(Args&&... args) {
  return detail::_concat<char>(detail::default_char_separator, std::forward(args)...);
}

/// @return a std::wstring obtained by streaming \c args to a
/// std::wostringstream
template <typename ... Args>
std::wstring
wconcat(Args&&... args) {
  return detail::_concat<wchar_t>(detail::default_wchar_separator,std::forward<Args>(args)...);
}

/// @return a std::string obtained by streaming \c args to a
/// std::ostringstream, separated by \c separator
template <typename... Args>
std::string concatsep(const char* separator, Args&&... args) {
  return detail::_concat<char>(separator, std::forward<Args>(args)...);
}

/// @return a std::wstring obtained by streaming \c args to a
/// std::wostringstream, separated by \c separator
template <typename ... Args>
std::wstring
wconcatsep(const wchar_t* separator, Args&&... args) {
  return detail::_concat<wchar_t>(separator, std::forward<Args>(args)...);
}

/// shorthand for concatsep , with \c separator=","
/// @return a std::string obtained by streaming \c args to a
/// std::ostringstream
template <typename ... Args>
std::string
concatcm(Args&&... args) {
  return concatsep(",", std::forward<Args>(args)...);
}

}
}

#endif  // SRC_MPQC_UTIL_MISC_STRING_H_
