// shamelessly borrowed from MADNESS (GPLed)
// written by Justus Calvin (justus.c79@gmail.com)

#ifndef _util_container_stdarray_h
#define _util_container_stdarray_h

#include <scconfig.h>

#if defined(SC_USE_ARRAY)
#  include <array>
#elif defined(SC_USE_TR1_ARRAY)
#  include <tr1/array>
#elif defined(SC_USE_BOOST_TR1_ARRAY_HPP)
#  include <boost/tr1/array.hpp>
#else
#  define SC_HAS_STD_ARRAY 1
#  include <util/container/stdarray_bits.h>
   namespace std {
      using namespace sc::tr1::array;
   }
#endif

namespace std {
#if defined(SC_HAS_STD_TR1_ARRAY) && !defined(SC_HAS_STD_ARRAY)
#   define SC_HAS_STD_ARRAY 1
    using ::std::tr1::array;
    using ::std::tr1::swap;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::get;
#endif
}

#endif
