// shamelessly borrowed from MADNESS (GPLed)
// written by Justus Calvin (justus.c79@gmail.com)

#ifndef _util_container_stdarray_h
#define _util_container_stdarray_h

#include <mpqc_config.h>

#ifdef HAVE_STD_ARRAY
#include <array>
#else
#include <util/container/stdarray_bits.h>
namespace std {
  using namespace sc::tr1::array;
}
#endif

#endif // _util_container_stdarray_h

