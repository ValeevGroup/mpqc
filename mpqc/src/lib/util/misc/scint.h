// This provides C99-like standard integer types.  It is based on boost.org
// code which has been modified for inclusion in the SC Toolkit.

//  (C) Copyright boost.org 1999. Permission to copy, use, modify, sell
//  and distribute this software is granted provided this copyright
//  notice appears in all copies. This software is provided "as is" without
//  express or implied warranty, and with no claim as to its suitability for
//  any purpose.

#ifndef util_misc_scint_h
#define util_misc_scint_h

#ifdef HAVE_STDINT_H
#include <stdint.h>

#else

//  This is not a complete implementation of the 1999 C Standard stdint.h
//  header; it doesn't supply various macros which are not advisable for use in
//  C++ programs.

#include <limits.h> // implementation artifact; not part of interface

//  These are fairly safe guesses for some 16-bit, and most 32-bit and 64-bit
//  platforms.  For other systems, they will have to be hand tailored.
//  Because the fast types are assumed to be the same as the undecorated types,
//  it may be possible to hand tailor a more efficient implementation.

//  8-bit types  -------------------------------------------------------------//

# if UCHAR_MAX == 0xff
     typedef signed char     int8_t;
     typedef signed char     int_least8_t;
     typedef signed char     int_fast8_t;
     typedef unsigned char   uint8_t;
     typedef unsigned char   uint_least8_t;
     typedef unsigned char   uint_fast8_t;
# else
#    error defaults not correct; you must hand modify scint.h
# endif

//  16-bit types  ------------------------------------------------------------//

# if USHRT_MAX == 0xffff
     typedef short           int16_t;
     typedef short           int_least16_t;
     typedef short           int_fast16_t;
     typedef unsigned short  uint16_t;
     typedef unsigned short  uint_least16_t;
     typedef unsigned short  uint_fast16_t;
# else
#    error defaults not correct; you must hand modify scint.h
# endif

//  32-bit types  ------------------------------------------------------------//

# if UINT_MAX == 0xffffffff
     typedef int             int32_t;
     typedef int             int_least32_t;
     typedef int             int_fast32_t;
     typedef unsigned int    uint32_t;
     typedef unsigned int    uint_least32_t;
     typedef unsigned int    uint_fast32_t;
# elif ULONG_MAX == 0xffffffff
     typedef long            int32_t;
     typedef long            int_least32_t;
     typedef long            int_fast32_t;
     typedef unsigned long   uint32_t;
     typedef unsigned long   uint_least32_t;
     typedef unsigned long   uint_fast32_t;
# else
#    error defaults not correct; you must hand modify scint.h
# endif

//  64-bit types + intmax_t and uintmax_t  -----------------------------------//

#if defined(ULONGLONG_MAX) && !defined(ULLONG_MAX)
#    define ULLONG_MAX ULONGLONG_MAX
#endif

# ifdef ULLONG_MAX
#    if ULLONG_MAX == 18446744073709551615 // 2**64 - 1
     typedef long long            intmax_t;
     typedef unsigned long long   uintmax_t;
     typedef long long            int64_t;
     typedef long long            int_least64_t;
     typedef long long            int_fast64_t;
     typedef unsigned long long   uint64_t;
     typedef unsigned long long   uint_least64_t;
     typedef unsigned long long   uint_fast64_t;
#    else
#       error defaults not correct; you must hand modify scint.h
#    endif
# elif ULONG_MAX != 0xffffffff

#    if ULONG_MAX == 18446744073709551615 // 2**64 - 1
     typedef long                 intmax_t;
     typedef unsigned long        uintmax_t;
     typedef long                 int64_t;
     typedef long                 int_least64_t;
     typedef long                 int_fast64_t;
     typedef unsigned long        uint64_t;
     typedef unsigned long        uint_least64_t;
     typedef unsigned long        uint_fast64_t;
#    else
#       error defaults not correct; you must hand modify scint.h
#    endif
# else // assume no 64-bit integers
#    error 64 bit integer types are required
     typedef int32_t              intmax_t;
     typedef uint32_t             uintmax_t;
# endif

#endif

typedef int_least64_t sc_int_least64_t;
typedef uint_least64_t sc_uint_least64_t;

#endif
