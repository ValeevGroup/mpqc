// This provides C99-like standard integer types.  It is based on boost.org
// code which has been modified for inclusion in the SC Toolkit.

//  (C) Copyright boost.org 1999. Permission to copy, use, modify, sell
//  and distribute this software is granted provided this copyright
//  notice appears in all copies. This software is provided "as is" without
//  express or implied warranty, and with no claim as to its suitability for
//  any purpose.

#ifndef util_misc_scint_h
#define util_misc_scint_h

#include <mpqc_config.h>

#if HAVE_STDINT_H

#include <stdint.h>

namespace sc {

typedef int8_t         sc_int8_t;
typedef int_least8_t   sc_int_least8_t;
typedef int_fast8_t    sc_int_fast8_t;
typedef uint8_t        sc_uint8_t;
typedef uint_least8_t  sc_uint_least8_t;
typedef uint_fast8_t   sc_uint_fast8_t;
                       
typedef int16_t        sc_int16_t;
typedef int_least16_t  sc_int_least16_t;
typedef int_fast16_t   sc_int_fast16_t;
typedef uint16_t       sc_uint16_t;
typedef uint_least16_t sc_uint_least16_t;
typedef uint_fast16_t  sc_uint_fast16_t;
                       
typedef int32_t        sc_int32_t;
typedef int_least32_t  sc_int_least32_t;
typedef int_fast32_t   sc_int_fast32_t;
typedef uint32_t       sc_uint32_t;
typedef uint_least32_t sc_uint_least32_t;
typedef uint_fast32_t  sc_uint_fast32_t;
                       
typedef intmax_t       sc_intmax_t;
typedef uintmax_t      sc_uintmax_t;
typedef int64_t        sc_int64_t;
typedef int_least64_t  sc_int_least64_t;
typedef int_fast64_t   sc_int_fast64_t;
typedef uint64_t       sc_uint64_t;
typedef uint_least64_t sc_uint_least64_t;
typedef uint_fast64_t  sc_uint_fast64_t;

}

#else

//  This is not a complete implementation of the 1999 C Standard stdint.h
//  header; it doesn't supply various macros which are not advisable for use in
//  C++ programs.

#include <limits.h> // implementation artifact; not part of interface

namespace sc {

//  These are fairly safe guesses for some 16-bit, and most 32-bit and 64-bit
//  platforms.  For other systems, they will have to be hand tailored.
//  Because the fast types are assumed to be the same as the undecorated types,
//  it may be possible to hand tailor a more efficient implementation.

//  8-bit types  -------------------------------------------------------------//

# if UCHAR_MAX == 0xff
     typedef signed char     sc_int8_t;
     typedef signed char     sc_int_least8_t;
     typedef signed char     sc_int_fast8_t;
     typedef unsigned char   sc_uint8_t;
     typedef unsigned char   sc_uint_least8_t;
     typedef unsigned char   sc_uint_fast8_t;
# else
#    error defaults not correct; you must hand modify scint.h
# endif

//  16-bit types  ------------------------------------------------------------//

# if USHRT_MAX == 0xffff
     typedef short           sc_int16_t;
     typedef short           sc_int_least16_t;
     typedef short           sc_int_fast16_t;
     typedef unsigned short  sc_uint16_t;
     typedef unsigned short  sc_uint_least16_t;
     typedef unsigned short  sc_uint_fast16_t;
# else
#    error defaults not correct; you must hand modify scint.h
# endif

//  32-bit types  ------------------------------------------------------------//

# if UINT_MAX == 0xffffffff
     typedef int             sc_int32_t;
     typedef int             sc_int_least32_t;
     typedef int             sc_int_fast32_t;
     typedef unsigned int    sc_uint32_t;
     typedef unsigned int    sc_uint_least32_t;
     typedef unsigned int    sc_uint_fast32_t;
# elif ULONG_MAX == 0xffffffff
     typedef long            sc_int32_t;
     typedef long            sc_int_least32_t;
     typedef long            sc_int_fast32_t;
     typedef unsigned long   sc_uint32_t;
     typedef unsigned long   sc_uint_least32_t;
     typedef unsigned long   sc_uint_fast32_t;
# else
#    error defaults not correct; you must hand modify scint.h
# endif

//  64-bit types + intmax_t and uintmax_t  -----------------------------------//

#if defined(ULONGLONG_MAX) && !defined(ULLONG_MAX)
#    define ULLONG_MAX ULONGLONG_MAX
#endif

# ifdef ULLONG_MAX
//#    if ULLONG_MAX == 18446744073709551615 // 2**64 - 1
#    if ULONGLONG_MAX == (0xffffffffffffffffuLL) // uLL reqd for xlC
     typedef long long            sc_intmax_t;
     typedef unsigned long long   sc_uintmax_t;
     typedef long long            sc_int64_t;
     typedef long long            sc_int_least64_t;
     typedef long long            sc_int_fast64_t;
     typedef unsigned long long   sc_uint64_t;
     typedef unsigned long long   sc_uint_least64_t;
     typedef unsigned long long   sc_uint_fast64_t;
#    else
#       error defaults not correct; you must hand modify scint.h
#    endif
# elif ULONG_MAX != 0xffffffff

#    if ULONG_MAX == 18446744073709551615 // 2**64 - 1
     typedef long                 sc_intmax_t;
     typedef unsigned long        sc_uintmax_t;
     typedef long                 sc_int64_t;
     typedef long                 sc_int_least64_t;
     typedef long                 sc_int_fast64_t;
     typedef unsigned long        sc_uint64_t;
     typedef unsigned long        sc_uint_least64_t;
     typedef unsigned long        sc_uint_fast64_t;
#    else
#       error defaults not correct; you must hand modify scint.h
#    endif
# else // assume no 64-bit integers
#    error 64 bit integer types are required
     typedef sc_int32_t              sc_intmax_t;
     typedef sc_uint32_t             sc_uintmax_t;
# endif

}

#endif

#endif
