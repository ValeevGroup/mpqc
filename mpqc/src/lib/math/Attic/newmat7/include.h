
#ifndef _math_newmat7_include_h
#define _math_newmat7_include_h

//$$ include.h           include files required by various versions of C++

//#define Glock                         // for Glockenspiel on the PC
//#define ATandT                        // for AT&T C++ on a Sun

//#define SETUP_C_SUBSCRIPTS              // allow element access via A[i][j]


#define TEMPS_DESTROYED_QUICKLY         // for compiler that delete
					// temporaries too quickly

//#define DO_FREE_CHECK                   // check news and deletes balance

#define USING_DOUBLE                    // elements of type double
//#define USING_FLOAT                   // elements of type float

#define Version21                       // version 2.1 or later


#ifdef _MSC_VER                         // Microsoft
   #include <stdlib.h>

   typedef int jmp_buf[9];
   extern "C"
   {
      int __cdecl setjmp(jmp_buf);
      void __cdecl longjmp(jmp_buf, int);
   }

   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
#endif

#ifdef __ZTC__                          // Zortech
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <stream.hpp>
      #define flush ""                  // doesn't have io manipulators
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
#endif

#ifdef __BCPLUSPLUS__                   // Borland
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #define SystemV                   // optional in Borland
      #include <values.h>               // Borland has both float and values
   #endif
   #undef __TURBOC__                    // also defined in Borland
#endif

#ifdef __TURBOC__                       // Turbo
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <iostream.h>
      #include <iomanip.h>
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #define SystemV                   // optional in Turbo
      #include <values.h>
   #endif
#endif

#ifdef ATandT                           // AT&T
#include <stdlib.h>
#ifdef WANT_STREAM
#include <iostream.h>
#include <iomanip.h>
#endif
#ifdef WANT_MATH
#include <math.h>
#define SystemV                         // must use System V on my Sun
#include <values.h>                     //    as float.h is not present
#endif
#endif

#ifdef __GNUG__                         // Gnu C++
   #include <stdlib.h>
   #ifdef WANT_STREAM
      #include <stream.h>               // no iomanip in G++
      #define flush ""
   #endif
   #ifdef WANT_MATH
      #include <math.h>
      #include <float.h>
   #endif
   #ifndef TEMPS_DESTROYED_QUICKLY
      #define TEMPS_DESTROYED_QUICKLY
   #endif
#endif

#ifdef Glock                            // Glockenspiel
   extern "C" { #include <stdlib.h> }
   #ifdef WANT_STREAM
      #include <stream.hxx>
      #include <iomanip.hxx>
   #endif
   #ifdef WANT_MATH
      extern "C" { #include <math.h> }
      extern "C" { #include <float.h> }
   #endif
   #define NO_LONG_NAMES                // very long names don't work
#endif



#ifdef USING_FLOAT                      // set precision type to float
typedef float Real;
typedef double long_Real;
#endif

#ifdef USING_DOUBLE                     // set precision type to double
typedef double Real;
typedef long double long_Real;
#endif

#endif

