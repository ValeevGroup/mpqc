
/* xdr.h -- definition of the external data representation classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      February, 1993
 */

/*
 * as I understand sun's xdr format, things are translated to big-endian
 * byte ordering, which is the byte ordering on sparcs.  For ints we've
 * got it made because we can just use the standard inet ntoh and hton
 * macros/functions.  We'll have to write our own float and double
 * conversions
 */

#ifndef _libQC_xdr_h
#define _libQC_xdr_h

extern "C" {
#ifdef __GNUC__
#include <limits.h>  // use the gnu limits.h
#else
#include <sys/limits.h>
#endif
}

#if defined(sparc) || defined(sgi) || defined(rs6000)
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif

class QCXDR 
{
  public:
   // these perform the translation to/from xdr format
   // single element translations...the value of the passed variable is not
   // changed
    char   translate(char);
    short  translate(unsigned short);
    int    translate(unsigned int);
    long   translate(unsigned long);
    float  translate(float);
    double translate(double);
    void*  translate(void*); // for pointers

   // array transformations...the elements of the array are changed
    void   translate(char*,int);
    void   translate(unsigned short*,int);
    void   translate(unsigned int*,int);
    void   translate(unsigned long*,int);
    void   translate(float*,int);
    void   translate(double*,int);

   // functions for actually performing byte swapping
    char   byte_swap(char);
    short  byte_swap(unsigned short);
    int    byte_swap(unsigned int);
    long   byte_swap(unsigned long);
    float  byte_swap(float);
    double byte_swap(double);
    void*  byte_swap(void*); // for pointers
};


inline char QCXDR::byte_swap(char c) { return c; }

#if (USHRT_MAX == 0xffff)
inline short QCXDR::byte_swap(unsigned short s)
{
  unsigned short r = s << 8;
  r |= (s & 0xff00) >> 8;
  return r;
}
#endif

#if (UINT_MAX == 0xffffffff) /* 4 byte words */
inline int QCXDR::byte_swap(unsigned int i)
{
  unsigned int r = i << 24;
  r |= (i & 0xff00) << 8;
  r |= (i & 0xff0000) >> 8;
  r |= (i & 0xff000000) >> 24;
  return r;
}
#endif

#if (ULONG_MAX == 0xffffffff)
inline long QCXDR::byte_swap(unsigned long l)
{
  unsigned long r = l << 24;
  r |= (l & 0xff00) << 8;
  r |= (l & 0xff0000) >> 8;
  r |= (l & 0xff000000) >> 24;
  return r;
}

// for now I'll assume that a 4byte ulong means a 4 byte float...

inline float QCXDR::byte_swap(float f)
{
  union foo {
    unsigned long l;
    float f;
    } fu;
  fu.f=f;
  fu.l=byte_swap(fu.l);
  return fu.f;
}

// and an 8 byte double
inline double QCXDR::byte_swap(double d)
{
  union foo {
    double d;
    unsigned long l[2];
  } din,dout;

  din.d=d;
  dout.l[0]=byte_swap(din.l[1]);
  dout.l[1]=byte_swap(din.l[0]);

  return dout.d;
}

// I'll also assume pointers are the same length as longs

inline void * QCXDR::byte_swap(void *p)
{
  return (void*) byte_swap((unsigned long) p);
}
  
#endif /* 4byte words */


#if BIGENDIAN

inline char   QCXDR::translate(char c) { return c; }
inline short  QCXDR::translate(unsigned short s) { return s; }
inline int    QCXDR::translate(unsigned int i) { return i; }
inline long   QCXDR::translate(unsigned long l) { return l; }
inline float  QCXDR::translate(float f) { return f; }
inline double QCXDR::translate(double d) { return d; }
inline void * QCXDR::translate(void * p) { return p; }

inline void QCXDR::translate(char*,int) { }
inline void QCXDR::translate(unsigned short*,int) { }
inline void QCXDR::translate(unsigned int*,int) { }
inline void QCXDR::translate(unsigned long*,int) { }
inline void QCXDR::translate(float*,int) { }
inline void QCXDR::translate(double*,int) { }

#else

inline char   QCXDR::translate(char c) { return byte_swap(c); }
inline short  QCXDR::translate(unsigned short s) { return byte_swap(s); }
inline int    QCXDR::translate(unsigned int i) { return byte_swap(i); }
inline long   QCXDR::translate(unsigned long l) { return byte_swap(l); }
inline float  QCXDR::translate(float f) { return byte_swap(f); }
inline double QCXDR::translate(double d) { return byte_swap(d); }
inline void*  QCXDR::translate(void* p) { return (void*) byte_swap(p); }

inline void QCXDR::translate(char*c,int n)
{
  if (!c) return;
  for (int k=0; k < n; k++) c[k] = byte_swap(c[k]);
}

inline void QCXDR::translate(unsigned short*s,int n)
{
  if (!s) return;
  for (int k=0; k < n; k++) s[k] = byte_swap(s[k]);
}

inline void QCXDR::translate(unsigned int*i,int n)
{
  if (!i) return;
  for (int k=0; k < n; k++) i[k] = byte_swap(i[k]);
}

inline void QCXDR::translate(unsigned long*l,int n)
{
  if (!l) return;
  for (int k=0; k < n; k++) l[k] = byte_swap(l[k]);
}

inline void QCXDR::translate(float*f,int n)
{
  if (!f) return;
  for (int k=0; k < n; k++) f[k] = byte_swap(f[k]);
}

inline void QCXDR::translate(double*d,int n)
{
  if (!d) return;
  for (int k=0; k < n; k++) d[k] = byte_swap(d[k]);
}

#endif /* BIGENDIAN */

#endif
