
#ifdef __GNUC__
#pragma interface
#endif

#ifdef USE_INLINE
#define INLINE inline
#else
#define INLINE
#endif

INLINE QCXDR::QCXDR() {}

INLINE char QCXDR::byte_swap(char c) { return c; }

#if (USHRT_MAX == 0xffff)
INLINE short QCXDR::byte_swap(unsigned short s)
{
  unsigned short r = s << 8;
  r |= (s & 0xff00) >> 8;
  return r;
}
#endif

#if (UINT_MAX == 0xffffffff) /* 4 byte words */
INLINE int QCXDR::byte_swap(unsigned int i)
{
  unsigned int r = i << 24;
  r |= (i & 0xff00) << 8;
  r |= (i & 0xff0000) >> 8;
  r |= (i & 0xff000000) >> 24;
  return r;
}
#endif

#if (ULONG_MAX == 0xffffffff)
INLINE long QCXDR::byte_swap(unsigned long l)
{
  unsigned long r = l << 24;
  r |= (l & 0xff00) << 8;
  r |= (l & 0xff0000) >> 8;
  r |= (l & 0xff000000) >> 24;
  return r;
}

// for now I'll assume that a 4byte ulong means a 4 byte float...

INLINE float QCXDR::byte_swap(float f)
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
INLINE double QCXDR::byte_swap(double d)
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

INLINE void * QCXDR::byte_swap(void *p)
{
  return (void*) byte_swap((unsigned long) p);
}
  
#endif /* 4byte words */


#if BIGENDIAN

INLINE char   QCXDR::translate(char c) { return c; }
INLINE short  QCXDR::translate(unsigned short s) { return s; }
INLINE int    QCXDR::translate(unsigned int i) { return i; }
INLINE long   QCXDR::translate(unsigned long l) { return l; }
INLINE float  QCXDR::translate(float f) { return f; }
INLINE double QCXDR::translate(double d) { return d; }
INLINE void * QCXDR::translate(void * p) { return p; }

INLINE void QCXDR::translate(char*,int) { }
INLINE void QCXDR::translate(unsigned short*,int) { }
INLINE void QCXDR::translate(unsigned int*,int) { }
INLINE void QCXDR::translate(unsigned long*,int) { }
INLINE void QCXDR::translate(float*,int) { }
INLINE void QCXDR::translate(double*,int) { }

#else

INLINE char   QCXDR::translate(char c) { return byte_swap(c); }
INLINE short  QCXDR::translate(unsigned short s) { return byte_swap(s); }
INLINE int    QCXDR::translate(unsigned int i) { return byte_swap(i); }
INLINE long   QCXDR::translate(unsigned long l) { return byte_swap(l); }
INLINE float  QCXDR::translate(float f) { return byte_swap(f); }
INLINE double QCXDR::translate(double d) { return byte_swap(d); }
INLINE void*  QCXDR::translate(void* p) { return (void*) byte_swap(p); }

INLINE void QCXDR::translate(char*c,int n)
{
  if (!c) return;
  for (int k=0; k < n; k++) c[k] = byte_swap(c[k]);
}

INLINE void QCXDR::translate(unsigned short*s,int n)
{
  if (!s) return;
  for (int k=0; k < n; k++) s[k] = byte_swap(s[k]);
}

INLINE void QCXDR::translate(unsigned int*i,int n)
{
  if (!i) return;
  for (int k=0; k < n; k++) i[k] = byte_swap(i[k]);
}

INLINE void QCXDR::translate(unsigned long*l,int n)
{
  if (!l) return;
  for (int k=0; k < n; k++) l[k] = byte_swap(l[k]);
}

INLINE void QCXDR::translate(float*f,int n)
{
  if (!f) return;
  for (int k=0; k < n; k++) f[k] = byte_swap(f[k]);
}

INLINE void QCXDR::translate(double*d,int n)
{
  if (!d) return;
  for (int k=0; k < n; k++) d[k] = byte_swap(d[k]);
}

#endif /* BIGENDIAN */

#undef INLINE
