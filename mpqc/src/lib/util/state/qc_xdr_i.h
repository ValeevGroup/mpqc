
#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_INLINE
#define INLINE inline
#else
#define INLINE
#endif

inline void byte_copy2(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  r[0] = d[0];
  r[1] = d[1];
}

inline void byte_copy4(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  r[0] = d[0];
  r[1] = d[1];
  r[2] = d[2];
  r[3] = d[3];
}

inline void byte_copy8(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  r[0] = d[0];
  r[1] = d[1];
  r[2] = d[2];
  r[3] = d[3];
  r[4] = d[4];
  r[5] = d[5];
  r[6] = d[6];
  r[7] = d[7];
}

inline void byte_copy16(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  byte_copy8(d,r);
  byte_copy8(&d[8],&r[8]);
}

inline void byte_swap2(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  r[1] = d[0];
  r[0] = d[1];
}

inline void byte_swap4(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  r[3] = d[0];
  r[2] = d[1];
  r[1] = d[2];
  r[0] = d[3];
}

inline void byte_swap8(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  r[7] = d[0];
  r[6] = d[1];
  r[5] = d[2];
  r[4] = d[3];
  r[3] = d[4];
  r[2] = d[5];
  r[1] = d[6];
  r[0] = d[7];
}

inline void byte_swap16(void*data, void*result)
{
  char*d = (char*) data;
  char*r = (char*) result;
  byte_swap8(d,&r[8]);
  byte_swap8(&d[8],r);
}

INLINE QCXDR::QCXDR() {}

INLINE char QCXDR::byte_swap(char c) { return c; }

INLINE short QCXDR::byte_swap(unsigned short s)
{
  unsigned short r;
  if (sizeof(unsigned short) == 2) byte_swap2(&s, &r);
  else if (sizeof(unsigned short) == 4) byte_swap4(&s, &r);
  else if (sizeof(unsigned short) == 8) byte_swap8(&s, &r);
  else if (sizeof(unsigned short) == 16) byte_swap16(&s, &r);
  return r;
}

INLINE int QCXDR::byte_swap(unsigned int i)
{
  unsigned int r;
  if (sizeof(unsigned int) == 2) byte_swap2(&i, &r);
  else if (sizeof(unsigned int) == 4) byte_swap4(&i, &r);
  else if (sizeof(unsigned int) == 8) byte_swap8(&i, &r);
  else if (sizeof(unsigned int) == 16) byte_swap16(&i, &r);
  return r;
}

INLINE long QCXDR::byte_swap(unsigned long l)
{
  unsigned long r;
  if (sizeof(unsigned long) == 2) byte_swap2(&l, &r);
  else if (sizeof(unsigned long) == 4) byte_swap4(&l, &r);
  else if (sizeof(unsigned long) == 8) byte_swap8(&l, &r);
  else if (sizeof(unsigned long) == 16) byte_swap16(&l, &r);
  return r;
}

// for now I'll assume that a 4byte ulong means a 4 byte float...

INLINE void QCXDR::byte_swap(float*f)
{
#if defined(sgi) && !defined(__GNUG__)
  fprintf(stderr,
          "QCXDR::byte_swap(float): not available because of compiler bug\n");
  abort();
  return;
#else
  char r[16];
  if (sizeof(float) == 2) { byte_swap2(f, r); byte_copy2(r, f); }
  else if (sizeof(float) == 4) { byte_swap4(f, r); byte_copy4(r, f); }
  else if (sizeof(float) == 8) { byte_swap8(f, r); byte_copy8(r, f); }
  else if (sizeof(float) == 16) { byte_swap16(f, r); byte_copy16(r, f); }
#endif
}

// and an 8 byte double
INLINE void QCXDR::byte_swap(double *d)
{
  char r[16];
  if (sizeof(double) == 2) { byte_swap2(d, r); byte_copy2(r, d); }
  else if (sizeof(double) == 4) { byte_swap4(d, r); byte_copy4(r, d); }
  else if (sizeof(double) == 8) { byte_swap8(d, r); byte_copy8(r, d); }
  else if (sizeof(double) == 16) { byte_swap16(d, r); byte_copy16(r, d); }
}

// I'll also assume pointers are the same length as longs

INLINE void * QCXDR::byte_swap(void *p)
{
  return (void*) byte_swap((unsigned long) p);
}
  

#if BIGENDIAN

INLINE char   QCXDR::translate(char c) { return c; }
INLINE short  QCXDR::translate(unsigned short s) { return s; }
INLINE int    QCXDR::translate(unsigned int i) { return i; }
INLINE long   QCXDR::translate(unsigned long l) { return l; }
INLINE void   QCXDR::translate(float *f) { }
INLINE void   QCXDR::translate(double *d) { }
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
INLINE void   QCXDR::translate(float*f) { byte_swap(f); }
INLINE void   QCXDR::translate(double*d) { byte_swap(d); }
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
  for (int k=0; k < n; k++) byte_swap(&f[k]);
}

INLINE void QCXDR::translate(double*d,int n)
{
  if (!d) return;
  for (int k=0; k < n; k++) byte_swap(&d[k]);
}

#endif /* BIGENDIAN */

#undef INLINE
