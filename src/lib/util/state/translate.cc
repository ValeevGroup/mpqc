
// translate.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/state/state.h>
#include <util/state/translate.h>
#include <util/state/stateio.h>

using namespace sc;

////////////////////////////////////////////////////////////////////////

static inline void
swap(char *d, int c1, int c2)
{
  char tmp = d[c1];
  d[c1] = d[c2];
  d[c2] = tmp;
}

static inline void
swap(char *d, const char *e, int c1, int c2)
{
  d[c1] = e[c2];
  d[c2] = e[c1];
}

static inline void
byte_swap2(char*d)
{
  swap(d,0,1);
}

static inline void
byte_swap2(char*d, const char *e)
{
  swap(d,e,0,1);
}

static inline void
byte_swap4(char*d)
{
  swap(d,0,3);
  swap(d,1,2);
}

static inline void
byte_swap4(char*d, const char *e)
{
  swap(d,e,0,3);
  swap(d,e,1,2);
}

static inline void
byte_swap8(char*d)
{
  swap(d,0,7);
  swap(d,1,6);
  swap(d,2,5);
  swap(d,3,4);
}

static inline void
byte_swap8(char*d, const char *e)
{
  swap(d,e,0,7);
  swap(d,e,1,6);
  swap(d,e,2,5);
  swap(d,e,3,4);
}

static inline void
byte_swap16(char*d)
{
  swap(d,0,15);
  swap(d,1,14);
  swap(d,2,13);
  swap(d,3,12);
  swap(d,4,11);
  swap(d,5,10);
  swap(d,6, 9);
  swap(d,7, 8);
}

static inline void
byte_swap16(char*d, const char *e)
{
  swap(d,e,0,15);
  swap(d,e,1,14);
  swap(d,e,2,13);
  swap(d,e,3,12);
  swap(d,e,4,11);
  swap(d,e,5,10);
  swap(d,e,6, 9);
  swap(d,e,7, 8);
}

static inline void
byte_swap2(void*data,int n)
{
  char *d = (char*)data;
  for (int i=0; i<n; i++) byte_swap2(&d[i<<1]);
}

static inline void
byte_swap2(void*data,const void*edata,int n)
{
  char *d = (char*)data;
  char *e = (char*)edata;
  for (int i=0; i<n; i++) byte_swap2(&d[i<<1],&e[i<<1]);
}

static inline void
byte_swap4(void*data,int n)
{
  char *d = (char*)data;
  for (int i=0; i<n; i++) byte_swap4(&d[i<<2]);
}

static inline void
byte_swap4(void*data,const void*edata,int n)
{
  char *d = (char*)data;
  char *e = (char*)edata;
  for (int i=0; i<n; i++) byte_swap4(&d[i<<2],&e[i<<2]);
}

static inline void
byte_swap8(void*data,int n)
{
  char *d = (char*)data;
  for (int i=0; i<n; i++) byte_swap8(&d[i<<3]);
}

static inline void
byte_swap8(void*data,const void*edata,int n)
{
  char *d = (char*)data;
  char *e = (char*)edata;
  for (int i=0; i<n; i++) byte_swap8(&d[i<<3],&e[i<<3]);
}

static inline void
byte_swap16(void*data,int n)
{
  char *d = (char*)data;
  for (int i=0; i<n; i++) byte_swap16(&d[i<<4]);
}

static inline void
byte_swap16(void*data,const void*edata,int n)
{
  char *d = (char*)data;
  char *e = (char*)edata;
  for (int i=0; i<n; i++) byte_swap16(&d[i<<4],&e[i<<4]);
}

#define BST(T) \
static inline void byte_swap(T*d,int n) \
{ \
  if (sizeof(T) == 2) byte_swap2(d,n); \
  else if (sizeof(T) == 4) byte_swap4(d,n); \
  else if (sizeof(T) == 8) byte_swap8(d,n); \
  else if (sizeof(T) == 16) byte_swap16(d,n); \
}

#define BSTV(T) \
static inline void byte_swap(T*d,const void*e,int n) \
{ \
  if (sizeof(T) == 2) byte_swap2(d,e,n); \
  else if (sizeof(T) == 4) byte_swap4(d,e,n); \
  else if (sizeof(T) == 8) byte_swap8(d,e,n); \
  else if (sizeof(T) == 16) byte_swap16(d,e,n); \
}

#define BSVT(T) \
static inline void byte_swap(void*d,const T*e,int n) \
{ \
  if (sizeof(T) == 2) byte_swap2(d,e,n); \
  else if (sizeof(T) == 4) byte_swap4(d,e,n); \
  else if (sizeof(T) == 8) byte_swap8(d,e,n); \
  else if (sizeof(T) == 16) byte_swap16(d,e,n); \
}

BSVT(short);
BSVT(unsigned int);
BSVT(int);
BSVT(unsigned long);
BSVT(long);
BSVT(float);
BSVT(double);

BSTV(short);
BSTV(unsigned int);
BSTV(int);
BSTV(unsigned long);
BSTV(long);
BSTV(float);
BSTV(double);

BST(short);
BST(unsigned int);
BST(int);
BST(unsigned long);
BST(long);
BST(float);
BST(double);

////////////////////////////////////////////////////////////////////////

TranslateData::TranslateData() {}
TranslateData::~TranslateData() {}

char
TranslateData::format_code()
{
#if BIGENDIAN
  return 'b';
#else
  return 'l';
#endif
}

TranslateData *
TranslateData::vctor(char code)
{
  if (code == 'b') return new TranslateDataBigEndian;
  if (code == 'l') return new TranslateDataLittleEndian;
  return 0;
}

void TranslateData::to_native  (char *,   int n) {}
void TranslateData::to_external(char *,   int n) {}
void TranslateData::to_native  (short *,  int n) {}
void TranslateData::to_external(short *,  int n) {}
void TranslateData::to_native  (unsigned int *, int n) {}
void TranslateData::to_external(unsigned int *, int n) {}
void TranslateData::to_native  (int *,    int n) {}
void TranslateData::to_external(int *,    int n) {}
void TranslateData::to_native  (unsigned long *,   int n) {}
void TranslateData::to_external(unsigned long *,   int n) {}
void TranslateData::to_native  (long *,   int n) {}
void TranslateData::to_external(long *,   int n) {}
void TranslateData::to_native  (float *,  int n) {}
void TranslateData::to_external(float *,  int n) {}
void TranslateData::to_native  (double *, int n) {}
void TranslateData::to_external(double *, int n) {}
void TranslateData::to_native  (char *d,   const void *s,   int n)
{ memcpy(d,s,n); }
void TranslateData::to_external(void *d,   const char *s,   int n)
{ memcpy(d,s,n); }
void TranslateData::to_native  (short *d,  const void *s,   int n)
{ memcpy(d,s,n*sizeof(short)); }
void TranslateData::to_external(void *d,   const short *s,  int n)
{ memcpy(d,s,n*sizeof(short)); }
void TranslateData::to_native  (unsigned int *d, const void *s, int n)
{ memcpy(d,s,n*sizeof(unsigned int)); }
void TranslateData::to_external(void *d, const unsigned int *s, int n)
{ memcpy(d,s,n*sizeof(unsigned int)); }
void TranslateData::to_native  (int *d,    const void *s,   int n)
{ memcpy(d,s,n*sizeof(int)); }
void TranslateData::to_external(void *d,   const int *s,    int n)
{ memcpy(d,s,n*sizeof(int)); }
void TranslateData::to_native  (unsigned long *d,   const void *s,   int n)
{ memcpy(d,s,n*sizeof(unsigned long)); }
void TranslateData::to_external(void *d,   const unsigned long *s,   int n)
{ memcpy(d,s,n*sizeof(unsigned long)); }
void TranslateData::to_native  (long *d,   const void *s,   int n)
{ memcpy(d,s,n*sizeof(long)); }
void TranslateData::to_external(void *d,   const long *s,   int n)
{ memcpy(d,s,n*sizeof(long)); }
void TranslateData::to_native  (float *d,  const void *s,   int n)
{ memcpy(d,s,n*sizeof(float)); }
void TranslateData::to_external(void *d,   const float *s,  int n)
{ memcpy(d,s,n*sizeof(float)); }
void TranslateData::to_native  (double *d, const void *s,   int n)
{ memcpy(d,s,n*sizeof(double)); }
void TranslateData::to_external(void *d,   const double *s, int n)
{ memcpy(d,s,n*sizeof(double)); }

////////////////////////////////////////////////////////////////////////

TranslateDataByteSwap::TranslateDataByteSwap() {}
TranslateDataByteSwap::~TranslateDataByteSwap() {}

char
TranslateDataByteSwap::format_code()
{
#if BIGENDIAN
  return 'l';
#else
  return 'b';
#endif
}

void TranslateDataByteSwap::to_native  (char *d,   int n) {}
void TranslateDataByteSwap::to_external(char *d,   int n) {}
void TranslateDataByteSwap::to_native  (short *d,  int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_external(short *d,  int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (unsigned int *d,    int n)
{ byte_swap(d,n); }
void TranslateDataByteSwap::to_external(unsigned int *d,    int n)
{ byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (int *d,    int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_external(int *d,    int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (unsigned long *d,   int n)
{ byte_swap(d,n); }
void TranslateDataByteSwap::to_external(unsigned long *d,   int n)
{ byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (long *d,   int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_external(long *d,   int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (float *d,  int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_external(float *d,  int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (double *d, int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_external(double *d, int n) { byte_swap(d,n); }
void TranslateDataByteSwap::to_native  (char *d,   const void *s,   int n)
{ memcpy(d,s,n); }
void TranslateDataByteSwap::to_external(void *d,   const char *s,   int n)
{ memcpy(d,s,n); }
void TranslateDataByteSwap::to_native  (short *d,  const void *s,   int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d,   const short *s,  int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_native  (unsigned int *d, const void *s, int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d, const unsigned int *s, int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_native  (int *d,    const void *s,   int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d,   const int *s,    int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_native  (unsigned long *d, const void *s, int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d, const unsigned long *s, int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_native  (long *d,   const void *s,   int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d,   const long *s,   int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_native  (float *d,  const void *s,   int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d,   const float *s,  int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_native  (double *d, const void *s,   int n)
{ byte_swap(d,s,n); }
void TranslateDataByteSwap::to_external(void *d,   const double *s, int n)
{ byte_swap(d,s,n); }

////////////////////////////////////////////////////////////////////////

TranslateDataOut::TranslateDataOut(StateOut*so, TranslateData*t):
  so_(so),
  translate_(t)
{
}

TranslateDataOut::~TranslateDataOut()
{
  delete translate_;
}

inline int
TranslateDataOut::putv(const void*d,int s)
{
  return so_->put_array_void(d,s);
}

int
TranslateDataOut::put(const char*d,int s)
{
  const int bsize = bufsize;
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l);
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const short*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const unsigned int*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const int*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const unsigned long*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const long*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const float*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

int
TranslateDataOut::put(const double*d,int s)
{
  const int bsize = bufsize/sizeof(*d);
  int o=0,r=0;
  while (s) {
      int l = (s>bsize?bsize:s);
      translate_->to_external(buf_,&d[o],l);
      r += putv(buf_,l*sizeof(*d));
      s-=l;
      o+=l;
    }
  return r;
}

////////////////////////////////////////////////////////////////////////

TranslateDataIn::TranslateDataIn(StateIn*si,TranslateData *t):
  si_(si),
  translate_(t)
{
}

TranslateDataIn::~TranslateDataIn()
{
  delete translate_;
}

inline int
TranslateDataIn::getv(void*d,int s)
{
  return si_->get_array_void(d,s);
}

int
TranslateDataIn::get(char*d,int s)
{
  int r = getv(d,s);
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(short*d,int s)
{
  int r = getv(d,s*sizeof(short));
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(unsigned int*d,int s)
{
  int r = getv(d,s*sizeof(unsigned int));
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(int*d,int s)
{
  int r = getv(d,s*sizeof(int));
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(unsigned long*d,int s)
{
  int r = getv(d,s*sizeof(unsigned long));
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(long*d,int s)
{
  int r = getv(d,s*sizeof(long));
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(float*d,int s)
{
  int r = getv(d,s*sizeof(float));
  translate_->to_native(d,s);
  return r;
}

int
TranslateDataIn::get(double*d,int s)
{
  int r = getv(d,s*sizeof(double));
  translate_->to_native(d,s);
  return r;
}

////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
