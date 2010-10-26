
// translate.h -- data translation classes for StateIn and StateOut
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

#ifndef _util_state_translate_h
#define _util_state_translate_h

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#if defined(WORDS_BIGENDIAN)
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif

namespace sc {

/** Generic data translation.
 */
class TranslateData {
  public:
    TranslateData();
    virtual ~TranslateData();

    /// Returns a code for the type of format for the external data.
    virtual char format_code();

    /** A virtual constructor that choses a specialization based on
        the format code. */
    static TranslateData *vctor(char code);

    /** Translates to native format in-place.
        Similar routines exist for all the basic types. */
    virtual void to_native  (char *,   int n);
    /** Translates to external format in-place.
        Similar routines exist for all the basic types. */
    virtual void to_external(char *,   int n);
    virtual void to_native  (short *,  int n);
    virtual void to_external(short *,  int n);
    virtual void to_native  (unsigned int *, int n);
    virtual void to_external(unsigned int *, int n);
    virtual void to_native  (int *,    int n);
    virtual void to_external(int *,    int n);
    virtual void to_native  (unsigned long *,   int n);
    virtual void to_external(unsigned long *,   int n);
    virtual void to_native  (long *,   int n);
    virtual void to_external(long *,   int n);
    virtual void to_native  (float *,  int n);
    virtual void to_external(float *,  int n);
    virtual void to_native  (double *, int n);
    virtual void to_external(double *, int n);

    /** Translates to native format.
        Similar routines exist for all the basic types. */
    virtual void to_native  (char *target,   const void *source,   int n);
    /** Translates to external format.
        Similar routines exist for all the basic types. */
    virtual void to_external(void *target,   const char *source,   int n);
    virtual void to_native  (short *,  const void *,   int n);
    virtual void to_external(void *,   const short *,  int n);
    virtual void to_native  (unsigned int *,    const void *,   int n);
    virtual void to_external(void *,   const unsigned int *,    int n);
    virtual void to_native  (int *,    const void *,   int n);
    virtual void to_external(void *,   const int *,    int n);
    virtual void to_native  (unsigned long *,   const void *,   int n);
    virtual void to_external(void *,   const unsigned long *,   int n);
    virtual void to_native  (long *,   const void *,   int n);
    virtual void to_external(void *,   const long *,   int n);
    virtual void to_native  (float *,  const void *,   int n);
    virtual void to_external(void *,   const float *,  int n);
    virtual void to_native  (double *, const void *,   int n);
    virtual void to_external(void *,   const double *, int n);
};

/** Data translation to an external representation with bytes swapped.
 */
class TranslateDataByteSwap: public TranslateData {
  public:
    TranslateDataByteSwap();
    virtual ~TranslateDataByteSwap();

    /// Returns a code for the type of format for the external data.
    virtual char format_code();

    /// Overridden translation routines exist for all the basic types.
    virtual void to_native  (char *,   int n);
    /// Overridden translation routines exist for all the basic types.
    virtual void to_external(char *,   int n);
    virtual void to_native  (short *,  int n);
    virtual void to_external(short *,  int n);
    virtual void to_native  (unsigned int *, int n);
    virtual void to_external(unsigned int *, int n);
    virtual void to_native  (int *,    int n);
    virtual void to_external(int *,    int n);
    virtual void to_native  (unsigned long *,   int n);
    virtual void to_external(unsigned long *,   int n);
    virtual void to_native  (long *,   int n);
    virtual void to_external(long *,   int n);
    virtual void to_native  (float *,  int n);
    virtual void to_external(float *,  int n);
    virtual void to_native  (double *, int n);
    virtual void to_external(double *, int n);

    /// Overridden translation routines exist for all the basic types.
    virtual void to_native  (char *,   const void *,   int n);
    /// Overridden translation routines exist for all the basic types.
    virtual void to_external(void *,   const char *,   int n);
    virtual void to_native  (short *,  const void *,   int n);
    virtual void to_external(void *,   const short *,  int n);
    virtual void to_native  (unsigned int *,    const void *,   int n);
    virtual void to_external(void *,   const unsigned int *,    int n);
    virtual void to_native  (int *,    const void *,   int n);
    virtual void to_external(void *,   const int *,    int n);
    virtual void to_native  (unsigned long *,   const void *,   int n);
    virtual void to_external(void *,   const unsigned long *,   int n);
    virtual void to_native  (long *,   const void *,   int n);
    virtual void to_external(void *,   const long *,   int n);
    virtual void to_native  (float *,  const void *,   int n);
    virtual void to_external(void *,   const float *,  int n);
    virtual void to_native  (double *, const void *,   int n);
    virtual void to_external(void *,   const double *, int n);
};

#if BIGENDIAN
typedef TranslateDataByteSwap TranslateDataLittleEndian;
typedef TranslateData TranslateDataBigEndian;
#else
typedef TranslateDataByteSwap TranslateDataBigEndian;
typedef TranslateData TranslateDataLittleEndian;
#endif

class StateOut;

/** Convert data to other formats.
    The generated data is inserted into a StateOut object.
 */
class TranslateDataOut {
  private:
    StateOut *so_;
    TranslateData *translate_;
    // the translation buffer
    enum { bufsize = 8192 };
    char buf_[bufsize];
  protected:
    int putv(const void*d,int s);
  public:
    /** Write to s using the translation defined by t.
        The t argument will be deleted by this. */
    TranslateDataOut(StateOut*s, TranslateData*t);
    virtual ~TranslateDataOut();

    /** Translate and write the data. A similar member exists for
        each basic type. */
    virtual int put(const char*,int);
    virtual int put(const short*,int);
    virtual int put(const unsigned int*,int);
    virtual int put(const int*,int);
    virtual int put(const unsigned long*,int);
    virtual int put(const long*,int);
    virtual int put(const float*,int);
    virtual int put(const double*,int);

    /// Returns the translator.
    TranslateData *translator() { return translate_; }
};

class StateIn;

/** Convert data from other formats.
    The data is taken from a StateIn object.
 */
class TranslateDataIn {
  private:
    StateIn *si_;
    TranslateData *translate_;
  protected:
    int getv(void*d,int s);
  public:
    /** Input data will come from s.  The t argument will be deleted by this.
     */
    TranslateDataIn(StateIn*s, TranslateData *t);
    virtual ~TranslateDataIn();

    /** Read and translate data.  A similar member exists for each basic
        type. */
    virtual int get(char*,int);
    virtual int get(short*,int);
    virtual int get(unsigned int*,int);
    virtual int get(int*,int);
    virtual int get(unsigned long*,int);
    virtual int get(long*,int);
    virtual int get(float*,int);
    virtual int get(double*,int);

    /// Return the translator.
    TranslateData *translator() { return translate_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
