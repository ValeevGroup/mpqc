//
// state_text.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _util_state_state_text_h
#define _util_state_state_text_h

#include <util/state/state_file.h>

namespace sc {

/**  @ingroup CoreState
 *   Writes out state information in an almost human readable format.

 StateOutText is intended for debugging only.  The state information can
 read in again with StateInText.
 */
class StateOutText: public StateOutFile {
  private:
    // do not allow copy constructor or assignment
    StateOutText(const StateOutText&);
    void operator=(const StateOutText&);
  protected:
    int no_newline_;
    int no_array_;
    void no_newline();
    void no_array();
    void newline();
    void start_array();
    void end_array();
    int putobject(const Ref<SavableState> &);
    int putparents(const ClassDesc*);
  public:
    StateOutText();
    StateOutText(std::ostream& s);
    StateOutText(const char *);
    ~StateOutText();
    int putstring(const char*);
    int put_array_char(const char*,int);
    int put_array_uint(const unsigned int*,int);
    int put_array_int(const int*,int);
    int put_array_ulong(const unsigned long*,int);
    int put_array_long(const long*,int);
    int put_array_float(const float*,int);
    int put_array_double(const double*,int);
    int put(const ClassDesc*);
    int put(const std::string &);
    int put(char r);
    int put(unsigned int r);
    int put(int r);
    int put(unsigned long r);
    int put(long r);
    int put(bool r);
    int put(float r);
    int put(double r);
    int put(const char*,int);
    int put(const unsigned int*,int);
    int put(const int*,int);
    int put(const unsigned long*,int);
    int put(const long*,int);
    int put(const float*,int);
    int put(const double*,int);
  };

/**  @ingroup CoreState
 *   Reads state information written with StateOutText.
 */
class StateInText: public StateInFile {
  private:
    // do not allow copy constructor or assignment
    StateInText(const StateInText&);
    void operator=(const StateInText&);
  protected:
    int newlines_;
    int no_newline_;
    int no_array_;
    void no_newline();
    void no_array();

    int read(char*);
    int read(unsigned int&);
    int read(int&);
    int read(unsigned long&);
    int read(long&);
    int read(bool&);
    int read(float&);
    int read(double&);
    void newline();
    void start_array();
    void end_array();
    int  getobject(Ref<SavableState> &);

    void abort();
  public:
    StateInText();
    StateInText(std::istream& s);
    StateInText(const char *);
    StateInText(const Ref<KeyVal> &);
    ~StateInText();
    int getstring(char*&);
    int get_array_char(char*,int);
    int get_array_uint(unsigned int*,int);
    int get_array_int(int*,int);
    int get_array_ulong(unsigned long*,int);
    int get_array_long(long*,int);
    int get_array_float(float*,int);
    int get_array_double(double*,int);
    int get(const ClassDesc**);
    int get(std::string&);
    int get(char&r, const char *key = 0);
    int get(unsigned int&r, const char *key = 0);
    int get(int&r, const char *key = 0);
    int get(unsigned long&r, const char *key = 0);
    int get(long&r, const char *key = 0);
    int get(bool&r, const char *key = 0);
    int get(float&r, const char *key = 0);
    int get(double&r, const char *key = 0);
    int get(char*&);
    int get(unsigned int*&);
    int get(int*&);
    int get(unsigned long*&);
    int get(long*&);
    int get(float*&);
    int get(double*&);
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
