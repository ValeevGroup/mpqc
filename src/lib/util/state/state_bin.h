//
// state_bin.h
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

#ifndef _util_state_state_bin_h
#define _util_state_state_bin_h

#include <util/state/state_file.h>

namespace sc {

/**  @ingroup CoreState
 *   Save state to a binary file.
 */
class StateOutBin: public StateOutFile {
  private:
    int file_position_;
    // do not allow copy constructor or assignment
    StateOutBin(const StateOutBin&);
    void operator=(const StateOutBin&);
    /** This cannot be overridden, since it is called
        by this classes ctor (implicitly, through put_header()).
        This goes for some other members too. */
    int put_array_void(const void*,int);
  public:
    StateOutBin();
    StateOutBin(std::ostream&);
    StateOutBin(const char *);
    ~StateOutBin();

    int open(const char *name);
    void close();

    int use_directory();

    int tell();
    void seek(int loc);
    int seekable();
  };

/**  @ingroup CoreState
 *   Read objects written with StateOutBin.
 */
class StateInBin: public StateInFile {
  private:
    int file_position_;
    // do not allow copy constructor or assignment
    StateInBin(const StateInBin&);
    void operator=(const StateInBin&);
    /** These cannot be overridden, since they are called
        by this classes ctor (implicitly, through get_header()).
        This goes for other some members too. */
    int get_array_void(void*,int);
  public:
    StateInBin();
    StateInBin(const Ref<KeyVal> &);
    StateInBin(std::istream&);
    StateInBin(const char *);
    ~StateInBin();

    int open(const char *name);

    int use_directory();

    int tell();
    void seek(int loc);
    int seekable();
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
