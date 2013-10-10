//
// state_file.h
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

#ifndef _util_state_state_file_h
#define _util_state_state_file_h

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>

namespace sc {

/**  @ingroup CoreState
 *   Writes state information to files.
 */
class StateOutFile: public StateOut {
  private:
    // do not allow copy constructor or assignment
    StateOutFile(const StateOutFile&);
    void operator=(const StateOutFile&);
  protected:
    int opened_;
    std::streambuf *buf_;
  public:
    /// State information will be written to ExEnv::outn().
    StateOutFile();
    /// State information will be written to s.
    StateOutFile(std::ostream& s);
    /// State information will be written to name.
    StateOutFile(const char *name);

    ~StateOutFile();

    /// State information will be written to name.
    virtual int open(const char *name);
    /// Flush the output stream.
    virtual void flush();
    /// Close the output stream.
    virtual void close();
  };

/**  @ingroup CoreState
 *   Reads state information from a file.
 */
class StateInFile: public StateIn {
  private:
    // do not allow copy constructor or assignment
    StateInFile(const StateInFile&);
    void operator=(const StateInFile&);
  protected:
    int opened_;
    std::streambuf *buf_;
  public:
    /// State information will be obtained from cin.
    StateInFile();
    /// State information will be obtained from fp.
    StateInFile(std::istream& s);
    /// State information will be obtained from name.
    StateInFile(const char *name);

    ~StateInFile();

    /// State information will be obtained from name.
    virtual int open(const char *name);
    /// Close the output file.
    virtual void close();
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
