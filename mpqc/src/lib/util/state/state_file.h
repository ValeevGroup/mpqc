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

#ifdef __GNUC__
#pragma interface
#endif

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

#include <util/state/state.h>

//. The \clsnmref{StateOutFile} provides a \clsnmref{StateOut}
//. which writes to files.  It is still abstract---one of its
//. derived classes, \clsnmref{StateOutFileText} or
//. \clsnmref{StateOutFileBin}, must be used to obtain a
//. \clsnmref{StateOut} object.  The
//. \clsnmref{StateOutFileText} class writes in a text format
//. and the \clsnmref{StateOutFileBin} writes in a binary
//. format.
class StateOutFile: public StateOut {
  private:
    // do not allow copy constructor or assignment
    StateOutFile(const StateOutFile&);
    void operator=(const StateOutFile&);
  protected:
    int opened_;
    streambuf *buf_;
  public:
    //. State information will be written to \srccd{stdout}.
    StateOutFile();
    //. State information will be written to \vrbl{fp}.
    StateOutFile(ostream& s);
    //. State information will be written to \filnm{name}.
    StateOutFile(const char *name);

    ~StateOutFile();

    //. State information will be written to \filnm{name}.
    virtual int open(const char *name);
    //. Miscellaneous file operations.
    virtual void flush();
    virtual void close();
  };

//. The \clsnm{StateInFile} provides a \clsnmref{StateIn} which
//. reads from files.  It is still abstract---one of its
//. derived classes, \clsnmref{StateInFileText} or
//. \clsnmref{StateInFileBin}, must be used to obtain a
//. \clsnmref{StateIn} object.  The \clsnmref{StateInFileText} class
//. reads with a text format and the \clsnmref{StateInFileBin}
//. reads with a binary format.
class StateInFile: public StateIn {
  private:
    // do not allow copy constructor or assignment
    StateInFile(const StateInFile&);
    void operator=(const StateInFile&);
  protected:
    int opened_;
    streambuf *buf_;
  public:
    //. State information will be obtained from \srccd{stdin}.
    StateInFile();
    //. State information will be obtained from \vrbl{fp}.
    StateInFile(istream& s);
    //. State information will be obtained from \filnm{name}.
    StateInFile(const char *name);

    ~StateInFile();

    //. State information will be obtained from \filnm{name}.
    virtual int open(const char *name);
    //. Miscellaneous file operations.
    virtual void close();
  };

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
