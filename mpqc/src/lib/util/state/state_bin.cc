//
// state_bin.cc
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

#include <util/class/class.h>
#include <util/state/state.h>

StateOutBin::StateOutBin() :
  StateOutFile()
{
}

StateOutBin::StateOutBin(ostream& s) :
  StateOutFile(s)
{
}

StateOutBin::StateOutBin(const char *path) :
  StateOutFile(path)
{
}

StateOutBin::~StateOutBin()
{
}

StateInBin::StateInBin() :
  StateInFile()
{
}

StateInBin::StateInBin(istream& s) :
  StateInFile(s)
{
}

StateInBin::StateInBin(const char *path) :
  StateInFile(path)
{
}

StateInBin::~StateInBin()
{
}

////////////////////////////////////////////////////////////////

int StateOutBin::put_array_void(const void*p,int size)
{
  if (buf_->sputn((const char *)p,size) != size) {
      cerr << "StateOutBin::put_array_void: failed" << endl;
      abort();
    }
  return size;
}

int StateInBin::get_array_void(void*p,int size)
{
  if (buf_->sgetn((char*)p,size) != size) {
      cerr << "StateInBin::get_array_void: failed" << endl;
      abort();
    }
  return size;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
