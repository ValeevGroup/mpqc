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

StateOutBin::StateOutBin(FILE* fp) :
  StateOutFile(fp)
{
}

StateOutBin::StateOutBin(const char *path, const char * mode) :
  StateOutFile(path,mode)
{
}

StateOutBin::~StateOutBin()
{
}

StateInBin::StateInBin() :
  StateInFile()
{
}

StateInBin::StateInBin(FILE* fp) :
  StateInFile(fp)
{
}

StateInBin::StateInBin(const char *path, const char * mode) :
  StateInFile(path,mode)
{
}

StateInBin::~StateInBin()
{
}

////////////////////////////////////////////////////////////////

int StateOutBin::put_array_void(const void*p,int size)
{
  return fwrite(p,1,size,fp_);
}

int StateInBin::get_array_void(void*p,int size)
{
  return fread(p,1,size,fp_);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
