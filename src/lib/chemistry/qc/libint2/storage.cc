//
// storage.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/libint2/int2e.h>
#include <chemistry/qc/libint2/storage.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

void
Int2eLibint2::init_storage(size_t size)
{
  if (size) storage_ = size;
  else storage_ = 0;
}

void
Int2eLibint2::done_storage()
{
  storage_ = 0;
}


void
Int2eLibint2::check_storage_() const
{
  //  if (storage_used_ > storage_) {
  //    ExEnv::err() << scprintf("Not enough storage allowed for integrals evaluators:") << endl;
  //    ExEnv::err() << scprintf("storage_ = %ld",storage_) << endl;
  //    ExEnv::err() << scprintf("storage_used_ = %ld",storage_used_) << endl;
  //    abort();
  //  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
