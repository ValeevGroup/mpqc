//
// reftest.cc
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

#include <util/ref/ref.h>
#include <util/ref/reftestx.h>

//template class Ref<X>;

int
main()
{
  {
  Ref<X> x1 = new X;
  RefX x2 = new X;
#if REF_MANAGE
  x2->unmanage();
#endif
#if REF_CHECK_STACK
    {
      X z;
      RefX zz(&z);
    }
#endif

  x2 = x1.pointer();
  if (x1 != x1) abort();
  if (x2 != x2) abort();

  int i;
  for (i=1000000; i; i--) {
      x1->reference();
    }
  for (i=1000000; i; i--) {
      x1->dereference();
    }
  for (i=1000000; i; i--) {
      Ref<X> y = x1;
    }

  cerr << "nx = " << X::nx << " (inner scope)" << endl;
  }

   cerr << "nx = " << X::nx << " (outer scope)" << endl;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
