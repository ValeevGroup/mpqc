//
// gnutest.cc
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

#include <util/container/set.h>

int
main()
{
  int i;
  AVLSet<int> iset;

  for (i=10; i<20; i++) {
      court << "adding " << i << endl;
      iset.add(i);
    }

  for (i=1; i<15; i++) {
      court << "adding " << i << endl;
      iset.add(i);
    }

  int ten=10, fifteen=15;
  iset.del(ten);
  iset.del(fifteen);

  cout << "iset:" << endl;
  Pix I;
  for (I=iset.first(); I; iset.next(I)) {
      cout << " " << iset(I);
    }
  cout << endl;

  ///////////////////////////////////////////////////////////////

  Arrayset<int> aset;

  for (I=aset.first(); I; aset.next(I)) {
      cout << " " << aset(I);
    }
  cout << endl;

  for (i=10; i<20; i++) {
      cout << "adding " << i << endl;
      aset.add(i);
    }

  for (i=1; i<15; i++) {
      cout << "adding " << i << endl;
      aset.add(i);
    }

  aset.del(ten);
  aset.del(fifteen);

  cout << "aset:";
  for (I=aset.first(); I; aset.next(I)) {
      cout << " " << aset(I);
    }
  cout << endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
