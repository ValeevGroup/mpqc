//
// tregtime.cc
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

#include <util/misc/regtime.h>

int
main()
{
  int i;
  RefRegionTimer tim = new RegionTimer("top", 1, 1);
  tim->enter("main");

  tim->enter("x");
  double x = 0.0;
  for (i=0; i<10000000; i++) {
      x += 0.0001;
    }
  tim->enter("subx");
  sleep(2);
  tim->exit("subx");
  cout << " x = " << x << endl;
  tim->exit("x");
  tim->enter("a");
  double a = 0.0;
  for (i=0; i<10000000; i++) {
      a += 0.0001;
    }
  tim->enter("subx");
  sleep(1);
  tim->exit("subx");
  cout << " a = " << a << endl;
  tim->exit("a");
  tim->enter("y");
  double y = 0.0;
  for (i=0; i<10000000; i++) {
      y += 0.0001;
    }
  cout << " y = " << y << endl;
  tim->change("z", "y");
  double z = 0.0;
  for (i=0; i<10000000; i++) {
      z += 0.0001;
    }
  cout << " z = " << z << endl;
  tim->exit();
  tim->exit("main");

  tim->print();

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
