//
// formiotest.cc
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

#include <util/misc/formio.h>

using namespace std;
using namespace sc;

int
main(int argc, char* argv[])
{
  int elem = ios::xalloc();
  cout << " elem = " << elem << endl;

  cout << indent << "l0" << endl;
  cout << incindent;
  cout << indent << "l1" << endl;
  cout << incindent;
  cout << indent << "l2" << endl;
  cout << indent << "l2" << endl;
  long ind = SCFormIO::getindent(cout);
  cout << indent << "xyz = " << skipnextindent;
  SCFormIO::setindent(cout,SCFormIO::getindent(cout) + 6);
  cout << indent << "lxyz0" << endl;
  cout << indent << "lxyz1" << endl;
  cout << indent << "lxyz2" << endl;
  SCFormIO::setindent(cout,ind);
  cout << decindent;
  cout << indent << "l1" << endl;
  cout << decindent;
  cout << indent << "l0" << endl;

  cout << indent << scprintf("%3d %10.5f",10,3.14) << endl;
}

