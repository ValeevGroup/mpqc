//
// psiotest.cc
//
// Copyright (C) 2004 Edward Valeev.
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

// a simple program to test the class stuff

#include <iostream>
#include <stdio.h>

#include <util/psi3/libpsio/psio.h>

using namespace std;
using namespace psi3::libpsio;

main()
{
  psio_init();
  psio_open(32,0);

  double A = .1111;  
  psio_write_entry(32, ":A", (char *)&A, sizeof(double));
  cout << "Wrote entry :A : value = " << A << endl;
  psio_close(32,1);

  psio_open(32,1);
  psio_read_entry(32, ":A", (char *)&A, sizeof(double));
  cout << "Read entry :A : value = " << A << endl;
  cout << "Table of contents of file 32:" << endl;
  psio_tocprint(32,stdout);
  psio_close(32,1);

  psio_done();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
