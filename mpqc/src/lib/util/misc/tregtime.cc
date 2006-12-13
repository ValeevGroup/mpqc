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

#include <unistd.h>

#include <util/misc/formio.h>
#include <util/misc/regtime.h>

#include <util/keyval/keyval.h>
static int (sc::KeyVal::*force_keyval_link)(const char*) = &sc::KeyVal::exists;

using namespace std;
using namespace sc;

int
main()
{
  int i;
  Ref<RegionTimer> tim = new RegionTimer("top", 1, 1);
  tim->enter("main");

  tim->enter("x");
  double x = 0.0;
  for (i=0; i<10000000; i++) {
      x += 0.0001;
    }
  tim->enter("subx");
  sleep(2);
  tim->exit("subx");
  ExEnv::outn() << indent << " x = " << x << endl;
  tim->exit("x");
  tim->enter("a");
  double a = 0.0;
  for (i=0; i<10000000; i++) {
      a += 0.0001;
    }
  tim->enter("subx");
  sleep(1);
  tim->exit("subx");
  ExEnv::outn() << indent << " a = " << a << endl;
  tim->exit("a");
  tim->enter("y");
  double y = 0.0;
  for (i=0; i<10000000; i++) {
      y += 0.0001;
    }
  ExEnv::outn() << indent << " y = " << y << endl;
  tim->change("z", "y");
  double z = 0.0;
  for (i=0; i<10000000; i++) {
      z += 0.0001;
    }
  ExEnv::outn() << " z = " << z << endl;
  tim->exit();
  tim->exit("main");

  Ref<RegionTimer> mtim = new RegionTimer("merged_top", 1, 1);
  mtim->enter("merged_a");
  for (int k=0; k<300000000; k++) { z += 0.0001; }
  mtim->change("merged_c");
  mtim->enter("merged_c_suba");
  for (int k=0; k<100000000; k++) { z += 0.0001; }
  mtim->change("merged_c_subb");
  mtim->exit("merged_c_subb");
  mtim->change("merged_b");
  for (int k=0; k<200000000; k++) { z += 0.0001; }
  mtim->exit("merged_b");

  tim->enter("merge_point");
  tim->merge(mtim);
  tim->exit("merge_point");

  tim->enter("merge_point*2");
  tim->merge(mtim);
  mtim->enter("merged_z");
  mtim->exit("merged_z");
  tim->merge(mtim);
  tim->exit("merge_point*2");

  RegionTimer::set_default_regiontimer(tim);

  Timer timertest("timertest");
  Timer r1("r1");
  r1.reset("r2");
  Timer r3("r3");
  r3.reset();
  r1.reset();

  {
    Timer x1("x1");
    Timer x2("x2");
    Timer x3("x3");
    // destructors are called in the reverse of the order of declaration
  }

  Timer y1("y1");
  y1.reset();

  mtim->print();

  tim->print();

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
