//
// bugtest.cc
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

#include <iostream>
#include <util/keyval/keyval.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/state/state_text.h>

using namespace std;
using namespace sc;

double bugtest_global_var_;
double t;

void
y()
{
  t = 1.0/bugtest_global_var_;
}

void
x(const Ref<Debugger> &d)
{
  d->traceback();
  y();
}

int
main(int argc, char **argv)
{
  const char *infile = SRCDIR "/bugtest.in";

  Ref<KeyVal> keyval = new ParsedKeyVal(infile);

  Ref<Debugger> d;

  d << keyval->describedclassvalue("debug");
  if (d.null()) {
      d = new Debugger();
      d->handle_defaults();
      d->set_prefix(999);
      d->set_traceback_on_signal(1);
      d->set_debug_on_signal(1);
      d->set_exit_on_signal(0);
    }
  d->set_exec(argv[0]);

  d->traceback("no particular problem");
  d->debug("no particular problem");
  x(d);

  StateOutText o("state.dat");
  SavableState::save_state(d.pointer(),o);
  o.flush();

  StateInText i("state.dat");
  d << SavableState::restore_state(i);

  d = 0;
  cout << indent << "bugtest: done" << endl;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
