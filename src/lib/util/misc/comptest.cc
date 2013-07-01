//
// comptest.cc
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

#include <iomanip>
#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif
#include <util/misc/compute.h>

using namespace std;
using namespace sc;

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class NCResult<int>;
#endif

class A: public Compute {
  public:
    A();
    Resultint i;
    AccResultdouble a;
    void compute();
    void print();
};

A::A():i(this),a(this)
{
  a.set_desired_accuracy(0.1);
  a.set_actual_accuracy(0.1);
  a.result_noupdate() = 0.0;
  a.computed() = 1;
}

void
A::compute()
{
  cout << "computing";
  if (i.needed()) {
      i.result_noupdate() = 5;
      i.computed() = 1;
      cout << " i";
    }
  if (a.needed()) {
      a.result_noupdate() += 0.001;
      a.computed() = 1;
      cout << " a";
      a.set_actual_accuracy(a.desired_accuracy());
    }
  cout << endl;
}

void
A::print()
{
  cout << "A: i = " << (int) i << ", a = "
       << setw(5) << setprecision(3) << (double)a << endl;
}

int
main(int argc, char* argv[])
{
  A a;

  cout << "should not compute a" << endl;
  a.print();

  a.a.set_desired_accuracy(0.01);

  cout << "should compute a" << endl;
  a.print();

  a.a.set_desired_accuracy(0.1);

  cout << "should not compute a" << endl;
  a.print();

  a.a.set_desired_accuracy(0.01);

  cout << "should not compute a" << endl;
  a.print();

  a.a.set_desired_accuracy(0.001);

  cout << "should compute a" << endl;
  a.print();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
