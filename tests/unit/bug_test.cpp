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
#include <sstream>

#include "catch.hpp"

#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/bug.h"
#include "mpqc/util/core/formio.h"

using namespace std;
using namespace mpqc;

double bugtest_global_var_;
double t;

void
y()
{
  t = 1.0/bugtest_global_var_;
}

void
x(const std::shared_ptr<Debugger> &d)
{
  d->traceback();
  y();
}

struct cout_redirect {
    cout_redirect( std::streambuf * new_buffer )
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

private:
    std::streambuf * old;
};

TEST_CASE("Debugger", "[debugger]") {
  const char *infile = SRCDIR "/bug_test.json";
  ifstream input(infile);
  REQUIRE(input.good());

  KeyVal keyval;
  REQUIRE_NOTHROW(keyval.read_json(input));

  // redirect cout to a stringstream until the end of this test
  std::ostringstream ss;
  cout_redirect _x(ss.rdbuf());

  {
    auto d = std::make_shared<Debugger>();
    REQUIRE(d != nullptr);
    REQUIRE_NOTHROW(d->handle_defaults());
    REQUIRE_NOTHROW(d->set_prefix(999));
    REQUIRE_NOTHROW(d->set_traceback_on_signal(1));
    REQUIRE_NOTHROW(d->set_debug_on_signal(1));
    REQUIRE_NOTHROW(d->set_exit_on_signal(0));
  }

  std::shared_ptr<Debugger> d = keyval.class_ptr<Debugger>("debug");
  REQUIRE_NOTHROW(d->set_exec("unit_tests"));

  d->traceback("no particular problem");
  d->debug("no particular problem");
  REQUIRE_NOTHROW(x(d));

  REQUIRE_NOTHROW(d = 0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
