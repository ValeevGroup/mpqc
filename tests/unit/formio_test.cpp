//
// formiotest.cpp
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

#include "catch.hpp"

#include "mpqc/util/misc/formio.h"
#include <sstream>

using namespace std;
using namespace mpqc;

TEST_CASE("FormIO", "[formio]") {
  std::ostringstream oss;

  int elem = ios::xalloc();
  oss << "elem = " << elem << endl;

  oss << indent << "l0" << endl;
  oss << incindent;
  oss << indent << "l1" << endl;
  oss << incindent;
  oss << indent << "l2" << endl;
  oss << indent << "l2" << endl;
  long ind = FormIO::getindent(oss);
  oss << indent << "xyz = " << skipnextindent;
  FormIO::setindent(oss,FormIO::getindent(oss) + 6);
  oss << indent << "lxyz0" << endl;
  oss << indent << "lxyz1" << endl;
  oss << indent << "lxyz2" << endl;
  FormIO::setindent(oss,ind);
  oss << decindent;
  oss << indent << "l1" << endl;
  oss << decindent;
  oss << indent << "l0" << endl;

  oss << indent << mpqc::printf("%3d %10.5f",10,3.14) << endl;

  std::string ref_output("elem = 5\n\
l0\n\
  l1\n\
    l2\n\
    l2\n\
    xyz = lxyz0\n\
          lxyz1\n\
          lxyz2\n\
  l1\n\
l0\n\
 10    3.14000\n\
");

  REQUIRE(ref_output == oss.str());
}

