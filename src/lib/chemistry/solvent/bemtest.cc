//
// bemtest.cc
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

#include <stdio.h>
#include <util/misc/ieee.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molshape.h>
#include <chemistry/solvent/bem.h>

using namespace sc;

int
main(int argc, char* argv[])
{
  // Abort on floating point errors.
  ieee_trap_errors();

  Ref<KeyVal> keyval = new ParsedKeyVal(SRCDIR "/bemtest.in");

  Ref<BEMSolvent> solvent;
  solvent << keyval->describedclassvalue("solvent");

  solvent->init();
  solvent->init_system_matrix();
  solvent->done();

  solvent->init();
  solvent->init_system_matrix();
  solvent->done();

  ConnollyShape::print_counts();
  CS2Sphere::print_counts();

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
