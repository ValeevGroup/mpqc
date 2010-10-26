//
// symmetrize.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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
#include <string.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>

using namespace std;
using namespace sc;

int
main(int argc, char *argv[])
{
  int i;

  if (argc < 2) {
    ExEnv::err0() << "usage: " << argv[0]
         << " input_file { keyword { tolerance } }" << endl;
    ExEnv::err0() << "  default keyword = \"molecule\"" << endl;
    ExEnv::err0() << "  default tolerance = \"1.0e-4\"" << endl;
    return 1;
  }

  char *infile = argv[1];
  Ref<KeyVal> kv(new ParsedKeyVal(infile));

  const char *keyword = argc>2?argv[2]:"molecule";
  Ref<Molecule> mol; mol << kv->describedclassvalue(keyword);

  const char *ctol = argc>3?argv[3]:"1.0e-4";
  double tol = atof(ctol);

  ExEnv::out0() << "Original molecule:" << endl;
  mol->print();
  
  Ref<PointGroup> highestpg = mol->highest_point_group(tol);
  ExEnv::out0() << "Point Group is " << highestpg->symbol() << endl;

  mol->set_point_group(highestpg, 10*tol);

  ExEnv::out0() << "Molecule at center of mass in highest point group:" << endl;
  mol->print();

  mol->cleanup_molecule();
  ExEnv::out0() << "cleaned molecule\n";
  mol->print();
  
  int nunique = mol->nunique();

  mol->transform_to_principal_axes();
  ExEnv::out0() << "cleaned molecule transformed to principle axes\n";
  mol->print();

  ExEnv::out0() << "resymmetrized molecule\n";
  mol->symmetrize();
  mol->print();

  mol->transform_to_symmetry_frame();
  ExEnv::out0() << "cleaned molecule transformed to symmetry frame\n";
  mol->print();

  ExEnv::out0() << scprintf("\nnunique=%d: ",nunique);
  for (i=0; i < nunique; i++) ExEnv::out0() << scprintf(" %d",mol->unique(i)+1);
  ExEnv::out0() << endl;
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
