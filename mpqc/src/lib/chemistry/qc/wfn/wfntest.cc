//
// wfntest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <chemistry/qc/wfn/obwfn.h>

// Force linkages:
const ClassDesc &fl0 = HCoreWfn::class_desc_;

main(int argc, char *argv[])
{
  char *input = (argc > 1) ? argv[1] : SRCDIR "/wfntest.kv";

  RefKeyVal rpkv(new ParsedKeyVal(input));
  
  // the output stream is standard out
  ostream &o = cout;

  RefOneBodyWavefunction wfn = rpkv->describedclassvalue("wavefunction");
  if (wfn.null()) {
    cerr << node0 << "wfn is null\n";
    exit(1);
  }

  wfn->overlap()->print("overlap");
  wfn->core_hamiltonian()->print("Hcore");
  wfn->hcore_guess()->print("guess vector");

  //wfn->print(o);
  //o << endl;

  RefOneBodyWavefunction oldwfn = rpkv->describedclassvalue("pwavefunction");
  
  RefSCMatrix evecs = wfn->projected_eigenvectors(oldwfn);

  evecs.print("projected wavefunction");

  StateOutText so("wfn.ckpt");
  wfn.save_state(so);
  so.close();

  RefMolecularEnergy me;
  StateInText si("wfn.ckpt");
  me.restore_state(si);
  
  me->print(o);
  o << me->value();
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
