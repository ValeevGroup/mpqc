//
// wfntest.cc
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
#include <util/state/stateio.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>

#include <chemistry/qc/wfn/obwfn.h>

using namespace std;

// Force linkages:
static ForceLink<HCoreWfn> fl0;

main(int argc, char *argv[])
{
  const char *input = (argc > 1) ? argv[1] : SRCDIR "/wfntest.kv";

  Ref<KeyVal> rpkv(new ParsedKeyVal(input));
  
  // the output stream is standard out
  ostream &o = cout;

  Ref<OneBodyWavefunction> wfn;
  wfn << rpkv->describedclassvalue("wavefunction");
  if (wfn.null()) {
    cerr << node0 << "wfn is null\n";
    exit(1);
  }

  wfn->overlap()->print("overlap");
  wfn->core_hamiltonian()->print("Hcore");
  wfn->hcore_guess()->print("guess vector");
  wfn->density()->print("density");
  wfn->natural_orbitals()->print("natural orbitals");
  wfn->natural_density()->print("natural density");

  //wfn->print(o);
  //o << endl;

  Ref<OneBodyWavefunction> oldwfn;
  oldwfn << rpkv->describedclassvalue("pwavefunction");
  
  RefSCMatrix evecs = wfn->projected_eigenvectors(oldwfn);

  evecs.print("projected wavefunction");

  StateOutBin so("wfn.ckpt");
  SavableState::save_state(wfn.pointer(),so);
  so.close();

  Ref<MolecularEnergy> me;
  StateInBin si("wfn.ckpt");
  me << SavableState::restore_state(si);
  
  me->print(o);
  o << me->value();
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
