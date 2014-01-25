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

#include <chemistry/qc/wfn/density.h>
#include <chemistry/qc/wfn/obwfn.h>

using namespace std;
using namespace sc;

// Force linkages:
static ForceLink<HCoreWfn> fl0;
static ForceLink<SuperpositionOfAtomicDensities> fl1;

#include <chemistry/molecule/linkage.h>

void
density_test(const Ref<Wavefunction> &wfn, double resolution)
{
  wfn->ao_density()->print("AO Density Matrix");
  wfn->overlap()->print("SO Overlap Matrix");
  wfn->natural_density()->print("Natural Density");
  wfn->natural_orbitals()->print("Natural Orbitals");

  Ref<ElectronDensity> ed = new ElectronDensity(wfn);
  Ref<BatchElectronDensity> bed = new BatchElectronDensity(wfn);

  SCVector3 r(0,0,0);
  for (double x = 0.0; x < 2.1; x += 0.25) {
    r[0] = x;
    ed->set_x(x);
    bed->set_x(x);
    SCVector3 edg;
    ed->get_gradient(edg);
    ExEnv::out0() << scprintf(" ED: %12.8f; ", ed->value())
                  << edg
                  << std::endl;
    SCVector3 bedg;
    bed->get_gradient(bedg);
    ExEnv::out0() << scprintf("BED: %12.8f; ", bed->value())
                  << bedg
                  << std::endl;
  }

  SCVector3 upper, lower;
  bed->boundingbox(DBL_EPSILON,DBL_MAX,lower,upper);
  ExEnv::out0() << "BED BB: " << lower << ", " << upper << std::endl;
  ed->boundingbox(DBL_EPSILON,DBL_MAX,lower,upper);
  ExEnv::out0() << " ED BB: " << lower << ", " << upper << std::endl;

  SCVector3 boxsize = upper - lower;
  ExEnv::out0() << "boxsize = " << boxsize << std::endl;
  int nx = int(boxsize[0]/resolution);
  int ny = int(boxsize[1]/resolution);
  int nz = int(boxsize[2]/resolution);
  ExEnv::out0() << "evaluating " << nx*ny*nz << " points" << std::endl;

  double nele_bed = 0.0;
  for (int i=0; i<nx; i++) {
    SCVector3 r(i*resolution + lower[0],0.,0.);
    for (int j=0; j<ny; j++) {
      r[1] = j*resolution + lower[1];
      for (int k=0; k<nz; k++) {
        r[2] = k*resolution + lower[2];
        bed->set_x(r);
        nele_bed += bed->value();
      }
    }
  }
  nele_bed *= resolution*resolution*resolution;
  ExEnv::out0() << scprintf("BED Nele = %12.8f",nele_bed) << std::endl;

  double nele_ed = 0.0;
  for (int i=0; i<nx; i++) {
    SCVector3 r(i*resolution + lower[0],0.,0.);
    for (int j=0; j<ny; j++) {
      r[1] = j*resolution + lower[1];
      for (int k=0; k<nz; k++) {
        r[2] = k*resolution + lower[2];
        ed->set_x(r);
        nele_ed += ed->value();
      }
    }
  }
  nele_ed *= resolution*resolution*resolution;
  ExEnv::out0() << scprintf(" ED Nele = %12.8f",nele_ed) << std::endl;
}

int
main(int argc, char *argv[])
{
  const char *input = (argc > 1) ? argv[1] : SRCDIR "/wfntest.kv";

  Ref<KeyVal> rpkv(new ParsedKeyVal(input));
  
  // the output stream is standard out
  ostream &o = cout;

  Ref<OneBodyWavefunction> wfn;
  wfn << rpkv->describedclassvalue("wavefunction");
  if (wfn == 0) {
    ExEnv::err0() << "wfn is null\n";
    exit(1);
  }

  double resolution = rpkv->doublevalue("resolution");
  density_test(wfn,resolution);

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

  wfn->ao_density().print("projected AO density");

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
