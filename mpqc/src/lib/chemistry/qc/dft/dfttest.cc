//
// dfttest.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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
#include <util/group/pregtime.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/dft/integrator.h>

#include <chemistry/qc/dft/linkage.h>
#include <chemistry/qc/scf/linkage.h>

int
main(int argc, char**argv)
{
  char *input =      (argc > 1)? argv[1] : SRCDIR "/dfttest.in";

  // open keyval input
  RefKeyVal keyval(new ParsedKeyVal(input));

  RefMessageGrp grp;
  grp = keyval->describedclassvalue("message");
  if (grp.nonnull()) MessageGrp::set_default_messagegrp(grp);
  else grp = MessageGrp::get_default_messagegrp();

  RefRegionTimer tim;
  tim = new ParallelRegionTimer(grp,"dfttest",1,0);
  RegionTimer::set_default_regiontimer(tim);
  tim->enter("input");
  
  if (keyval->exists("matrixkit"))
    SCMatrixKit::set_default_matrixkit(keyval
                                       ->describedclassvalue("matrixkit"));

  RefDenFunctional functional = keyval->describedclassvalue("functional");
  RefDenIntegrator integrator = keyval->describedclassvalue("integrator");
  RefWavefunction  wfn        = keyval->describedclassvalue("wfn");
  RefWavefunction  hfacm      = keyval->describedclassvalue("hfacm");

  integrator->set_wavefunction(wfn);

  tim->exit("input");

  tim->enter("energy");
  cout << "energy = " << wfn->energy() << endl;
  tim->exit("energy");

  tim->enter("integration");
  integrator->integrate(functional);
  tim->exit("integration");

  tim->enter("hfacm");
  cout << scprintf("hfacm e = % 14.10f", hfacm->energy()) << endl;
  tim->exit("hfacm");

  //tim->print(cout);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
