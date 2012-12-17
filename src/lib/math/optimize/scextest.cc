//
// scextest.cc
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

#include <util/keyval/keyval.h>
#include <util/state/stateio.h>
#include <util/state/state_text.h>
#include <math/scmat/local.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrap.h>
#include <math/optimize/scextrapmat.h>

using namespace std;
using namespace sc;

// Force linkages:
#ifndef __PIC__
static ForceLink<DIIS> fl0;
#endif

int
main(int argc, char* argv[])
{
  int i;
  
  Ref<KeyVal> keyval = new ParsedKeyVal( SRCDIR "/scextest.in");

  Ref<SelfConsistentExtrapolation> extrap;
  extrap << keyval->describedclassvalue("scextrap");

  RefSCDimension dim = new SCDimension(3, "test_dim");
  Ref<SCMatrixKit> kit = new LocalSCMatrixKit;

  RefSymmSCMatrix datamat(dim,kit);
  datamat.assign(0.0);
  datamat->shift_diagonal(2.0);

  RefDiagSCMatrix val(dim,kit);
  RefSCMatrix vec(dim,dim,kit);

  // solve f(x) = x

  i = 0;
  while (i < 100 && !extrap->converged()) {
      datamat.diagonalize(val,vec);
      for (int j=0; j<datamat.dim().n(); j++) {
          double v = val.get_element(j);
          val.set_element(j, sqrt(v));
        }
      RefSymmSCMatrix newdatamat(dim,kit);
      newdatamat.assign(0.0);
      newdatamat.accumulate_transform(vec, val);
      RefSymmSCMatrix errormat = newdatamat - datamat;

      datamat.assign(newdatamat);
      Ref<SCExtrapData> data = new SymmSCMatrixSCExtrapData(datamat);
      Ref<SCExtrapError> error = new SymmSCMatrixSCExtrapError(errormat);

      ExEnv::out0() << "Iteration " << i << ":" << endl;

      datamat.print("Datamat:");
      errormat.print("Errormat:");

      extrap->extrapolate(data, error);

      datamat.print("Extrap Datamat");

      i++;
    }

  StateOutText s("scextest.ckpt");
  SavableState::save_state(extrap.pointer(),s);
  s.close();

  StateInText si("scextest.ckpt");
  Ref<SelfConsistentExtrapolation> e2;
  e2 << SavableState::restore_state(si);
  
  si.close();
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
