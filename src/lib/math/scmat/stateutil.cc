//
// stateutil.cc
//
// Copyright (C) 2009 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>

using namespace sc;

template<>
void
sc::FromStateIn<sc::RefSCVector>(sc::RefSCVector& t, StateIn& si, int& count) {
  int nonnull;
  count += si.get(nonnull);
  if (nonnull) {
    RefSCDimension dim;  dim << SavableState::restore_state(si);
    const bool blockeddim = dim->blocks()->subdim(0).nonnull();  // has subdims? use blocked kit
    Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
    Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
    t = kit->vector(dim);
    t.restore(si);
  }
}

template<>
void
sc::FromStateIn<sc::RefSCMatrix>(sc::RefSCMatrix& t, StateIn& si, int& count) {
  int nonnull;
  count += si.get(nonnull);
  if (nonnull) {
    RefSCDimension rowdim;  rowdim << SavableState::restore_state(si);
    RefSCDimension coldim;  coldim << SavableState::restore_state(si);
    const bool blockeddim = rowdim->blocks()->subdim(0).nonnull()
                         && coldim->blocks()->subdim(0).nonnull();  // has subdims? use blocked kit
    Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
    Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
    t = kit->matrix(rowdim,coldim);
    t.restore(si);
  }
}

template<>
void
sc::FromStateIn<sc::RefSymmSCMatrix>(sc::RefSymmSCMatrix& t, StateIn& si, int& count) {
  int nonnull;
  count += si.get(nonnull);
  if (nonnull) {
    RefSCDimension dim;  dim << SavableState::restore_state(si);
    const bool blockeddim = dim->blocks()->subdim(0).nonnull();  // has subdims? use blocked kit
    Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
    Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
    t = kit->symmmatrix(dim);
    t.restore(si);
  }
}

template<>
void
sc::FromStateIn<sc::RefDiagSCMatrix>(sc::RefDiagSCMatrix& t, StateIn& si, int& count) {
  int nonnull;
  count += si.get(nonnull);
  if (nonnull) {
    RefSCDimension dim;  dim << SavableState::restore_state(si);
    const bool blockeddim = dim->blocks()->subdim(0).nonnull();  // has subdims? use blocked kit
    Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
    Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
    t = kit->diagmatrix(dim);
    t.restore(si);
  }
}

template<>
void
sc::ToStateOut<sc::RefSCVector>(const sc::RefSCVector& t, StateOut& so, int& count) {
  if (t.null())
    count += so.put(0);
  else {
    count += so.put(1);
    SavableState::save_state(t.dim().pointer(),so);
    const_cast<RefSCVector&>(t).save(so);
  }
}

template<>
void
sc::ToStateOut<sc::RefSCMatrix>(const sc::RefSCMatrix& t, StateOut& so, int& count) {
  if (t.null())
    count += so.put(0);
  else {
    count += so.put(1);
    SavableState::save_state(t.rowdim().pointer(),so);
    SavableState::save_state(t.coldim().pointer(),so);
    const_cast<RefSCMatrix&>(t).save(so);
  }
}

template<>
void
sc::ToStateOut<sc::RefSymmSCMatrix>(const sc::RefSymmSCMatrix& t, StateOut& so, int& count) {
  if (t.null())
    count += so.put(0);
  else {
    count += so.put(1);
    SavableState::save_state(t.dim().pointer(),so);
    const_cast<RefSymmSCMatrix&>(t).save(so);
  }
}

template<>
void
sc::ToStateOut<sc::RefDiagSCMatrix>(const sc::RefDiagSCMatrix& t, StateOut& so, int& count) {
  if (t.null())
    count += so.put(0);
  else {
    count += so.put(1);
    SavableState::save_state(t.dim().pointer(),so);
    const_cast<RefDiagSCMatrix&>(t).save(so);
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
