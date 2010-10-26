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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>

using namespace sc;

void
sc::detail::FromStateIn<sc::RefSCMatrix>::get(sc::RefSCMatrix& t, StateIn& si, int& count) {
  RefSCDimension rowdim;  rowdim << SavableState::restore_state(si);
  RefSCDimension coldim;  coldim << SavableState::restore_state(si);
  const bool blockeddim = rowdim->blocks().nonnull();
  Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
  Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
  t = kit->matrix(rowdim,coldim);
  t.restore(si);
}

void
sc::detail::FromStateIn<sc::RefSymmSCMatrix>::get(sc::RefSymmSCMatrix& t, StateIn& si, int& count) {
  RefSCDimension dim;  dim << SavableState::restore_state(si);
  const bool blockeddim = dim->blocks().nonnull();
  Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
  Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
  t = kit->symmmatrix(dim);
  t.restore(si);
}

void
sc::detail::FromStateIn<sc::RefDiagSCMatrix>::get(sc::RefDiagSCMatrix& t, StateIn& si, int& count) {
  RefSCDimension dim;  dim << SavableState::restore_state(si);
  const bool blockeddim = dim->blocks().nonnull();
  Ref<SCMatrixKit> defaultkit = SCMatrixKit::default_matrixkit();
  Ref<SCMatrixKit> kit = blockeddim ? new BlockedSCMatrixKit(defaultkit) : defaultkit;
  t = kit->diagmatrix(dim);
  t.restore(si);
}

void
sc::detail::ToStateOut<sc::RefSCMatrix>::put(const sc::RefSCMatrix& t, StateOut& so, int& count) {
  SavableState::save_state(t.rowdim().pointer(),so);
  SavableState::save_state(t.coldim().pointer(),so);
  (const_cast<SCMatrix*>(t.pointer()))->save(so);
}

void
sc::detail::ToStateOut<sc::RefSymmSCMatrix>::put(const sc::RefSymmSCMatrix& t, StateOut& so, int& count) {
  SavableState::save_state(t.dim().pointer(),so);
  (const_cast<SymmSCMatrix*>(t.pointer()))->save(so);
}

void
sc::detail::ToStateOut<sc::RefDiagSCMatrix>::put(const sc::RefDiagSCMatrix& t, StateOut& so, int& count) {
  SavableState::save_state(t.dim().pointer(),so);
  (const_cast<DiagSCMatrix*>(t.pointer()))->save(so);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
