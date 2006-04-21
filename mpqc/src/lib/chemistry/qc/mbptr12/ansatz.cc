//
// ansatz.cc
//
// Copyright (C) 2006 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <chemistry/qc/mbptr12/ansatz.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>

using namespace sc;

static ClassDesc LinearR12Ansatz_cd(
  typeid(LinearR12Ansatz),"LinearR12Ansatz",1,"virtual public SavableState",
  create<LinearR12Ansatz>, create<LinearR12Ansatz>, create<LinearR12Ansatz>);

LinearR12Ansatz::LinearR12Ansatz() : projector_(LinearR12::Projector_2), diag_(false) {}

LinearR12Ansatz::LinearR12Ansatz(const Ref<KeyVal>& keyval)
{
  projector_ = (LinearR12::Projector)keyval->intvalue("projector",KeyValValueint(2));
  diag_ = keyval->booleanvalue("diag",KeyValValueboolean((int)false));
}

LinearR12Ansatz::LinearR12Ansatz(StateIn& s) :
  SavableState(s)
{
  int p; s.get(p); projector_ = (LinearR12::Projector)p;
  int d; s.get(d); diag_ = (bool)d;
}

LinearR12Ansatz::~LinearR12Ansatz() {}

void
LinearR12Ansatz::save_data_state(StateOut& s)
{
  s.put((int)projector_);
  s.put((int)diag_);
}

void
LinearR12Ansatz::print(std::ostream& o) const
{
  o << indent << "LinearR12Ansatz:" << std::endl;
  o << incindent;

  o << indent << "Projector: ";
  switch(projector_) {
    case LinearR12::Projector_2: o << "2"; break;
    case LinearR12::Projector_3: o << "3"; break;
  }
  o << std::endl;
  
  o << indent << "Ansatz: " << (diag_ ? "diagonal" : "orbital-invariant") << std::endl;
  o << decindent;
}

LinearR12::Projector
LinearR12Ansatz::projector() const { return projector_; }

bool
LinearR12Ansatz::diag() const { return diag_; }

