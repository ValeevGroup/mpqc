//
// ansatz.cc
//
// Copyright (C) 2006 Edward Valeev
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

#include <chemistry/qc/mbptr12/ansatz.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/class/scexception.h>

using namespace sc;

static ClassDesc LinearR12Ansatz_cd(
  typeid(LinearR12Ansatz),"LinearR12Ansatz",3,"virtual public SavableState",
  create<LinearR12Ansatz>, create<LinearR12Ansatz>, create<LinearR12Ansatz>);

LinearR12Ansatz::LinearR12Ansatz() :
  projector_(LinearR12::Projector_2),
  diag_(false),
  amplitudes_(LinearR12::GeminalAmplitudeAnsatz_fullopt),
  wof_(false),
  orbital_product_(LinearR12::OrbProd_ij) {
}

LinearR12Ansatz::LinearR12Ansatz(const Ref<KeyVal>& keyval)
{
  projector_ = (LinearR12::Projector)keyval->intvalue("projector",KeyValValueint(2));
  const bool default_wof = (projector_ == LinearR12::Projector_0);
  wof_ = keyval->booleanvalue("wof",KeyValValueboolean((int)default_wof));

  diag_ = keyval->booleanvalue("diag",KeyValValueboolean((int)false));

  // amplitudes should be fixed by default if the diagonal ansatz is used
  const std::string default_amplitudes_str = diag_ ? std::string("fixed") : std::string("optimized");
  const std::string amplitudes_str = keyval->stringvalue("amplitudes",
                                                         KeyValValuestring(default_amplitudes_str));
  if (amplitudes_str == std::string("fixed"))
    amplitudes_ = LinearR12::GeminalAmplitudeAnsatz_fixed;
  else if (amplitudes_str == std::string("optimized"))
    amplitudes_ = LinearR12::GeminalAmplitudeAnsatz_fullopt;
  else if (amplitudes_str == std::string("scaledfixed"))
    amplitudes_ = LinearR12::GeminalAmplitudeAnsatz_scaledfixed;
  else
    throw InputError("Invalid value for keyword \"amplitudes\"",__FILE__,__LINE__);
  if ( (diag_==false) && (amplitudes_!=LinearR12::GeminalAmplitudeAnsatz_fullopt) ){
    throw InputError("LinearR12Ansatz::LinearR12Ansatz -- amplitudes can only be fixed if diag is true",__FILE__,__LINE__);
  }

  std::string op = keyval->stringvalue("orbital_product",KeyValValuestring("ij"));
  if (op == "ij")
    orbital_product_ = LinearR12::OrbProd_ij;
  else if (op == "pq")
    orbital_product_ = LinearR12::OrbProd_pq;
  else
    throw InputError("LinearR12Ansatz::LinearR12Ansatz -- invalid value for orbital_product",__FILE__,__LINE__);
  if (orbital_product_ != LinearR12::OrbProd_ij && diag_) {
    throw InputError("LinearR12Ansatz::LinearR12Ansatz -- diagonal ansatz only allowed when orbital_product = ij",__FILE__,__LINE__);
  }
}

LinearR12Ansatz::LinearR12Ansatz(StateIn& s) :
  SavableState(s)
{
  int p; s.get(p); projector_ = (LinearR12::Projector)p;
  int d; s.get(d); diag_ = (bool)d;
  int a; s.get(a); amplitudes_ = (LinearR12::GeminalAmplitudeAnsatz)a;
  if (s.version(::class_desc<LinearR12Ansatz>()) >= 3) {
    int w; s.get(w); wof_ = (bool)w;
  }
  if (s.version(::class_desc<LinearR12Ansatz>()) >= 2) {
    int o; s.get(o); orbital_product_ = (LinearR12::OrbitalProduct)o;
  }
}

LinearR12Ansatz::~LinearR12Ansatz() {}

void
LinearR12Ansatz::save_data_state(StateOut& s)
{
  s.put((int)projector_);
  s.put((int)diag_);
  s.put((int)amplitudes_);
  s.put((int)wof_);
  s.put((int)orbital_product_);
}

void
LinearR12Ansatz::print(std::ostream& o) const
{
  o << indent << "LinearR12Ansatz:" << std::endl;
  o << incindent;

  o << indent << "Orbital Product Space: ";
  switch(orbital_product_) {
    case LinearR12::OrbProd_ij: o << "ij"; break;
    case LinearR12::OrbProd_pq: o << "pq"; break;
  }
  o << std::endl;

  o << indent << "Projector: ";
  switch(projector_) {
    case LinearR12::Projector_0: o << "0  , i.e. 1"; break;
    case LinearR12::Projector_1: o << "1  , i.e. (1-P1)(1-P2)"; break;
    case LinearR12::Projector_2: o << "2  , i.e. (1-O1)(1-O2)-V1V2"; break;
    case LinearR12::Projector_3: o << "3  , i.e. 1-P1P2"; break;
  }
  o << std::endl;

  std::string amplitudes_str;
  switch (amplitudes_) {
    case LinearR12::GeminalAmplitudeAnsatz_fullopt:
      amplitudes_str = std::string("optimized");
    case LinearR12::GeminalAmplitudeAnsatz_fixed:
      amplitudes_str = std::string("fixed");
    case LinearR12::GeminalAmplitudeAnsatz_scaledfixed:
      amplitudes_str = std::string("scaled fixed");
  }
  o << indent << "Ansatz: " << (diag_ ? "diagonal" : "orbital-invariant")
              << " with " << amplitudes_str << " amplitudes" << std::endl;
  o << indent << "WOF: " << (wof_ ? "true" : "false") << std::endl;
  o << decindent;
}

LinearR12::Projector
LinearR12Ansatz::projector() const { return projector_; }

bool
LinearR12Ansatz::diag() const { return diag_; }

LinearR12::GeminalAmplitudeAnsatz
LinearR12Ansatz::amplitudes() const { return amplitudes_; }

bool
LinearR12Ansatz::wof() const { return wof_; }

LinearR12::OrbitalProduct
LinearR12Ansatz::orbital_product() const { return orbital_product_; }
