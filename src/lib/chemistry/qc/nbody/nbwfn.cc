//
// nbwfn.cc
//
// Copyright (C) 2012 Edward Valeev
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

#include <chemistry/qc/nbody/nbwfn.h>

using namespace sc;

ClassDesc
ManyBodyWavefunction::class_desc_(typeid(ManyBodyWavefunction),
                     "ManyBodyWavefunction",
                     1,               // version
                     "public Wavefunction", // must match parent
                     0,               // change to create<ManyBodyWavefunction> if this class is DefaultConstructible
                     0, // change to 0 if this class is not KeyValConstructible
                     0  // change to 0 if this class is not StateInConstructible
                     );

ManyBodyWavefunction::ManyBodyWavefunction(const Ref<KeyVal>& keyval) : Wavefunction(keyval) {

  // if world not given, make this the center of a new World
  world_ << keyval->describedclassvalue("world", KeyValValueRefDescribedClass(0));
  if (world_ == 0)
    world_ = new WavefunctionWorld(keyval);
  if (world_ == 0)
    throw InputError("CI requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
                     __FILE__, __LINE__, "world");
  if (world_->wfn() == 0) world_->set_wfn(this);

  // tell reference to live in my world, unless user said otherwise (probably a bad idea)
  if (not keyval->exists("reference:world")) {
    Ref<KeyVal> ref_kv = new PrefixKeyVal(keyval, "reference");
    if (ref_kv == 0)
      throw InputError("ManyBodyWavefunction::ManyBodyWavefunction -- reference keyword missing",
                       __FILE__, __LINE__, "reference", "0", this->class_desc());

    const char* ref_classname = keyval->classname("reference");
    if (ref_classname == 0)
      throw InputError("ManyBodyWavefunction::ManyBodyWavefunction -- invalid reference specification",
                       __FILE__, __LINE__, "reference", "0");

    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("world", world_.pointer());
    Ref<KeyVal> newref_kv = new AggregateKeyVal(akv, ref_kv);

    refwfn_ << newref_kv->describedclass(ref_classname);
  }
  else
    refwfn_ << keyval->describedclassvalue("reference", KeyValValueRefDescribedClass(0));
  if (refwfn_ == 0)
    throw InputError("ManyBodyWavefunction::ManyBodyWavefunction -- invalid reference object given",
                     __FILE__, __LINE__, "reference", "0", this->class_desc());
  if (refwfn_->world() != world_)
    throw InputError("ManyBodyWavefunction and its reference should live in the same WavefunctionWorld",
                     __FILE__,__LINE__);

  this->set_desired_value_accuracy(desired_value_accuracy());
}

ManyBodyWavefunction::ManyBodyWavefunction(StateIn& si) : Wavefunction(si) {
  world_ << SavableState::restore_state(si);
  refwfn_ << SavableState::restore_state(si);
}

ManyBodyWavefunction::~ManyBodyWavefunction() {
  // this may be necessary if this is a templated class
  const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
}

void
ManyBodyWavefunction::save_data_state(StateOut& so) {
  SavableState::save_state(world_.pointer(), so);
  SavableState::save_state(refwfn_.pointer(), so);
}

void ManyBodyWavefunction::purge() {
  MolecularEnergy::purge();
  refwfn_->purge();
  world_->obsolete();
}

void ManyBodyWavefunction::obsolete() {
  Function::obsolete();
  refwfn_->obsolete();
  world_->obsolete();
}

void ManyBodyWavefunction::symmetry_changed() {
  Wavefunction::symmetry_changed();
  //refwfn_->symmetry_changed();
}

void sc::ManyBodyWavefunction::set_desired_value_accuracy(double acc) {
  Function::set_desired_value_accuracy(acc);
  refwfn_->set_desired_value_accuracy(acc * ref_to_corr_acc());
}

void
sc::ManyBodyWavefunction::print(std::ostream&o) const {
  using std::endl;
  o << indent << "ManyBodyWavefunction:" << endl;
  o << incindent;
  Wavefunction::print(o);
  o << indent << "Reference Wavefunction:" << endl;
  o << incindent; refwfn_->print(o); o << decindent << endl;
  o << decindent;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
