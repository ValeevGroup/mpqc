//
// state.cc
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

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <util/misc/formio.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/stateio.h>

using namespace std;
using namespace sc;

#define DEBUG 0

/////////////////////////////////////////////////////////////////

static ClassDesc SavableState_cd(
    typeid(SavableState),"SavableState",1,"public DescribedClass");

SavableState::SavableState()
{
}

SavableState::SavableState(const SavableState&)
{
}

SavableState& SavableState::operator=(const SavableState&)
{
  return *this;
}

SavableState::SavableState(StateIn&si)
{
  // In case si is looking for the next pointer, let it know i
  // have one.
  reference();
  Ref<SavableState> th(this);
  si.haveobject(th);
  th = 0;
  dereference();

  // The following gets the version of this class and all of the
  // parent classes.  This is only needed for restoring objects
  // that were saved with save_object_state and don't necessarily
  // have all of their version information already restored.
  if (si.need_classdesc()) {
      const ClassDesc* tcd;
      si.get(&tcd);
    }
}

SavableState::~SavableState()
{
}

void
SavableState::save_state(StateOut&so)
{
  save_state(this,so);
}

void
SavableState::save_state(SavableState*th,StateOut&so)
{
  so.putobject(th);
}

SavableState*
SavableState::restore_state(StateIn& si)
{
  return dir_restore_state(si,0,0);
}

SavableState*
SavableState::key_restore_state(StateIn& si,
                                const char *keyword)
{
  return dir_restore_state(si,0,keyword);
}

SavableState*
SavableState::dir_restore_state(StateIn&si, const char *objectname,
                                const char *keyword)
{
  Ref<KeyVal> old_override;
  Ref<SavableState> overriding_value;
  int p = si.push_key(keyword);
  const int can_override_objects = 0;
  if (can_override_objects && keyword && si.override()) {
      overriding_value << si.override()->describedclassvalue(si.key());
      old_override = si.override();
      if (overriding_value) {
          si.set_override(0);
        }
    }
  // restore the pointer
  Ref<SavableState> ss;
  if (objectname) si.dir_getobject(ss, objectname);
  else si.getobject(ss);
  if (overriding_value) {
      ExEnv::out0() << indent
           << "overriding \"" << si.key() << "\": object of type ";
      if (ss == 0) ExEnv::out0() << "(null)";
      else ExEnv::out0() << ss->class_name();
      ExEnv::out0() << " -> object of type "
           << overriding_value->class_name()
           << endl;
      ss = overriding_value;
    }
  SavableState *ret = ss.pointer();
  if (ret) {
      ret->reference();
      ss = 0;
      ret->dereference();
    }
  if (old_override) {
      si.set_override(old_override);
    }
  si.pop_key(p);
  return ret;
}

void
SavableState::save_object_state(StateOut&so)
{
  save_vbase_state(so);
  save_data_state(so);
}

void
SavableState::save_vbase_state(StateOut&so)
{
  SavableState::save_data_state(so);
}
void
SavableState::save_data_state(StateOut& so)
{
  if (so.need_classdesc()) so.put(class_desc());
}

/////////////////////////////////////////////////////////////////////////////

DummySavableState::DummySavableState() {}

DummySavableState::DummySavableState(StateIn& si) {
}

void
DummySavableState::save_data_state(StateOut& so) {
}

ClassDesc
DummySavableState::class_desc_(typeid(this_type),
                               "DummySavableState",
                               1,
                               "virtual public SavableState",
                               create<this_type>,
                               0,
                               create<this_type> );

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
