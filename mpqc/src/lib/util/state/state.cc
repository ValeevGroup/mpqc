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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <util/misc/formio.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/stateio.h>

#define DEBUG 0

/////////////////////////////////////////////////////////////////

DescribedClass_REF_def(SavableState);

#define CLASSNAME SavableState
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SavableState::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

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
  RefSavableState th(this);
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
  RefKeyVal old_override;
  RefSavableState overriding_value;
  int p = si.push_key(keyword);
  const int can_override_objects = 0;
  if (can_override_objects && keyword && si.override().nonnull()) {
      overriding_value = si.override()->describedclassvalue(si.key());
      old_override = si.override();
      if (overriding_value.nonnull()) {
          si.set_override(0);
        }
    }
  // restore the pointer
  RefSavableState ss;
  if (objectname) si.dir_getobject(ss, objectname);
  else si.getobject(ss);
  if (overriding_value.nonnull()) {
      ExEnv::out() << node0 << indent
           << "overriding \"" << si.key() << "\": object of type ";
      if (ss.null()) ExEnv::out() << node0 << "(null)";
      else ExEnv::out() << node0 << ss->class_name();
      ExEnv::out() << node0 << " -> object of type "
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
  if (old_override.nonnull()) {
      si.set_override(old_override);
    }
  si.pop_key(p);
  return ret;
}

void
SavableState::save_object_state_(StateOut&so, const ClassDesc *cd)
{
  if (class_desc() != cd) {
      ExEnv::err() <<  "Warning:"
           << cd->name()
           << "::save_object_state: "
           << "exact type not known -- object not saved" << endl;
      return;
    }
  save_vbase_state(so);
  save_data_state(so);
}

void
SavableState::save_object_state(StateOut&)
{
  ExEnv::err() << "SavableState::save_object_state(StateOut&):" << endl
       << " only can be used when exact type is known" << endl
       << " otherwise use save_state(StateOut&)" << endl
       << " (object not saved)" << endl;
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

////////////////////////////////////////////////////////////////////////
// SSRefBase members

void
SSRefBase::save_data_state(StateOut&s)
{
  save_state(s);
}

void
SSRefBase::save_state(StateOut&so)
{
  SavableState::save_state(sspointer(),so);
}

void
SSRefBase::check_castdown_result(void* t, SavableState *ss,
                                 const ClassDesc *cd)
{
  if (!t && ss) {
      ExEnv::err() << node0
           << "SSRef::restore_state() got type \"" << ss->class_name()
           << "\""
           << " but expected \""
           << cd->name()
           << "\""
           << endl;
        abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
