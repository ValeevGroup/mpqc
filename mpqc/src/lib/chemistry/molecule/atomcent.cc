//
// atomcent.cc
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

#include <string.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/atomcent.h>

DescribedClass_REF_def(AtomicCenter);

#define CLASSNAME AtomicCenter
#define PARENTS public SavableState
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
AtomicCenter::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AtomicCenter::AtomicCenter() :
  label_(0)
{
}

AtomicCenter::AtomicCenter(const AtomicCenter&ac) :
  p(ac.p),element_(ac.element_),label_(0)
{
  if (ac.label_) {
    label_ = new char[strlen(ac.label_)+1];
    strcpy(label_,ac.label_);
  }
}

AtomicCenter::AtomicCenter(const char*symbol,double x,double y,double z, 
                           const char *lab) :
  p(x,y,z),
  element_(symbol),
  label_(0)
{
  if (lab) {
    label_ = new char[strlen(lab)+1];
    strcpy(label_,lab);
  }
}

AtomicCenter::~AtomicCenter()
{
  if (label_) delete[] label_; label_=0;
}

AtomicCenter& AtomicCenter::operator=(const AtomicCenter&ac)
{
  p = ac.p;
  element_ = ac.element_;
  if (ac.label_) {
    label_ = new char[strlen(ac.label_)+1];
    strcpy(label_,ac.label_);
  }

  return *this;
}

void AtomicCenter::save_data_state(StateOut& so)
{
  p.save_object_state(so);
  element_.save_object_state(so);
  so.putstring(label_);
}

AtomicCenter::AtomicCenter(StateIn& si):
  SavableState(si),
  p(si),
  element_(si)
{
  si.getstring(label_);
}

void AtomicCenter::print(ostream& os)
{
  os << node0 << indent
     << scprintf("%2s",element().symbol());
  point().print(os);
}

double
AtomicCenter::mass() const
{
  return element_.mass();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
