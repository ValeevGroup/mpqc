//
// ref.cc --- implementation of the reference counting classes
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

#include <util/ref/ref.h>
#include <util/misc/exenv.h>

using namespace std;

void
RefCount::error(const char * w) const
{
  ExEnv::err() << "RefCount: ERROR: " << w << endl;
  abort();
}

void
RefCount::too_many_refs() const
{
  error("Too many refs.");
}

void
RefCount::not_enough_refs() const
{
  error("Ref count dropped below zero.");
}

#if REF_CHECKSUM
void
RefCount::bad_checksum() const
{
  error("Bad checksum.");
}
#endif

RefCount::~RefCount()
{
#if REF_CHECKSUM
  if (_reference_count_ == 0 && _checksum_ == 0) {
      error("Deleting a deleted or overwritten object.");
    }
#endif
#if REF_MANAGE
  if (managed() && nreference()) {
      error("Deleting a referenced object.");
    }
#endif
#if REF_CHECKSUM
  _reference_count_ = 0;
  _checksum_ = 0;
#endif
}

///////////////////////////////////////////////////////////////////////

void
RefBase::warn ( const char * msg) const
{
  ExEnv::err() << "WARNING: " << msg << endl;
}
void
RefBase::warn_ref_to_stack() const
{
  warn("Ref: creating a reference to stack data");
}
void
RefBase::warn_skip_stack_delete() const
{
  warn("Ref: skipping delete of object on the stack");
}
void
RefBase::warn_bad_ref_count() const
{
  warn("Ref: bad reference count in referenced object\n");
}
void
RefBase::ref_info(RefCount*p, ostream& os) const
{
  if (p)
      os << "nreference() = " << p->nreference() << endl;
  else
      os << "reference is null" << endl;
}

void
RefBase::require_nonnull() const
{
  if (parentpointer() == 0) {
      ExEnv::err() << "RefBase: needed a nonnull pointer but got null"
           << endl;
      abort();
    }
}

RefBase::~RefBase()
{
}

void
RefBase::check_pointer() const
{
  if (parentpointer() && parentpointer()->nreference() <= 0) {
      warn_bad_ref_count();
    }
}

void
RefBase::ref_info(ostream& os) const
{
  RefBase::ref_info(parentpointer(),os);
}

void
RefBase::reference(RefCount *p)
{
  if (p) {
#if REF_CHECK_STACK
      if (DO_REF_CHECK_STACK(p)) {
          DO_REF_UNMANAGE(p);
          warn_ref_to_stack();
        }
#endif
      p->reference();
    }
}

void
RefBase::dereference(RefCount *p)
{
  if (p && p->dereference()<=0) {
      delete p;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
