//
// object.cc
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
#include <util/misc/formio.h>
#include <util/render/render.h>
#include <util/render/object.h>

using namespace std;
using namespace sc;

static ClassDesc RenderedObject_cd(
  typeid(RenderedObject),"RenderedObject",1,"public DescribedClass",
  0, 0, 0);

RenderedObject::RenderedObject(const Ref<Material>& material):
  name_(),
  material_(material)
{
}

RenderedObject::RenderedObject(const Ref<KeyVal>& keyval)
{
  name_ = keyval->stringvalue("name");
  material_ << keyval->describedclassvalue("material");
  appearance_ << keyval->describedclassvalue("appearance");
  transform_ << keyval->describedclassvalue("transform");
}

RenderedObject::~RenderedObject()
{
}

void
RenderedObject::set_name(const char *name)
{
  if (name) name_ = name;
  else name_.resize(0);
}

void
RenderedObject::print(ostream& os) const
{
  os << "RenderedObject:" << endl;
  if (material_.nonnull()) {
      os << scprintf("  material = 0x%x\n", material_.pointer());
    }
  if (appearance_.nonnull()) {
      os << scprintf("  appearance = 0x%x\n", appearance_.pointer());
    }
  if (transform_.nonnull()) {
      os << scprintf("  transform = 0x%x\n", transform_.pointer());
    }
  os.flush();
}
  

static ClassDesc RenderedObjectSet_cd(
  typeid(RenderedObjectSet),"RenderedObjectSet",1,"public RenderedObject",
  0, create<RenderedObjectSet>, 0);

RenderedObjectSet::RenderedObjectSet(int capacity)
{
  capacity_ = capacity;
  n_ = 0;
  array_ = new Ref<RenderedObject>[capacity_];
}

RenderedObjectSet::RenderedObjectSet(const Ref<KeyVal>& keyval):
  RenderedObject(keyval)
{
  capacity_ = keyval->count("objects");
  if (keyval->error() != KeyVal::OK) {
      ExEnv::errn() << "RenderedObjectSet: error counting objects" << endl;
      abort();
    }
  n_ = capacity_;
  array_ = new Ref<RenderedObject>[capacity_];
  for (int i=0; i<n_; i++) {
      array_[i] << keyval->describedclassvalue("objects",i);
      if (keyval->error() != KeyVal::OK) {
          ExEnv::errn() << "RenderedObjectSet: error reading objects" << endl;
          abort();
        }
    }
}

RenderedObjectSet::~RenderedObjectSet()
{
  delete[] array_;
}

void
RenderedObjectSet::add(const Ref<RenderedObject>& object)
{
  if (capacity_ == n_) {
      capacity_ += 10;
      Ref<RenderedObject> *tmp = new Ref<RenderedObject>[capacity_];
      for (int i=0; i<n_; i++) {
          tmp[i] = array_[i];
        }
      delete[] array_;
      array_ = tmp;
    }
  array_[n_] = object;
  n_++;
}

void
RenderedObjectSet::render(const Ref<Render>& render)
{
  render->set(this);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
