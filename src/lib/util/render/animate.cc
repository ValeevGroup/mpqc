//
// animate.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#include <util/misc/formio.h>
#include <util/render/animate.h>
#include <util/render/object.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// AnimatedObject

static ClassDesc AnimatedObject_cd(
  typeid(AnimatedObject),"AnimatedObject",1,"public DescribedClass",
  0, 0, 0);

AnimatedObject::AnimatedObject()
{
}

AnimatedObject::AnimatedObject(const Ref<KeyVal>& keyval)
{
  name_ = keyval->stringvalue("name");
}

AnimatedObject::~AnimatedObject()
{
}

void
AnimatedObject::set_name(const char *name)
{
  if (name) name_ = name;
  else name_.resize(0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
