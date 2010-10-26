//
// sphere.cc
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

#include <util/render/render.h>
#include <util/render/object.h>
#include <util/render/sphere.h>

using namespace sc;

static ClassDesc RenderedSphere_cd(
  typeid(RenderedSphere),"RenderedSphere",1,"public RenderedObject",
  0, create<RenderedSphere>, 0);

RenderedSphere::RenderedSphere()
{
}

RenderedSphere::RenderedSphere(const Ref<Material>& material):
  RenderedObject(material)
{
}

RenderedSphere::RenderedSphere(const Ref<KeyVal>& keyval):
  RenderedObject(keyval)
{
}

RenderedSphere::~RenderedSphere()
{
}

void
RenderedSphere::render(const Ref<Render>& render)
{
  render->sphere(this);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
