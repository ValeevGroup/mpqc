//
// render.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <util/render/render.h>
#include <util/render/object.h>
#include <util/render/find.h>
#include <util/render/polygons.h>
#include <util/render/polysphere.h>

#define CLASSNAME Render
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
Render::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

Render::Render()
{
}

Render::Render(const RefKeyVal& keyval)
{
}

Render::~Render()
{
}

void
Render::push_material(const RefMaterial& m)
{
  material_stack_.push(m);
}

void
Render::push_appearance(const RefAppearance& a)
{
  appearance_stack_.push(a);
}

void
Render::push_transform(const RefTransform& t)
{
  transform_stack_.push(t);
}

RefMaterial
Render::pop_material()
{
  return material_stack_.pop();
}

RefAppearance
Render::pop_appearance()
{
  return appearance_stack_.pop();
}

RefTransform
Render::pop_transform()
{
  return transform_stack_.pop();
}

void
Render::render(const RefRenderedObject& object)
{
  if (object->material().nonnull()) push_material(object->material());
  if (object->transform().nonnull()) push_transform(object->transform());
  if (object->appearance().nonnull()) push_appearance(object->appearance());
  object->render(this);
  if (object->material().nonnull()) pop_material();
  if (object->transform().nonnull()) pop_transform();
  if (object->appearance().nonnull()) pop_appearance();
}

void
Render::set(const RefRenderedObjectSet& set)
{
  for (int i=0; i<set->n(); i++) {
      render(set->element(i));
    }
}

// This renders spheres by creating a RenderedPolygon object
void
Render::sphere(const RefRenderedSphere& sphere)
{
  // find the level of accuracy which should be used to render the sphere
  int level = 1;
  find_int_parameter_in_appearance_stack(appearance_stack_,
                                         &Appearance::level,
                                         level);
  RefRenderedPolygons poly(new RenderedPolygons);

  polysphere(level, poly.pointer());

  render(poly.pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
