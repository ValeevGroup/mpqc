//
// object.h
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

#ifndef _util_render_object_h
#define _util_render_object_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <util/keyval/keyval.h>
#include <util/render/material.h>
#include <util/render/appearance.h>
#include <util/render/transform.h>

DescribedClass_REF_fwddec(Render);

class RenderedObject: public DescribedClass {
#   define CLASSNAME RenderedObject
#   include <util/class/classda.h>
  protected:
    char* name_;
    RefMaterial material_;
    RefAppearance appearance_;
    RefTransform transform_;

    friend class Render;
  public:
    RenderedObject(const RefMaterial& = 0);
    RenderedObject(const RefKeyVal&);
    ~RenderedObject();
    const char* name() const { return name_; }
    void set_name(const char *);
    RefMaterial material() const { return material_; }
    RefAppearance appearance() const { return appearance_; }
    RefTransform transform() const { return transform_; }
    void material(const RefMaterial&m) { material_ = m; }
    void appearance(const RefAppearance&a) { appearance_ = a; }
    void transform(const RefTransform&t) { transform_ = t; }

    virtual void print(ostream& = ExEnv::out()) const;

    // to be called only by derivatives of Render
    virtual void render(const RefRender&) = 0;
};
DescribedClass_REF_dec(RenderedObject);

class RenderedObjectSet: public RenderedObject {
#   define CLASSNAME RenderedObjectSet
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    int capacity_;
    int n_;
    RefRenderedObject *array_;
  protected:
    void render(const RefRender&);
  public:
    RenderedObjectSet(int capacity = 10);
    RenderedObjectSet(const RefKeyVal&);
    ~RenderedObjectSet();
    int n() const { return n_; }
    const RefRenderedObject& element(int i) const { return array_[i]; }
    virtual void add(const RefRenderedObject&);
};
DescribedClass_REF_dec(RenderedObjectSet);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
