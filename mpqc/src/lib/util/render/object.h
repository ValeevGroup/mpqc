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

#include <iostream>

#include <util/keyval/keyval.h>
#include <util/render/material.h>
#include <util/render/appearance.h>
#include <util/render/transform.h>

namespace sc {

class Render;

class RenderedObject: public DescribedClass {
  protected:
    std::string name_;
    Ref<Material> material_;
    Ref<Appearance> appearance_;
    Ref<Transform> transform_;

    friend class Render;
  public:
    RenderedObject(const Ref<Material>& = 0);
    RenderedObject(const Ref<KeyVal>&);
    ~RenderedObject();
    const char* name() const { return name_.c_str(); }
    void set_name(const char *);
    Ref<Material> material() const { return material_; }
    Ref<Appearance> appearance() const { return appearance_; }
    Ref<Transform> transform() const { return transform_; }
    void material(const Ref<Material>&m) { material_ = m; }
    void appearance(const Ref<Appearance>&a) { appearance_ = a; }
    void transform(const Ref<Transform>&t) { transform_ = t; }

    virtual void print(std::ostream& = ExEnv::out0()) const;

    // to be called only by derivatives of Render
    virtual void render(const Ref<Render>&) = 0;
};


class RenderedObjectSet: public RenderedObject {
  private:
    int capacity_;
    int n_;
    Ref<RenderedObject> *array_;
  protected:
    void render(const Ref<Render>&);
  public:
    RenderedObjectSet(int capacity = 10);
    RenderedObjectSet(const Ref<KeyVal>&);
    ~RenderedObjectSet();
    int n() const { return n_; }
    const Ref<RenderedObject>& element(int i) const { return array_[i]; }
    virtual void add(const Ref<RenderedObject>&);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
