//
// render.h
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

#ifndef _util_render_render_h
#define _util_render_render_h

#include <util/class/class.h>
#include <util/render/appearance.h>
#include <util/render/material.h>
#include <util/render/transform.h>
#include <util/render/stack.h>

namespace sc {

class RenderedObject;
class AnimatedObject;
class RenderedObjectSet;
class RenderedSphere;
class RenderedPolygons;
class RenderedPolylines;

class Render: public DescribedClass {
  protected:
    Ref<Material> default_material_;
    Ref<Appearance> default_appearance_;
    Ref<Transform> default_transform_;

    Stack<Ref<Material> > material_stack_;
    Stack<Ref<Appearance> > appearance_stack_;
    Stack<Ref<Transform> > transform_stack_;

    virtual void push_material(const Ref<Material>& m);
    virtual void push_appearance(const Ref<Appearance>& a);
    virtual void push_transform(const Ref<Transform>& t);
    virtual Ref<Material> pop_material();
    virtual Ref<Appearance> pop_appearance();
    virtual Ref<Transform> pop_transform();

  public:
    Render();
    Render(const Ref<KeyVal>&);
    virtual ~Render();

    Ref<Material> default_material() { return default_material_; }
    Ref<Appearance> default_appearance() { return default_appearance_; }
    Ref<Transform> default_transform() { return default_transform_; }
    void default_material(const Ref<Material>& m) { default_material_ = m; }
    void default_appearance(const Ref<Appearance>& a) {default_appearance_ = a;}
    void default_transform(const Ref<Transform>& t) {default_transform_ = t;}

    virtual void clear() = 0;

    virtual void render(const Ref<RenderedObject>&);
    virtual void animate(const Ref<AnimatedObject> &);

    virtual void set(const Ref<RenderedObjectSet>&);
    virtual void sphere(const Ref<RenderedSphere>&);
    virtual void polygons(const Ref<RenderedPolygons>&) = 0;
    virtual void polylines(const Ref<RenderedPolylines>&) = 0;
};


class FileRender: public Render {
  protected:
    char* filename_;
    char* basename_;
    std::streambuf *sbuf_;
    int delete_sbuf_;
    int depth_;

    char *get_filename(const char *objectname);
    void open_sbuf(const char *objectname);
    void close_sbuf();
  public:
    FileRender(const char * filename);
    FileRender(std::ostream &o = ExEnv::out0());
    FileRender(const Ref<KeyVal>&);
    virtual ~FileRender();

    void clear();

    virtual void set_filename(const char *name);
    virtual void set_basename(const char *name);
    virtual const char *file_extension();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
