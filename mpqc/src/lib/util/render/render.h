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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_render_render_h
#define _util_render_render_h

#include <util/class/class.h>
#include <util/render/appearance.h>
#include <util/render/material.h>
#include <util/render/transform.h>
#include <util/render/stack.h>

DescribedClass_REF_fwddec(AnimatedObject);
DescribedClass_REF_fwddec(RenderedObject);
DescribedClass_REF_fwddec(RenderedObjectSet);
DescribedClass_REF_fwddec(RenderedSphere);
DescribedClass_REF_fwddec(RenderedPolygons);
DescribedClass_REF_fwddec(RenderedPolylines);

class Render: public DescribedClass {
#   define CLASSNAME Render
#   include <util/class/classda.h>
  protected:
    RefMaterial default_material_;
    RefAppearance default_appearance_;
    RefTransform default_transform_;

    Stack<RefMaterial> material_stack_;
    Stack<RefAppearance> appearance_stack_;
    Stack<RefTransform> transform_stack_;

    virtual void push_material(const RefMaterial& m);
    virtual void push_appearance(const RefAppearance& a);
    virtual void push_transform(const RefTransform& t);
    virtual RefMaterial pop_material();
    virtual RefAppearance pop_appearance();
    virtual RefTransform pop_transform();

  public:
    Render();
    Render(const RefKeyVal&);
    virtual ~Render();

    RefMaterial default_material() { return default_material_; }
    RefAppearance default_appearance() { return default_appearance_; }
    RefTransform default_transform() { return default_transform_; }
    void default_material(const RefMaterial& m) { default_material_ = m; }
    void default_appearance(const RefAppearance& a) {default_appearance_ = a;}
    void default_transform(const RefTransform& t) {default_transform_ = t;}

    virtual void clear() = 0;

    virtual void render(const RefRenderedObject&);
    virtual void animate(const RefAnimatedObject &);

    virtual void set(const RefRenderedObjectSet&);
    virtual void sphere(const RefRenderedSphere&);
    virtual void polygons(const RefRenderedPolygons&) = 0;
    virtual void polylines(const RefRenderedPolylines&) = 0;
};
DescribedClass_REF_dec(Render);

class FileRender: public Render {
#   define CLASSNAME FileRender
#   include <util/class/classd.h>
  protected:
    char* filename_;
    char* basename_;
    streambuf *sbuf_;
    int delete_sbuf_;
    int depth_;

    char *get_filename(const char *objectname);
    void open_sbuf(const char *objectname);
    void close_sbuf();
  public:
    FileRender(const char * filename);
    FileRender(ostream &o = ExEnv::out());
    FileRender(const RefKeyVal&);
    virtual ~FileRender();

    void clear();

    virtual void set_filename(const char *name);
    virtual void set_basename(const char *name);
    virtual const char *file_extension();
};
DescribedClass_REF_dec(FileRender);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
