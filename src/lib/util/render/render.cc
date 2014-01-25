//
// render.cc
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

#include <util/misc/formio.h>

#include <util/render/render.h>
#include <util/render/object.h>
#include <util/render/animate.h>
#include <util/render/find.h>
#include <util/render/polygons.h>
#include <util/render/polysphere.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// Render

static ClassDesc Render_cd(
  typeid(Render),"Render",1,"public DescribedClass",
  0, 0, 0);

Render::Render()
{
}

Render::Render(const Ref<KeyVal>& keyval)
{
}

Render::~Render()
{
}

void
Render::push_material(const Ref<Material>& m)
{
  material_stack_.push(m);
}

void
Render::push_appearance(const Ref<Appearance>& a)
{
  appearance_stack_.push(a);
}

void
Render::push_transform(const Ref<Transform>& t)
{
  transform_stack_.push(t);
}

Ref<Material>
Render::pop_material()
{
  return material_stack_.pop();
}

Ref<Appearance>
Render::pop_appearance()
{
  return appearance_stack_.pop();
}

Ref<Transform>
Render::pop_transform()
{
  return transform_stack_.pop();
}

void
Render::render(const Ref<RenderedObject>& object)
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
Render::set(const Ref<RenderedObjectSet>& set)
{
  for (int i=0; i<set->n(); i++) {
      render(set->element(i));
    }
}

// This renders spheres by creating a RenderedPolygon object
void
Render::sphere(const Ref<RenderedSphere>& sphere)
{
  // find the level of accuracy which should be used to render the sphere
  int level = 1;
  find_int_parameter_in_appearance_stack(appearance_stack_,
                                         &Appearance::level,
                                         level);
  Ref<RenderedPolygons> poly(new RenderedPolygons);

  polysphere(level, poly.pointer());

  render(poly.pointer());
}

void
Render::animate(const Ref<AnimatedObject> &animated_object)
{
  for (int i=0; i<animated_object->nobject(); i++) {
      Ref<RenderedObject> object = animated_object->object(i);
      render(object);
    }
}

/////////////////////////////////////////////////////////////////////////////
// FileRender

static ClassDesc FileRender_cd(
  typeid(FileRender),"FileRender",1,"public Render",
  0, 0, 0);

FileRender::FileRender(const char * filename)
{
  filename_ = 0;
  basename_ = 0;
  depth_ = 0;
  sbuf_ = 0;
  delete_sbuf_ = 0;
  set_filename(filename);
}

FileRender::FileRender(ostream &o)
{
  filename_ = 0;
  basename_ = 0;
  depth_ = 0;
  sbuf_ = o.rdbuf();
  delete_sbuf_ = 0;
}

FileRender::FileRender(const Ref<KeyVal>& keyval):
  Render(keyval)
{
  std::string filename = keyval->stringvalue("filename");
  std::string basename = keyval->stringvalue("basename");
  if (filename.empty() && basename.empty()) {
      const char *cbasename = SCFormIO::default_basename();
      if (cbasename) {
          basename = cbasename;
        }
    }
  depth_ = 0;
  sbuf_ = 0;
  delete_sbuf_ = 0;
  filename_ = 0;
  basename_ = 0;
  if (!basename.empty()) {
      set_basename(basename.c_str());
    }
  // filename overrides basename
  if (!filename.empty()) {
      set_filename(filename.c_str());
    }
}

void
FileRender::set_filename(const char *filename)
{
  delete[] basename_;
  delete[] filename_;
  filename_ = 0;
  basename_ = 0;
  if (filename) filename_ = strcpy(new char[strlen(filename)+1],filename);
}

void
FileRender::set_basename(const char *basename)
{
  delete[] filename_;
  delete[] basename_;
  filename_ = 0;
  basename_ = 0;
  if (basename) basename_ = strcpy(new char[strlen(basename)+1],basename);
}

FileRender::~FileRender()
{
  delete[] filename_;
  delete[] basename_;
  if (delete_sbuf_) delete sbuf_;
}

void
FileRender::clear()
{
}

char *
FileRender::get_filename(const char *objectname)
{
  char *file = 0;

  if (filename_) {
      // if a file name is given then it is the entire name of the file
      file = strcpy(new char[strlen(filename_) + 1],filename_);
    }
  else if (basename_) {
      // if we have a base name, it is used to construct a filename
      const char *ext = file_extension();
      int lenobjectname;
      if (objectname == 0) lenobjectname = 0;
      else lenobjectname = strlen(objectname);
      if (lenobjectname) lenobjectname++;
      file = new char[strlen(basename_)+lenobjectname+strlen(ext)+1];
      strcpy(file, basename_);
      if (lenobjectname) {
          strcat(file, ".");
          strcat(file, objectname);
        }
      strcat(file, ext);
    }
  else {
      if (!objectname) objectname = "renderedobject";
      const char *ext = file_extension();
      file = new char[strlen(objectname)+strlen(ext)+1];
      strcpy(file, objectname);
      strcat(file, ext);
    }

  return file;
}

void
FileRender::open_sbuf(const char *objectname)
{
  if (depth_++) return;

  char *file = get_filename(objectname);

  if (file) {
      filebuf *fbuf = new filebuf();
      fbuf->open(file,ios::out);
      if (!fbuf->is_open()) {
          ExEnv::errn() << scprintf("FileRender: couldn't open \"%s\"\n", filename_);
          abort();
        }
      sbuf_ = fbuf;
      delete_sbuf_ = 1;
      delete[] file;
    }
}

void
FileRender::close_sbuf()
{
  if (--depth_ > 0) return;

  if (delete_sbuf_) {
      delete sbuf_;
      sbuf_ = 0;
      delete_sbuf_ = 0;
    }
}

const char *
FileRender::file_extension()
{
  return ".ren";
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
