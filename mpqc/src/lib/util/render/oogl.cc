//
// oogl.cc
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
#include <util/keyval/keyval.h>
#include <util/render/oogl.h>
#include <util/render/object.h>
#include <util/render/sphere.h>
#include <util/render/polygons.h>
#include <util/render/polylines.h>
#include <util/render/material.h>

#define CLASSNAME OOGLRender
#define PARENTS public Render
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
OOGLRender::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Render::_castdown(cd);
  return do_castdowns(casts,cd);
}

OOGLRender::OOGLRender(const char * filename)
{
  filename_ = strcpy(new char[strlen(filename)+1], filename);
  oogl_spheres_ = 0;
  sbuf_ = 0;
  delete_sbuf_ = 0;
  clear();
}

OOGLRender::OOGLRender(ostream &o)
{
  oogl_spheres_ = 0;
  filename_ = 0;
  sbuf_ = o.rdbuf();
  delete_sbuf_ = 0;
  clear();
}

OOGLRender::OOGLRender(const RefKeyVal& keyval):
  Render(keyval)
{
  filename_ = keyval->pcharvalue("filename");
  oogl_spheres_ = keyval->booleanvalue("oogl_spheres");
  if (filename_) {
      sbuf_ = 0;
      delete_sbuf_ = 0;
    }
  else {
      sbuf_ = cout.rdbuf();
      delete_sbuf_ = 0;
    }
  clear();
}

OOGLRender::~OOGLRender()
{
  delete[] filename_;
  if (delete_sbuf_) delete sbuf_;
}

void
OOGLRender::clear()
{
  if (filename_) {
      if (delete_sbuf_) {
          delete sbuf_;
        }
      filebuf *fbuf = new filebuf();
      fbuf->open(filename_,ios::out);
      if (!fbuf->is_open()) {
          cerr << scprintf("OOGLRender: couldn't open \"%s\"\n", filename_);
          abort();
        }
      sbuf_ = fbuf;
      delete_sbuf_ = 1;
    }
}

void
OOGLRender::render(const RefRenderedObject& object)
{
  ostream o(sbuf_);
  o << "{" << endl;
  if (object->name()) {
      o << "define " << object->name() << endl;
    }
  if (object->transform().nonnull()) {
      o << "= INST" << endl;
      o << "transform {" << endl;
      for (int i=0; i<4; i++) {
          for (int j=0; j<4; j++) {
              o << scprintf(" %10.4f", object->transform()->transform()[j][i]);
            }
          o << endl;
        }
      o << "}" << endl;
      o << "geom {" << endl;
    }
  if (object->material().nonnull()
      ||object->appearance().nonnull()) {
      o << "appearance {" << endl;
      if (object->material().nonnull()) {
          o << "material {" << endl;
          if (object->material()->ambient().is_set()) {
              if (object->material()->ambient().overrides()) o << "*";
              o << scprintf("ambient %10.4f %10.4f %10.4f",
                            object->material()->ambient().value().red(),
                            object->material()->ambient().value().green(),
                            object->material()->ambient().value().blue())
                << endl;
            }
          if (object->material()->diffuse().is_set()) {
              if (object->material()->diffuse().overrides()) o << "*";
              o << scprintf("diffuse %10.4f %10.4f %10.4f",
                            object->material()->diffuse().value().red(),
                            object->material()->diffuse().value().green(),
                            object->material()->diffuse().value().blue())
                << endl;
            }
          o << "}" << endl;
        }
      o << "}" << endl;
    }

  Render::render(object);

  if (object->transform().nonnull()) {
      o << "}" << endl;
    }
  o << "}" << endl;
}

void
OOGLRender::set(const RefRenderedObjectSet& set)
{
  ostream o(sbuf_);
  o << "LIST" << endl;
  for (int i=0; i<set->n(); i++) {
      render(set->element(i));
    }
}

void
OOGLRender::sphere(const RefRenderedSphere& sphere)
{
  if (oogl_spheres_) {
      ostream o(sbuf_);
      o << " = SPHERE 1.0 0.0 0.0 0.0" << endl;
    }
  else {
      Render::sphere(sphere);
    }
}

void
OOGLRender::polygons(const RefRenderedPolygons& poly)
{
  ostream o(sbuf_);
  if (poly->have_vertex_rgb()) {
      o << " = COFF" << endl;
    }
  else {
      o << " = OFF" << endl;
    }
  o << poly->nvertex() << " "
    << poly->nface() << " 0" << endl;
  int i;
  for (i=0; i<poly->nvertex(); i++) {
      o << scprintf(" %10.4f %10.4f %10.4f",
                    poly->vertex(i,0),
                    poly->vertex(i,1),
                    poly->vertex(i,2));
      if (poly->have_vertex_rgb()) {
          // The 1.0 is alpha
          o << scprintf(" %10.4f %10.4f %10.4f 1.0",
                        poly->vertex_rgb(i,0),
                        poly->vertex_rgb(i,1),
                        poly->vertex_rgb(i,2));
        }
      o << endl;
    }
  for (i=0; i<poly->nface(); i++) {
      o << " " << poly->nvertex_in_face(i);
      for (int j=0; j<poly->nvertex_in_face(i); j++) {
          o << " " << poly->face(i,j);
        }
      o << endl;
    }
}

void
OOGLRender::polylines(const RefRenderedPolylines& poly)
{
  int i;
  ostream o(sbuf_);

  int nvertex= 0;
  for (i=0; i<poly->npolyline(); i++) nvertex += poly->nvertex_in_polyline(i);
  o << " = VECT" << endl;
  o << poly->npolyline()
    << nvertex
    << (poly->have_vertex_rgb()? nvertex:0)
    << endl;
  for (i=0; i<poly->npolyline(); i++) {
      o << " " << poly->nvertex_in_polyline(i);
    }
  o << endl;
  if (poly->have_vertex_rgb()) {
      for (i=0; i<poly->npolyline(); i++) {
          o << " " << poly->nvertex_in_polyline(i);
        }
    }
  else {
      for (i=0; i<poly->npolyline(); i++) {
          o << " 0";
        }
    }
  o << endl;
  for (i=0; i<poly->npolyline(); i++) {
      for (int j=0; j<poly->nvertex_in_polyline(i); j++) {
          int ivertex = poly->polyline(i,j);
          o << scprintf(" %10.4f %10.4f %10.4f",
                        poly->vertex(ivertex,0),
                        poly->vertex(ivertex,1),
                        poly->vertex(ivertex,2));
        }
    }
  o << endl;
  if (poly->have_vertex_rgb()) {
      for (i=0; i<poly->npolyline(); i++) {
          for (int j=0; j<poly->nvertex_in_polyline(i); j++) {
              int ivertex = poly->polyline(i,j);
              o << scprintf(" %10.4f %10.4f %10.4f 1.0",
                            poly->vertex_rgb(ivertex,0),
                            poly->vertex_rgb(ivertex,1),
                            poly->vertex_rgb(ivertex,2));
            }
        }
      o << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

