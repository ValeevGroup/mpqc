
#include <util/keyval/keyval.h>
#include "oogl.h"
#include "object.h"
#include "sphere.h"
#include "polygons.h"
#include "material.h"

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
  fp_ = 0;
  clear();
}

OOGLRender::OOGLRender(FILE * fp)
{
  filename_ = 0;
  fp_ = fp;
  clear();
}

OOGLRender::OOGLRender(const RefKeyVal& keyval):
  Render(keyval)
{
  filename_ = keyval->pcharvalue("filename");
  oogl_spheres_ = keyval->booleanvalue("oogl_spheres");
  if (!filename_) fp_ = stdout;
  else fp_ = 0;
  clear();
}

OOGLRender::~OOGLRender()
{
  if (filename_) {
      delete[] filename_;
      if (fp_) fclose(fp_);
    }
}

void
OOGLRender::clear()
{
  if (filename_) {
      if (fp_) {
          fclose(fp_);
        }
      fp_ = fopen(filename_, "w");
      if (!fp_) {
          fprintf(stderr,"OOGLRender: couldn't open \"%s\"\n", filename_);
          abort();
        }
    }
}

void
OOGLRender::render(const RefRenderedObject& object)
{
  fprintf(fp_, "{\n");
  if (object->name()) {
      fprintf(fp_, "define %s\n", object->name());
    }
  if (object->transform().nonnull()) {
      fprintf(fp_, "= INST\n");
      fprintf(fp_, "transform {\n");
      for (int i=0; i<4; i++) {
          for (int j=0; j<4; j++) {
              fprintf(fp_, " %10.4f", object->transform()->transform()[j][i]);
            }
          fprintf(fp_,"\n");
        }
      fprintf(fp_,"}\n");
      fprintf(fp_,"geom {\n");
    }
  if (object->material().nonnull()
      ||object->appearance().nonnull()) {
      fprintf(fp_, "appearance {\n");
      if (object->material().nonnull()) {
          fprintf(fp_, "material {\n");
          if (object->material()->ambient().is_set()) {
              if (object->material()->ambient().overrides()) fprintf(fp_, "*");
              fprintf(fp_, "ambient %10.4f %10.4f %10.4f\n",
                      object->material()->ambient().value().red(),
                      object->material()->ambient().value().green(),
                      object->material()->ambient().value().blue());
            }
          if (object->material()->diffuse().is_set()) {
              if (object->material()->diffuse().overrides()) fprintf(fp_, "*");
              fprintf(fp_, "diffuse %10.4f %10.4f %10.4f\n",
                      object->material()->diffuse().value().red(),
                      object->material()->diffuse().value().green(),
                      object->material()->diffuse().value().blue());
            }
          fprintf(fp_, "}\n");
        }
      fprintf(fp_, "}\n");
    }

  Render::render(object);

  if (object->transform().nonnull()) {
      fprintf(fp_,"}\n");
    }
  fprintf(fp_, "}\n");
}

void
OOGLRender::set(const RefRenderedObjectSet& set)
{
  fprintf(fp_,"LIST\n");
  for (int i=0; i<set->n(); i++) {
      render(set->element(i));
    }
}

void
OOGLRender::sphere(const RefRenderedSphere& sphere)
{
  if (oogl_spheres_) {
      fprintf(fp_," = SPHERE 1.0 0.0 0.0 0.0\n");
    }
  else {
      Render::sphere(sphere);
    }
}

void
OOGLRender::polygons(const RefRenderedPolygons& poly)
{
  if (poly->have_vertex_rgb()) {
      fprintf(fp_, " = COFF\n");
    }
  else {
      fprintf(fp_," = OFF\n");
    }
  fprintf(fp_, "%d %d 0\n", poly->nvertex(), poly->nface());
  for (int i=0; i<poly->nvertex(); i++) {
      fprintf(fp_, " %10.4f %10.4f %10.4f",
              poly->vertex(i,0),
              poly->vertex(i,1),
              poly->vertex(i,2));
      if (poly->have_vertex_rgb()) {
          // The 1.0 is alpha
          fprintf(fp_, " %10.4f %10.4f %10.4f 1.0",
                  poly->vertex_rgb(i,0),
                  poly->vertex_rgb(i,1),
                  poly->vertex_rgb(i,2));
        }
      fprintf(fp_, "\n");
    }
  for (i=0; i<poly->nface(); i++) {
      fprintf(fp_, " %d", poly->nvertex_in_face(i));
      for (int j=0; j<poly->nvertex_in_face(i); j++) {
          fprintf(fp_, " %d", poly->face(i,j));
        }
      fprintf(fp_, "\n");
    }
}

