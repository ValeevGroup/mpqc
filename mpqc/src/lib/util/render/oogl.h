
#ifndef _util_render_oogl_h
#define _util_render_oogl_h

#include <stdio.h>
#include <util/render/render.h>

class OOGLRender: public Render {
#   define CLASSNAME OOGLRender
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    char* filename_;
    FILE* fp_;
    int oogl_spheres_;
  public:
    OOGLRender(const char * filename);
    OOGLRender(FILE * fp = stdout);
    OOGLRender(const RefKeyVal&);
    virtual ~OOGLRender();

    void render(const RefRenderedObject&);
    void clear();

    void set(const RefRenderedObjectSet&);
    void sphere(const RefRenderedSphere&);
    void polygons(const RefRenderedPolygons&);
};
DescribedClass_REF_dec(OOGLRender);

#endif
