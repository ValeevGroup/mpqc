
#ifndef _util_render_sphere_h
#define _util_render_sphere_h

#include <util/render/object.h>

class RenderedSphere: public RenderedObject {
#   define CLASSNAME RenderedSphere
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    void render(const RefRender&);
  public:
    RenderedSphere(const RefMaterial&);
    RenderedSphere(const RefKeyVal&);
    ~RenderedSphere();
};
DescribedClass_REF_dec(RenderedSphere);

#endif
