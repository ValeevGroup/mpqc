
#ifndef _util_render_material_h
#define _util_render_material_h

#include <stdio.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/parameter.h>
#include <util/render/color.h>

class Material: public DescribedClass {
#   define CLASSNAME Material
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    Parameter<Color> diffuse_;
    Parameter<Color> ambient_;
  public:
    Material();
    Material(const RefKeyVal&);
    ~Material();
    Parameter<Color>& diffuse() { return diffuse_; }
    Parameter<Color>& ambient() { return ambient_; }
    void print(FILE*fp = stdout);
};
DescribedClass_REF_dec(Material);

#endif
