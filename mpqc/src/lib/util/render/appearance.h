
#ifndef _util_render_appearance_h
#define _util_render_appearance_h

#include <stdio.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/parameter.h>

class Appearance: public DescribedClass {
#   define CLASSNAME Appearance
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    Parameter<int> level_; // level of accuracy used to generate spheres, etc
  public:
    Appearance();
    Appearance(const RefKeyVal&);
    ~Appearance();
    Parameter<int>& level() { return level_; }

    void print(FILE*fp);
};
DescribedClass_REF_dec(Appearance);

#endif
