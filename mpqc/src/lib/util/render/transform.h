
#ifndef _util_render_transform_h
#define _util_render_transform_h

#include <stdio.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/algebra3.h>

class Transform: public DescribedClass {
#   define CLASSNAME Transform
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    mat4 transform_;
  public:
    Transform() { transform_ = identity3D(); }
    Transform(const RefKeyVal&);
    ~Transform();
    mat4& transform() { return transform_; }
    void print(FILE*fp = stdout);
};
DescribedClass_REF_dec(Transform);

#endif
