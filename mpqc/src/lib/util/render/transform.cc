
#include <stdio.h>
#include <stdlib.h>
#include "transform.h"

#define CLASSNAME Transform
#define HAVE_KEYVAL_CTOR
#define PARENTS public DescribedClass
#include <util/class/classi.h>
void *
Transform::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

Transform::Transform(const RefKeyVal& keyval)
{
  transform_ = identity3D();
  if (keyval->exists("translate")) {
      if (keyval->count("translate") != 3) {
          fprintf(stderr,"Transform: error in translation\n");
          abort();
        }
      vec3 r(keyval->doublevalue("translate",0),
             keyval->doublevalue("translate",1),
             keyval->doublevalue("translate",2));
      transform_ = transform_ * translation3D(r);
    }
  if (keyval->exists("rotate")) {
      if (keyval->count("rotate:axis") != 3
          || !keyval->exists("rotate:angle")) {
          fprintf(stderr,"Transform: error in rotation\n");
          abort();
        }
      vec3 axis(keyval->doublevalue("rotate:axis",0),
                keyval->doublevalue("rotate:axis",1),
                keyval->doublevalue("rotate:axis",2));
      transform_ = transform_
                 * rotation3D(axis, keyval->doublevalue("rotate:angle"));
    }
  if (keyval->exists("scale")) {
      double scalefactor = keyval->doublevalue("scale");
      transform_ = transform_
                 * scaling3D(vec3(scalefactor,scalefactor,scalefactor));
    }
}

Transform::~Transform()
{
}

void
Transform::print(FILE*fp)
{
  fprintf(fp, "Transform\n");
}
