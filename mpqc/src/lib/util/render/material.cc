
#include "material.h"

Color::Color(const RefKeyVal& keyval)
{
  red_ = keyval->doublevalue(0);
  green_ = keyval->doublevalue(1);
  blue_ = keyval->doublevalue(2);
}

#define CLASSNAME Material
#define HAVE_KEYVAL_CTOR
#define PARENTS public DescribedClass
#include <util/class/classi.h>
void *
Material::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

Material::Material()
{
  diffuse_.set(Color(0.5, 0.5, 0.5));
  ambient_.set(Color(0.5, 0.5, 0.5));
}

Material::Material(const RefKeyVal& keyval)
{
  if (keyval->exists("diffuse")) {
      diffuse_.set(Color(new PrefixKeyVal("diffuse", keyval)));
    }
  if (keyval->exists("ambient")) {
      ambient_.set(Color(new PrefixKeyVal("ambient", keyval)));
    }
}

Material::~Material()
{
}

void
Material::print(FILE*fp)
{
  fprintf(fp, "Material\n");
}
