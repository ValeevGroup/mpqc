
#include <util/render/material.h>

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
      Color c(new PrefixKeyVal("diffuse", keyval));
      diffuse_.set(c);
    }
  if (keyval->exists("ambient")) {
      Color c(new PrefixKeyVal("ambient", keyval));
      ambient_.set(c);
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
