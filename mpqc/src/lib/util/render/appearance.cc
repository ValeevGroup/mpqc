
#include <util/render/appearance.h>

#define CLASSNAME Appearance
#define HAVE_KEYVAL_CTOR
#define PARENTS public DescribedClass
#include <util/class/classi.h>
void *
Appearance::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

Appearance::Appearance()
{
  level_.set(1);
}

Appearance::Appearance(const RefKeyVal& keyval)
{
  int level = keyval->intvalue("level");
  if (keyval->error() == KeyVal::OK) level_.set(level);
}

Appearance::~Appearance()
{
}

void
Appearance::print(FILE*fp)
{
  fprintf(fp, "Appearance:\n");
  fprintf(fp, "  level is ");
  if (level_.is_set()) {
      fprintf(fp, "set to %d", level_.value());
      if (level_.overrides()) {
          fprintf(fp, " and overrides");
        }
    }
  else fprintf(fp, "not set");
  fprintf(fp, "\n");
}
