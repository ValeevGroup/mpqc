
#include <stdio.h>
#include <stdlib.h>
#include <util/render/render.h>
#include <util/render/object.h>

#define CLASSNAME RenderedObject
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
RenderedObject::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedObject::RenderedObject(const RefMaterial& material):
  material_(material),
  name_(0)
{
}

RenderedObject::RenderedObject(const RefKeyVal& keyval):
  material_(keyval->describedclassvalue("material")),
  transform_(keyval->describedclassvalue("transform")),
  appearance_(keyval->describedclassvalue("appearance")),
  name_(keyval->pcharvalue("name"))
{
}

RenderedObject::~RenderedObject()
{
  
  if (name_) delete[] name_;
}

void
RenderedObject::print(FILE*fp)
{
  fprintf(fp, "RenderedObject:\n");
  if (material_.nonnull()) {
      fprintf(fp, "  material = 0x%x\n", material_.pointer());
    }
  if (appearance_.nonnull()) {
      fprintf(fp, "  appearance = 0x%x\n", appearance_.pointer());
    }
  if (transform_.nonnull()) {
      fprintf(fp, "  transform = 0x%x\n", transform_.pointer());
    }
}
  

#define CLASSNAME RenderedObjectSet
#define HAVE_KEYVAL_CTOR
#define PARENTS public RenderedObject
#include <util/class/classi.h>
void *
RenderedObjectSet::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedObject::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedObjectSet::RenderedObjectSet(int capacity)
{
  capacity_ = capacity;
  n_ = 0;
  array_ = new RefRenderedObject[capacity_];
}

RenderedObjectSet::RenderedObjectSet(const RefKeyVal& keyval):
  RenderedObject(keyval)
{
  capacity_ = keyval->count("objects");
  if (keyval->error() != KeyVal::OK) {
      fprintf(stderr,"RenderedObjectSet: error counting objects\n");
      abort();
    }
  n_ = capacity_;
  array_ = new RefRenderedObject[capacity_];
  for (int i=0; i<n_; i++) {
      array_[i] = keyval->describedclassvalue("objects",i);
      if (keyval->error() != KeyVal::OK) {
          fprintf(stderr,"RenderedObjectSet: error reading objects\n");
          abort();
        }
    }
}

RenderedObjectSet::~RenderedObjectSet()
{
  delete[] array_;
}

void
RenderedObjectSet::add(const RefRenderedObject& object)
{
  if (capacity_ == n_) {
      capacity_ += 10;
      RefRenderedObject *tmp = new RefRenderedObject[capacity_];
      for (int i=0; i<n_; i++) {
          tmp[i] = array_[i];
        }
      delete[] array_;
      array_ = tmp;
    }
  array_[n_] = object;
  n_++;
}

void
RenderedObjectSet::render(const RefRender& render)
{
  render->set(this);
}

