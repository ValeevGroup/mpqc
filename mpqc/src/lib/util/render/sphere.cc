
#include <stdio.h>
#include <stdlib.h>

#include "render.h"
#include "object.h"
#include "sphere.h"

#define CLASSNAME RenderedSphere
#define HAVE_KEYVAL_CTOR
#define PARENTS public RenderedObject
#include <util/class/classi.h>
void *
RenderedSphere::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedObject::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedSphere::RenderedSphere()
{
}

RenderedSphere::RenderedSphere(const RefMaterial& material):
  RenderedObject(material)
{
}

RenderedSphere::RenderedSphere(const RefKeyVal& keyval):
  RenderedObject(keyval)
{
}

RenderedSphere::~RenderedSphere()
{
}

void
RenderedSphere::render(const RefRender& render)
{
  render->sphere(this);
}
