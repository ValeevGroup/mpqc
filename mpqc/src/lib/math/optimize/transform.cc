
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/transform.h>
#include <util/misc/formio.h>

////////////////////////////////////////////////////////////////////////////

REF_def(NonlinearTransform);

NonlinearTransform::~NonlinearTransform()
{
}

void
NonlinearTransform::transform_gradient(const RefSCVector& g)
{
  if (g.null()) return;
  g.assign(linear_transform_ * g);
}

void
NonlinearTransform::transform_hessian(const RefSymmSCMatrix& h)
{
  if (h.null()) return;
  cout << indent
       << "WARNING: NonlinearTransform::transform_hessian: "
       << "using linear transform\n";
  RefSymmSCMatrix newh = h->clone();
  newh.assign(0.0);
  newh->accumulate_transform(linear_transform_.pointer(), h.pointer());
  h.assign(newh);
}

void
NonlinearTransform::transform_ihessian(const RefSymmSCMatrix &ih)
{
  if (ih.null()) return;
  RefSymmSCMatrix h(ih.gi());
  transform_hessian(h);
  ih.assign(h.gi());
}

////////////////////////////////////////////////////////////////////////////

IdentityTransform::~IdentityTransform()
{
}

void
IdentityTransform::transform_coordinates(const RefSCVector& x)
{
}

void
IdentityTransform::transform_gradient(const RefSCVector& g)
{
}

void
IdentityTransform::transform_hessian(const RefSymmSCMatrix& h)
{
}

void
IdentityTransform::transform_ihessian(const RefSymmSCMatrix &ih)
{
}
