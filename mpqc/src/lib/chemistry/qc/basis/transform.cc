
#if defined(__GNUC__)
#pragma implementation
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <chemistry/qc/basis/transform.h>

SphericalTransform::SphericalTransform()
{
  n_ = 0;
  l_ = 0;
  components_ = 0;
}

SphericalTransform::SphericalTransform(int l) : l_(l)
{
  n_ = 0;
  components_ = 0;
}

void
SphericalTransform::init()
{
  int i = 0;

  if (l_==2) {
    add(0,0,2,  2.0 * sqrt(0.25), i);
    add(2,0,0, -1.0 * sqrt(0.25), i);
    add(0,2,0, -1.0 * sqrt(0.25), i);
    i++;
    add(2,0,0,  1.0 * sqrt(0.75), i);
    add(0,2,0, -1.0 * sqrt(0.75), i);
    i++;
    add(1,1,0,  1.0 * sqrt(3.0), i);
    i++;
    add(0,1,1,  1.0 * sqrt(3.0), i);
    i++;
    add(1,0,1,  1.0 * sqrt(3.0), i);

  } else if (l_==3) {
    // orthonormal functions
    add(0,0,3,  2.0 * sqrt(0.25), i);
    add(2,0,1, -3.0 * sqrt(0.25), i);
    add(0,2,1, -3.0 * sqrt(0.25), i);
    i++;
    add(1,0,2,  4.0 * sqrt(0.375), i);
    add(3,0,0, -1.0 * sqrt(0.375), i);
    add(1,2,0, -1.0 * sqrt(0.375), i);
    i++;
    add(0,1,2,  4.0 * sqrt(0.375), i);
    add(0,3,0, -1.0 * sqrt(0.375), i);
    add(2,1,0, -1.0 * sqrt(0.375), i);
    i++;
    add(2,0,1,  1.0 * sqrt(3.75), i);
    add(0,2,1, -1.0 * sqrt(3.75), i);
    i++;
    add(1,1,1,  1.0 * sqrt(15.0), i);
    i++;
    add(3,0,0,  1.0 * sqrt(0.625), i);
    add(1,2,0, -3.0 * sqrt(0.625), i);
    i++;
    add(0,3,0,  1.0 * sqrt(0.625), i);
    add(2,1,0, -3.0 * sqrt(0.625), i);

  } else if (l_==4) {
    // orthonormal functions
    add(0,0,4,  8.0 * sqrt(1.0/64.0), i);
    add(4,0,0,  3.0 * sqrt(1.0/64.0), i);
    add(0,4,0,  3.0 * sqrt(1.0/64.0), i);
    add(2,0,2,-24.0 * sqrt(1.0/64.0), i);
    add(0,2,2,-24.0 * sqrt(1.0/64.0), i);
    add(2,2,0,  6.0 * sqrt(1.0/64.0), i);
    i++;
    add(1,0,3,  4.0 * sqrt(0.625), i);
    add(3,0,1, -3.0 * sqrt(0.625), i);
    add(1,2,1, -3.0 * sqrt(0.625), i);
    i++;
    add(0,1,3,  4.0 * sqrt(0.625), i);
    add(0,3,1, -3.0 * sqrt(0.625), i);
    add(2,1,1, -3.0 * sqrt(0.625), i);
    i++;
    add(2,0,2,  6.0 * sqrt(0.3125), i);
    add(0,2,2, -6.0 * sqrt(0.3125), i);
    add(4,0,0, -1.0 * sqrt(0.3125), i);
    add(0,4,0,  1.0 * sqrt(0.3125), i);
    i++;
    add(1,1,2,  6.0 * sqrt(1.25), i);
    add(3,1,0, -1.0 * sqrt(1.25), i);
    add(1,3,0, -1.0 * sqrt(1.25), i);
    i++;
    add(1,2,1,  3.0 * sqrt(4.375), i);
    add(3,0,1, -1.0 * sqrt(4.375), i);
    i++;
    add(2,1,1,  3.0 * sqrt(4.375), i);
    add(0,3,1, -1.0 * sqrt(4.375), i);
    i++;
    add(2,2,0,  6.0 * sqrt(35.0/64.0), i);
    add(4,0,0, -1.0 * sqrt(35.0/64.0), i);
    add(0,4,0, -1.0 * sqrt(35.0/64.0), i);
    i++;
    add(3,1,0,  1.0 * sqrt(8.75), i);
    add(1,3,0, -1.0 * sqrt(8.75), i);

  } else {
    fprintf(stderr, "SphericalTransform: cannot handle l = %d\n", l_);
    abort();
  }
}

SphericalTransform::~SphericalTransform()
{
  if (components_) {
    delete[] components_;
    components_ = 0;
  }
}

void
SphericalTransform::add(int a, int b, int c, double coef, int pureindex)
{
  int i;

  SphericalTransformComponent *ncomp = new_components();

  for (i=0; i<n_; i++)
    ncomp[i] = components_[i];

  ncomp[i].init(a, b, c, coef, pureindex);

  delete[] components_;
  components_ = ncomp;
  n_++;
}

///////////////////////////////////////////////////////////////////////////

ISphericalTransform::ISphericalTransform() :
  SphericalTransform()
{
}

ISphericalTransform::ISphericalTransform(int l) :
  SphericalTransform(l)
{
}

void
ISphericalTransform::init()
{
  int i = 0;

  if (l_==2) {
    add(0,0,2,  2.0/3.0 , i);
    add(2,0,0, -1.0/3.0, i);
    add(0,2,0, -1.0/3.0, i);
    i++;
    add(2,0,0,  sqrt(1.0/3.0), i);
    add(0,2,0, -sqrt(1.0/3.0), i);
    i++;
    add(1,1,0,  sqrt(1.0/3.0), i);
    i++;
    add(0,1,1,  sqrt(1.0/3.0), i);
    i++;
    add(1,0,1,  sqrt(1.0/3.0), i);

  } else if (l_==3) {
    // orthonormal functions
    add(0,0,3,  2.0/11.0, i);
    add(2,0,1, -3.0/11.0, i);
    add(0,2,1, -3.0/11.0, i);
    i++;
    add(1,0,2,  4.0 * sqrt(1.0/150.0), i);
    add(3,0,0, -3.0 * sqrt(1.0/150.0), i);
    add(1,2,0, -1.0 * sqrt(1.0/150.0), i);
    i++;
    add(0,1,2,  4.0 * sqrt(1.0/150.0), i);
    add(0,3,0, -3.0 * sqrt(1.0/150.0), i);
    add(2,1,0, -1.0 * sqrt(1.0/150.0), i);
    i++;
    add(2,0,1,  sqrt(1.0/15.0), i);
    add(0,2,1, -sqrt(1.0/15.0), i);
    i++;
    add(1,1,1,  sqrt(1.0/15.0), i);
    i++;
    add(3,0,0,  sqrt(1.0/10.0), i);
    add(1,2,0, -sqrt(1.0/10.0), i);
    i++;
    add(0,3,0,  sqrt(1.0/10.0), i);
    add(2,1,0, -sqrt(1.0/10.0), i);

  } else if (l_==4) {
    // orthonormal functions
    add(0,0,4, 19.0/370.0, i);
    add(4,0,0,  9.0/370.0, i);
    add(0,4,0,  9.0/370.0, i);
    add(2,0,2,-57.0/370.0, i);
    add(0,2,2,-57.0/370.0, i);
    add(2,2,0,  3.0/370.0, i);
    i++;
    add(1,0,3,  1.0/19.0 * sqrt(10.0), i);
    add(3,0,1, -9.0/19.0 * sqrt(0.1), i);
    add(1,2,1, -3.0/19.0 * sqrt(0.1), i);
    i++;
    add(0,1,3,  1.0/19.0 * sqrt(10.0), i);
    add(0,3,1, -9.0/19.0 * sqrt(0.1), i);
    add(2,1,1, -3.0/19.0 * sqrt(0.1), i);
    i++;
    add(2,0,2, 12.0/37.0 * sqrt(0.2), i);
    add(0,2,2,-12.0/37.0 * sqrt(0.2), i);
    add(4,0,0, -2.0/37.0 * sqrt(0.2), i);
    add(0,4,0,  2.0/37.0 * sqrt(0.2), i);
    i++;
    add(1,1,2,  6.0/19.0 * sqrt(0.2), i);
    add(3,1,0, -1.0/19.0 * sqrt(0.2), i);
    add(1,3,0, -1.0/19.0 * sqrt(0.2), i);
    i++;
    add(1,0,3,  3.0/19.0 * sqrt(2.0/35.0), i);
    add(3,0,1,-13.0/19.0 * sqrt(1.0/70.0), i);
    add(1,2,1,  3.0/19.0 * sqrt(.7), i);
    i++;
    add(0,1,3,  3.0/19.0 * sqrt(2.0/35.0), i);
    add(0,3,1,-13.0/19.0 * sqrt(1.0/70.0), i);
    add(2,1,1,  3.0/19.0 * sqrt(0.7), i);
    i++;
    add(2,2,0, 93.0/(74.0*sqrt(35.0)), i);
    add(2,0,2,  9.0/(74.0*sqrt(35.0)), i);
    add(0,2,2,  9.0/(74.0*sqrt(35.0)), i);
    add(4,0,0,-17.0/(74.0*sqrt(35.0)), i);
    add(0,4,0,-17.0/(74.0*sqrt(35.0)), i);
    add(0,0,4, -3.0/(74.0*sqrt(35.0)), i);
    i++;
    add(3,1,0,  sqrt(1.0/35.0), i);
    add(1,3,0, -sqrt(1.0/35.0), i);

  } else {
    fprintf(stderr, "ISphericalTransform: cannot handle l = %d\n", l_);
    abort();
  }
}

///////////////////////////////////////////////////////////////////////////

SphericalTransformIter::SphericalTransformIter()
{
  transform_=0;
}

SphericalTransformIter::SphericalTransformIter(SphericalTransform*t)
{
  transform_ = t;
}
