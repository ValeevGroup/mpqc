
#if defined(__GNUC__)
#pragma implementation
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <math/symmetry/tform.h>

// initialize the transformations
static SymSphericalTransform trans2(2);
static SymSphericalTransform trans3(3);
static SymSphericalTransform trans4(4);

static ISymSphericalTransform itrans2(2);
static ISymSphericalTransform itrans3(3);
static ISymSphericalTransform itrans4(4);

inline int
sym_cartindex(int am, int a, int b)
{
  return (((((am+1)<<1)-a)*(a+1))>>1)-b-1;
}

////////////////////////////////////////////////////////////////////////

SymCartesianIter::SymCartesianIter(int l):l_(l)
{
}

SymCartesianIter::~SymCartesianIter()
{
}

int
SymCartesianIter::n()
{
  return ((l_>=0)?((((l_)+2)*((l_)+1))>>1):0);
}

void
SymCartesianIter::start()
{
  bfn_=a_=c_=0;
}

void
SymCartesianIter::next()
{
  if (c_<l_-a_) c_++; else {c_=0; a_++;} bfn_++;
}

SymCartesianIter::operator int()
{
  return a_<=l_;
}

int
SymCartesianIter::a()
{
  return a_;
}

int
SymCartesianIter::b()
{
  return l_-a_-c_;
}

int
SymCartesianIter::c()
{
  return c_;
}

int
SymCartesianIter::l(int i)
{
  return i?(i==1?(l_-a_-c_):c_):a_;
}

int
SymCartesianIter::l()
{
  return l_;
}

int
SymCartesianIter::bfn()
{
  return bfn_;
}

////////////////////////////////////////////////////////////////////////
// RedundantCartianIter

SymRedundantCartesianIter::SymRedundantCartesianIter(int l)
{
  l_ = l;
  axis_ = new int[l_];
}

SymRedundantCartesianIter::~SymRedundantCartesianIter()
{
  delete[] axis_;
}

int
SymRedundantCartesianIter::l(int axis)
{
  int i;
  int r = 0;
  for (i=0; i<l_; i++) if (axis_[i]==axis) r++;
  return r;
}

int
SymRedundantCartesianIter::bfn()
{
  int i = a();
  int j = b();
  int am = l();
  return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
}

SymRedundantCartesianIter::operator int()
{
  return !done_;
}

void
SymRedundantCartesianIter::start()
{
  if (l_==0) done_ = 1;
  else done_ = 0;
  int i;
  for (i=0; i<l_; i++) {
      axis_[i] = 0;
    }
}

void
SymRedundantCartesianIter::next()
{
  int i;
  for (i=0; i<l_; i++) {
      if (axis_[i] == 2) axis_[i] = 0;
      else {
          axis_[i]++;
          return;
        }
    }
  done_ = 1;
}

////////////////////////////////////////////////////////////////////////
// RedundantCartianIter

SymRedundantCartesianSubIter::SymRedundantCartesianSubIter(int l)
{
  l_ = l;
  axis_ = new int[l_];
}

SymRedundantCartesianSubIter::~SymRedundantCartesianSubIter()
{
  delete[] axis_;
}

int
SymRedundantCartesianSubIter::bfn()
{
  int i = a();
  int j = b();
  int am = l();
  return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
}

SymRedundantCartesianSubIter::operator int()
{
  return !done_;
}

void
SymRedundantCartesianSubIter::start(int a, int b, int c)
{
  if (l_ != a + b + c) {
      fprintf(stderr, "RedundantCartesianSubIter::start: bad args\n");
      abort();
    }
  if (l_==0) {
      done_ = 1;
      return;
    }
  else {
      done_ = 0;
    }
  e_[0] = a;
  e_[1] = b;
  e_[2] = c;
  int i;
  for (i=0; i<l_; i++) {
      axis_[i] = 0;
    }
  while (!done_ && !valid()) { advance(); }
}

void
SymRedundantCartesianSubIter::next()
{
  advance();
  while (!done_ && !valid()) { advance(); }
}

void
SymRedundantCartesianSubIter::advance()
{
  int i;
  for (i=0; i<l_; i++) {
      if (axis_[i] == 2) axis_[i] = 0;
      else {
          axis_[i]++;
          return;
        }
    }
  done_ = 1;
}

int
SymRedundantCartesianSubIter::valid()
{
  int t[3];
  int i;

  for (i=0; i<3; i++) t[i] = 0;
  for (i=0; i<l_; i++) t[axis_[i]]++;

  return t[0] == e_[0] && t[1] == e_[1] && t[2] == e_[2];
}

////////////////////////////////////////////////////////////////////////

void
SymSphericalTransform::Component::init(int a, int b, int c, double coef,
                                    int pureindex)
{
  a_ = a;
  b_ = b;
  c_ = c;
  coef_ = coef;
  pureindex_ = pureindex;
  cartindex_ = sym_cartindex(a+b+c,a,b);
}

SymSphericalTransform::SymSphericalTransform()
{
  n_ = 0;
  l_ = 0;
  components_ = 0;
}

SymSphericalTransform::SymSphericalTransform(int l)
{
  n_ = 0;
  l_ = l;
  components_ = 0;
  int i = 0;

  if (l==2) {
    add(0,0,2,  1.0 , i);
    add(2,0,0, -0.5 , i);
    add(0,2,0, -0.5 , i);
    i++;
    add(2,0,0,  0.5 * sqrt(3.0), i);
    add(0,2,0, -0.5 * sqrt(3.0), i);
    i++;
    add(1,1,0,  1.0 * sqrt(3.0), i);
    i++;
    add(0,1,1,  1.0 * sqrt(3.0), i);
    i++;
    add(1,0,1,  1.0 * sqrt(3.0), i);
  } else if (l==3) {
#if 0
      // orthonormal functions
    add(0,0,3,  1.0 , i);
    add(2,0,1, -1.5 , i);
    add(0,2,1, -1.5 , i);
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
#else
      // unnormalized and nonorthogonal
      // the pure to cartesian transform matrix is orthogonal, however
    add(0,0,3,  2.0, i);
    add(2,0,1, -3.0, i);
    add(0,2,1, -3.0, i);
    i++;           
    add(1,0,2,  4.0, i);
    add(3,0,0, -1.2, i);
    add(1,2,0, -0.4, i);
    i++;           
    add(0,1,2,  4.0, i);
    add(0,3,0, -1.2, i);
    add(2,1,0, -0.4, i);
    i++;           
    add(2,0,1,  1.0, i);
    add(0,2,1, -1.0, i);
    i++;           
    add(1,1,1,  1.0, i);
    i++;           
    add(3,0,0,  1.0, i);
    add(1,2,0, -3.0, i);
    i++;           
    add(0,3,0,  1.0, i);
    add(2,1,0, -3.0, i);
#endif
  } else if (l==4) {
#if 0
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
#else
      // unnormalized and nonorthogonal
      // the pure to cartesian transform matrix is orthogonal, however
      add(0,0,4,  8.0, i);
      add(4,0,0,  3.0 + 15.0/19.0, i);
      add(0,4,0,  3.0 + 15.0/19.0, i);
      add(2,0,2,-24.0, i);
      add(0,2,2,-24.0, i);
      add(2,2,0,  6.0 - 6.0 * 15.0/19.0, i);
      i++;
      add(1,0,3,  4.0, i);
      add(3,0,1, -3.6, i);
      add(1,2,1, -1.2, i);
      i++;
      add(0,1,3,  4.0, i);
      add(0,3,1, -3.6, i);
      add(2,1,1, -1.2, i);
      i++;
      add(2,0,2,  6.0, i);
      add(0,2,2, -6.0, i);
      add(4,0,0, -1.0, i);
      add(0,4,0,  1.0, i);
      i++;
      add(1,1,2,  6.0, i);
      add(3,1,0, -1.0, i);
      add(1,3,0, -1.0, i);
      i++;
      add(1,2,1,  3.0, i);
      add(3,0,1, -1.0, i);
      i++;
      add(2,1,1,  3.0, i);
      add(0,3,1, -1.0, i);
      i++;
      add(2,2,0,  6.0, i);
      add(4,0,0, -1.0, i);
      add(0,4,0, -1.0, i);
      i++;
      add(3,1,0,  1.0, i);
      add(1,3,0, -1.0, i);
#endif
    }
  else {
      fprintf(stderr, "SphericalTransform: cannot handle l = %d\n", l);
      abort();
    }
}

SymSphericalTransform::~SymSphericalTransform()
{
  delete[] components_;
}

void
SymSphericalTransform::add(int a, int b, int c, double coef, int pureindex)
{
  Component *ncomp = new Component[n_+1];
  int i;
  for (i=0; i<n_; i++) ncomp[i] = components_[i];
  ncomp[i].init(a, b, c, coef, pureindex);
  delete[] components_;
  components_ = ncomp;
  n_++;
}

SymSphericalTransformIter::SymSphericalTransformIter(SymSphericalTransform*t)
{
  transform_ = t;
}

SymSphericalTransformIter::SymSphericalTransformIter(int l, int inverse)
{
  if (l==2) {
      if (inverse) transform_ = &itrans2;
      else transform_ = &trans2;
    }
  else if (l==3) {
      if (inverse) transform_ = &itrans3;
      else transform_ = &trans3;
    }
  else if (l==4) {
      if (inverse) transform_ = &itrans4;
      else transform_ = &trans4;
    }
  else {
      fprintf(stderr, "SphericalTransformIter: cannot handle l = %d\n", l);
      abort();
    }
}

ISymSphericalTransform::ISymSphericalTransform(int l)
{
  l_ = l;
  
  // this transform is orthogonal, but not normalized
  SymSphericalTransform t(l);

  double *norm = new double[2*l+1];
  int i;
  for (i=0; i<2*l+1; i++) norm[i] = 0.0;

  SymSphericalTransformIter I(&t);

  for (I.begin(); I.ready(); I.next()) {
      norm[I.pureindex()] += I.coef()*I.coef();
    }

  for (i=0; i<2*l+1; i++) norm[i] = 1.0/norm[i];

  for (I.begin(); I.ready(); I.next()) {
      add(I.a(), I.b(), I.c(), I.coef()*norm[I.pureindex()], I.pureindex());
    }

  delete[] norm;
}

////////////////////////////////////////////////////////////////////////

void
SymRotation::done() {
  if (r) {
    for (int i=0; i < _n; i++) {
      if (r[i]) delete[] r[i];
    }
    delete[] r;
    r=0;
  }
  _n=0;
}

SymRotation::SymRotation(int a, SymmetryOperation& so, int pure):
  _am(0),
  _n(0),
  r(0)
{
  if (pure)
    init_pure(a,so);
  else
    init(a,so);
}

SymRotation::~SymRotation()
{
  done();
}

// Compute the transformation matrices for general cartesian shells
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
// Written by clj...blame him
void
SymRotation::init(int a, SymmetryOperation&so)
{
  done();

  _am=a;
  
  SymCartesianIter I(_am);
  SymRedundantCartesianIter J(_am);
  int lI[3];
  int k, iI;
  
  _n = I.n();
  r = new double*[_n];

  for (I.start(); I; I.next()) {
    r[I.bfn()] = new double[_n];
    memset(r[I.bfn()],0,sizeof(double)*_n);

    for (J.start(); J; J.next()) {
      double tmp = 1.0;

      for (k=0; k < 3; k++) {
        lI[k] = I.l(k);
      }
      
      for (k=0; k < _am; k++) {
        for (iI=0; lI[iI]==0; iI++);
        lI[iI]--;
        //tmp *= so(J.axis(k),iI);
        tmp *= so(iI,J.axis(k));
      }

      r[I.bfn()][J.bfn()] += tmp;
    }
  }
}

// Compute the transformation matrices for general pure am
// by summing contributions from the cartesian components
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
void
SymRotation::init_pure(int a, SymmetryOperation&so)
{
  done();

  _am=a;
  
  SymSphericalTransformIter I(_am);
  SymSphericalTransformIter J(_am, 1);
  SymRedundantCartesianSubIter K(_am);
  int lI[3];
  int m, iI;
  
  _n = I.n();

  r = new double*[_n];
  for (m=0; m<_n; m++) {
      r[m] = new double[_n];
      memset(r[m],0,sizeof(double)*_n);
    }

  for (I.start(); I; I.next()) {
    for (J.start(); J; J.next()) {
      double coef = I.coef()*J.coef();
      double tmp = 0.0;

      for (K.start(J.a(), J.b(), J.c()); K; K.next()) {
        double tmp2 = coef;
        for (m=0; m < 3; m++) {
          lI[m] = I.l(m);
        }
        
        for (m=0; m < _am; m++) {
          for (iI=0; lI[iI]==0; iI++);
          lI[iI]--;
          //tmp2 *= so(K.axis(m),iI);
          tmp2 *= so(iI,K.axis(m));
        }
        tmp += tmp2;
      }
      r[I.bfn()][J.bfn()] += tmp;
    }
  }
}

double
SymRotation::trace() const {
  double t=0;
  for (int i=0; i < _n; i++)
    t += r[i][i];
  return t;
}

void
SymRotation::print() const
{
  for (int i=0; i < _n; i++) {
    printf("%5d ",i+1);
    for (int j=0; j < _n; j++) {
      printf(" %10.7f",r[i][j]);
    }
    printf("\n");
  }
}
