//
// shape.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <util/misc/math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/isosurf/shape.h>

using namespace std;
using namespace sc;

static const double shape_infinity = 1.0e23;

// given a vector X find which of the points in the vector of
// vectors, A, is closest to it and return the distance
static double
closest_distance(SCVector3& X,SCVector3*A,int n,SCVector3*grad)
{
  SCVector3 T = X-A[0];
  double min = T.dot(T);
  int imin = 0;
  for (int i=1; i<n; i++) {
      T = X-A[i];
      double tmp = T.dot(T);
      if (tmp < min) {min = tmp; imin = i;}
    }
  if (grad) {
      T = X - A[imin];
      T.normalize();
      *grad = T;
    }
  return sqrt(min);
}

//////////////////////////////////////////////////////////////////////
// Shape

static ClassDesc Shape_cd(
  typeid(Shape),"Shape",1,"public Volume",
  0, 0, 0);

Shape::Shape():
  Volume()
{
}

Shape::Shape(const Ref<KeyVal>& keyval):
  Volume(keyval)
{
}

Shape::~Shape()
{
}

void
Shape::compute()
{
  SCVector3 r;
  get_x(r);
  if (gradient_needed()) {
      if (!gradient_implemented()) {
          ExEnv::errn() << "Shape::compute: gradient not implemented" << endl;
          abort();
        }
      SCVector3 v;
      set_value(distance_to_surface(r,&v));
      set_actual_value_accuracy(desired_value_accuracy());
      set_gradient(v);
      set_actual_gradient_accuracy(desired_gradient_accuracy());
    }
  else if (value_needed()) {
      set_value(distance_to_surface(r));
      set_actual_value_accuracy(desired_value_accuracy());
    }
  if (hessian_needed()) {
      ExEnv::errn() << "Shape::compute(): can't do hessian yet" << endl;
      abort();
    }
}

int
Shape::is_outside(const SCVector3&r) const
{
  if (distance_to_surface(r)>0.0) return 1;
  return 0;
}

// Shape overrides volume's interpolate so it always gets points on
// the outside of the shape are always returned.

// interpolate using the bisection algorithm
void
Shape::interpolate(const SCVector3& A,
                   const SCVector3& B,
                   double val,
                   SCVector3& result)
{
  if (val < 0.0) {
      failure("Shape::interpolate(): val is < 0.0");
    }

  set_x(A);
  double value0 = value() - val;

  set_x(B);
  double value1 = value() - val;

  if (value0*value1 > 0.0) {
      failure("Shape::interpolate(): values at endpoints don't bracket val");
    }
  else if (value0 == 0.0) {
      result = A;
      return;
    }
  else if (value1 == 0.0) {
      result = B;
      return;
    }

  SCVector3 BA = B - A;

  double length = BA.norm();
  int niter = (int) (log(length/interpolation_accuracy())/M_LN2);

  double f0 = 0.0;
  double f1 = 1.0;
  double fnext = 0.5;

  SCVector3 X = A + fnext*BA;
  set_x(X);
  double valuenext = value() - val;

  do {
      for (int i=0; i<niter; i++) {
          if (valuenext*value0 <= 0.0) {
              value1 = valuenext;
              f1 = fnext;
              fnext = (f0 + fnext)*0.5;
            }
          else {
              value0 = valuenext;
              f0 = fnext;
              fnext = (fnext + f1)*0.5;
            }
          X = A + fnext*BA;
          set_x(X);
          valuenext = value() - val;
        }
      niter = 1;
    } while (valuenext < 0.0);

  result = X;
}

int
Shape::value_implemented() const
{
  return 1;
}

//////////////////////////////////////////////////////////////////////
// SphereShape

static ClassDesc SphereShape_cd(
  typeid(SphereShape),"SphereShape",1,"public Shape",
  0, create<SphereShape>, 0);

SphereShape::SphereShape(const SCVector3&o,double r):
  _origin(o),
  _radius(r)
{
}

SphereShape::SphereShape(const SphereShape&s):
  _origin(s._origin),
  _radius(s._radius)
{
}

SphereShape::SphereShape(const Ref<KeyVal>& keyval):
  _origin(new PrefixKeyVal(keyval,"origin")),
  _radius(keyval->doublevalue("radius"))
{
}

SphereShape::~SphereShape()
{
}

double
SphereShape::distance_to_surface(const SCVector3&p,SCVector3*grad) const
{
  int i;
  double r2 = 0.0;
  for (i=0; i<3; i++) {
      double tmp = p[i] - _origin[i];
      r2 += tmp*tmp;
    }
  double r = sqrt(r2);
  double d = r - _radius;
  if (grad) {
      SCVector3 pv(p);
      SCVector3 o(_origin);
      SCVector3 unit = pv - o;
      unit.normalize();
      for (i=0; i<3; i++) grad->elem(i) = unit[i];
    }
  return d;
}

void SphereShape::print(ostream&o) const
{
  o << indent
    << scprintf("SphereShape: r = %8.4f o = (%8.4f %8.4f %8.4f)",
                radius(),origin()[0],origin()[1],origin()[2])
    << endl;
}

void
SphereShape::boundingbox(double valuemin, double valuemax,
                         SCVector3& p1,
                         SCVector3& p2)
{
  if (valuemax < 0.0) valuemax = 0.0;

  int i;
  for (i=0; i<3; i++) {
      p1[i] = _origin[i] - _radius - valuemax;
      p2[i] = _origin[i] + _radius + valuemax;
    }
}

int
SphereShape::gradient_implemented() const
{
  return 1;
}

////////////////////////////////////////////////////////////////////////
// UncappedTorusHoleShape

static ClassDesc UncappedTorusHoleShape_cd(
  typeid(UncappedTorusHoleShape),"UncappedTorusHoleShape",1,"public Shape",
  0, 0, 0);

UncappedTorusHoleShape::UncappedTorusHoleShape(double r,
                               const SphereShape& s1,
                               const SphereShape& s2):
_s1(s1),
_s2(s2),
_r(r)
{
}

UncappedTorusHoleShape*
UncappedTorusHoleShape::newUncappedTorusHoleShape(double r,
                                                  const SphereShape&s1,
                                                  const SphereShape&s2)
{
  // if the probe sphere fits between the two spheres, then there
  // is no need to include this shape
  SCVector3 A(s1.origin());
  SCVector3 B(s2.origin());
  SCVector3 BA = B - A;
  if (2.0*r <  BA.norm() - s1.radius() - s2.radius()) return 0;

  // check to see if the surface is reentrant
  double rrs1 = r+s1.radius();
  double rrs2 = r+s2.radius();
  SCVector3 R12 = ((SCVector3)s1.origin()) - ((SCVector3)s2.origin());
  double r12 = sqrt(R12.dot(R12));
  if (sqrt(rrs1*rrs1-r*r) + sqrt(rrs2*rrs2-r*r) < r12)
    return new ReentrantUncappedTorusHoleShape(r,s1,s2);

  // otherwise create an ordinary torus hole
  return new NonreentrantUncappedTorusHoleShape(r,s1,s2);
}

// Given a node, finds a sphere in the plane of n and the centers
// of _s1 and _s2 that touches the UncappedTorusHole.  There are two
// candidates, the one closest to n is chosen.
void
UncappedTorusHoleShape::in_plane_sphere(
    const SCVector3& n,
    SCVector3& P) const
{
  // the center of the sphere is given by the vector equation
  // P = A + a R(AB) + b U(perp),
  // where
  // A is the center of _s1
  // B is the center of _s2
  // R(AB) is the vector from A to B, R(AB) = B - A
  // U(perp) is a unit vect perp to R(AB) and lies in the plane of n, A, and B
  // The unknown scalars, a and b are given by solving the following
  // equations:
  // | P - A | = r(A) + _r, and
  // | P - B | = r(B) + _r,
  // which give
  // | a R(AB) + b U(perp) | = r(A) + _r, and
  // | (a-1) R(AB) + b U(perp) | = r(B) + _r.
  // These further simplify to
  // a^2 r(AB)^2 + b^2 = (r(A)+_r)^2, and
  // (a-1)^2 r(AB)^2 + b^2 = (r(B)+_r)^2.
  // Thus,
  // a = (((r(A)+_r)^2 - (r(B)+_r)^2 )/(2 r(AB)^2)) + 1/2
  // b^2 = (r(A)+r)^2 - a^2 r(AB)^2

  SCVector3 A = _s1.origin();
  SCVector3 B = _s2.origin();
  SCVector3 N = n;
  SCVector3 R_AB = B - A;
  SCVector3 R_AN = N - A;

  // vector projection of R_AN onto R_AB
  SCVector3 P_AN_AB = R_AB * (R_AN.dot(R_AB)/R_AB.dot(R_AB));
  // the perpendicular vector
  SCVector3 U_perp = R_AN - P_AN_AB;

  // if |U| is tiny, then any vector perp to AB will do
  double u = U_perp.dot(U_perp);
  if (u<1.0e-23) {
      SCVector3 vtry = R_AB;
      vtry[0] += 1.0;
      vtry = vtry - R_AB * (vtry.dot(R_AB)/R_AB.dot(R_AB));
      if (vtry.dot(vtry) < 1.0e-23) {
          vtry = R_AB;
          vtry[1] += 1.0;
          vtry = vtry - R_AB * (vtry.dot(R_AB)/R_AB.dot(R_AB));
        }
      U_perp = vtry;
    }

  U_perp.normalize();
  //ExEnv::outn() << "A: "; A.print();
  //ExEnv::outn() << "U_perp: "; U_perp.print();
  //ExEnv::outn() << "R_AB: "; R_AB.print();

  double r_AB = sqrt(R_AB.dot(R_AB));
  double r_A = _s1.radius();
  double r_B = _s2.radius();

  double r_Ar = r_A + _r;
  double r_Br = r_B + _r;
  double a = ((r_Ar*r_Ar - r_Br*r_Br)/(2.0*r_AB*r_AB)) + 0.5;
  double b = sqrt(r_Ar*r_Ar - a*a*r_AB*r_AB);

  //ExEnv::outn() << scprintf("r_Ar = %f, r_AB = %f\n",r_Ar,r_AB);
  //ExEnv::outn() << scprintf("a = %f, b = %f\n",a,b);

  P = A + a * R_AB + b * U_perp;
  //ExEnv::outn() << "a*R_AB: "; (a*R_AB).print();
  //ExEnv::outn() << "b*U_perp: "; (b*U_perp).print();
}

void
UncappedTorusHoleShape::print(ostream&o) const
{
  o << indent << "UncappedTorusHoleShape:" << endl;
  o << incindent;
  o << indent << "r = " << _r << endl;
  o << indent << "s1 = ";
  o << incindent << skipnextindent;
  _s1.print(o);
  o << decindent;
  o << indent << "s2 = ";
  o << incindent << skipnextindent;
  _s2.print(o);
  o << decindent;
  o << decindent;
}

void
UncappedTorusHoleShape::boundingbox(double valuemin, double valuemax,
                                    SCVector3& p1,
                                    SCVector3& p2)
{
  SCVector3 p11;
  SCVector3 p12;
  SCVector3 p21;
  SCVector3 p22;

  _s1.boundingbox(valuemin,valuemax,p11,p12);
  _s2.boundingbox(valuemin,valuemax,p21,p22);

  int i;
  for (i=0; i<3; i++) {
      if (p11[i] < p21[i]) p1[i] = p11[i];
      else p1[i] = p21[i];
      if (p12[i] > p22[i]) p2[i] = p12[i];
      else p2[i] = p22[i];
    }
}

int
UncappedTorusHoleShape::gradient_implemented() const
{
  return 1;
}

/////////////////////////////////////////////////////////////////////
// is in triangle

static int
is_in_unbounded_triangle(const SCVector3&XP,const SCVector3&AP,const SCVector3&BP)
{
  SCVector3 axis = BP.cross(AP);

  SCVector3 BP_perp = BP; BP_perp.rotate(M_PI_2,axis);
  double u = BP_perp.dot(XP)/BP_perp.dot(AP);
  if (u<0.0) return 0;

  SCVector3 AP_perp = AP; AP_perp.rotate(M_PI_2,axis);
  double w = AP_perp.dot(XP)/AP_perp.dot(BP);
  if (w<0.0) return 0;

  return 1;
}

/////////////////////////////////////////////////////////////////////
// ReentrantUncappedTorusHoleShape

static ClassDesc ReentrantUncappedTorusHoleShape_cd(
  typeid(ReentrantUncappedTorusHoleShape),"ReentrantUncappedTorusHoleShape",1,"public UncappedTorusHoleShape",
  0, 0, 0);

ReentrantUncappedTorusHoleShape::ReentrantUncappedTorusHoleShape(double r,
                                                 const SphereShape& s1,
                                                 const SphereShape& s2):
  UncappedTorusHoleShape(r,s1,s2)
{
  rAP = r + s1.radius();
  rBP = r + s2.radius();
  BA = B() - A();

  // Find the points at the ends of the two cone-like objects, I[0] and I[1].
  // they are given by:
  //   I = z BA, where BA = B-A and I is actually IA = I - A
  //   r^2 = PI.PI, where PI = PA-I and P is the center of a probe sphere
  // this gives
  //  z^2 BA.BA - 2z PA.BA + PA.PA - r^2 = 0

  SCVector3 arbitrary; 
  arbitrary[0] = arbitrary[1] = arbitrary[2] = 0.0;
  SCVector3 P;
  in_plane_sphere(arbitrary,P);
  SCVector3 PA = P - A();

  double a = BA.dot(BA);
  double minus_b = 2.0 * PA.dot(BA);
  double c = PA.dot(PA) - r*r;
  double b2m4ac = minus_b*minus_b - 4*a*c;
  double sb2m4ac;
  if (b2m4ac >= 0.0) {
      sb2m4ac = sqrt(b2m4ac);
    }
  else if (b2m4ac > -1.0e-10) {
      sb2m4ac = 0.0;
    }
  else {
      ExEnv::errn() << "ReentrantUncappedTorusHoleShape:: imaginary point" << endl;
      abort();
    }
  double zA = (minus_b - sb2m4ac)/(2.0*a);
  double zB = (minus_b + sb2m4ac)/(2.0*a);
  I[0] = BA * zA + A();
  I[1] = BA * zB + A();
}
ReentrantUncappedTorusHoleShape::~ReentrantUncappedTorusHoleShape()
{
}
int
ReentrantUncappedTorusHoleShape::
  is_outside(const SCVector3&X) const
{
  SCVector3 Xv(X);

  SCVector3 P;
  in_plane_sphere(X,P);
  SCVector3 XP = Xv - P;
  double rXP = XP.norm();
  if (rXP > rAP || rXP > rBP) return 1;

  SCVector3 AP = A() - P;
  SCVector3 BP = B() - P;

  if (!is_in_unbounded_triangle(XP,AP,BP)) return 1;

  if (rXP < radius()) {
      return 1;
    }

  return 0;
}
double
ReentrantUncappedTorusHoleShape::
  distance_to_surface(const SCVector3&X,SCVector3*grad) const
{
  SCVector3 Xv(X);

  SCVector3 P;
  in_plane_sphere(X,P);
  SCVector3 XP = Xv - P;
  double rXP = XP.norm();
  if (rXP > rAP || rXP > rBP) return shape_infinity;

  SCVector3 AP = A() - P;
  SCVector3 BP = B() - P;

  if (!is_in_unbounded_triangle(XP,AP,BP)) return shape_infinity;

  SCVector3 I1P = I[0] - P;
  SCVector3 I2P = I[1] - P;

  if (!is_in_unbounded_triangle(XP,I1P,I2P)) {
      if (rXP < radius()) {
          if (grad) {
              SCVector3 unit(XP);
              unit.normalize();
              *grad = unit;
            }
          return radius() - rXP;
        }
      else return -1.0;
    }

  return closest_distance(Xv,(SCVector3*)I,2,grad);
}

int
ReentrantUncappedTorusHoleShape::gradient_implemented() const
{
  return 1;
}

/////////////////////////////////////////////////////////////////////
// NonreentrantUncappedTorusHoleShape

static ClassDesc NonreentrantUncappedTorusHoleShape_cd(
  typeid(NonreentrantUncappedTorusHoleShape),"NonreentrantUncappedTorusHoleShape",1,"public UncappedTorusHoleShape",
  0, 0, 0);

NonreentrantUncappedTorusHoleShape::
  NonreentrantUncappedTorusHoleShape(double r,
                                     const SphereShape& s1,
                                     const SphereShape& s2):
  UncappedTorusHoleShape(r,s1,s2)
{
  rAP = r + s1.radius();
  rBP = r + s2.radius();
  BA = B() - A();
}
NonreentrantUncappedTorusHoleShape::~NonreentrantUncappedTorusHoleShape()
{
}
double NonreentrantUncappedTorusHoleShape::
  distance_to_surface(const SCVector3&X,SCVector3* grad) const
{
  SCVector3 Xv(X);

  SCVector3 P;
  in_plane_sphere(X,P);
  SCVector3 PX = P - Xv;
  double rPX = PX.norm();
  if (rPX > rAP || rPX > rBP) return shape_infinity;

  SCVector3 PA = P - A();
  SCVector3 XA = Xv - A();

  SCVector3 axis = BA.cross(PA);

  SCVector3 BA_perp = BA; BA_perp.rotate(M_PI_2,axis);
  double u = BA_perp.dot(XA)/BA_perp.dot(PA);
  if (u<0.0 || u>1.0) return shape_infinity;

  SCVector3 PA_perp = PA; PA_perp.rotate(M_PI_2,axis);
  double w = PA_perp.dot(XA)/PA_perp.dot(BA);
  if (w<0.0 || w>1.0) return shape_infinity;

  double uw = u+w;
  if (uw<0.0 || uw>1.0) return shape_infinity;

  if (rPX < radius()) {
      if (grad) {
          SCVector3 unit(PX);
          unit.normalize();
          *grad = unit;
        }
      return radius() - rPX;
    }

  return -1;
}

int
NonreentrantUncappedTorusHoleShape::gradient_implemented() const
{
  return 1;
}

/////////////////////////////////////////////////////////////////////
// Uncapped5SphereExclusionShape

static ClassDesc Uncapped5SphereExclusionShape_cd(
  typeid(Uncapped5SphereExclusionShape),"Uncapped5SphereExclusionShape",1,"public Shape",
  0, 0, 0);

Uncapped5SphereExclusionShape*
Uncapped5SphereExclusionShape::
  newUncapped5SphereExclusionShape(double r,
                                   const SphereShape& s1,
                                   const SphereShape& s2,
                                   const SphereShape& s3)
{
  Uncapped5SphereExclusionShape* s =
    new Uncapped5SphereExclusionShape(r,s1,s2,s3);
  if (s->solution_exists()) {
      return s;
    }
  delete s;
  return 0;
}
static int verbose = 0;
Uncapped5SphereExclusionShape::
  Uncapped5SphereExclusionShape(double radius,
                                const SphereShape&s1,
                                const SphereShape&s2,
                                const SphereShape&s3):
  _s1(s1),
  _s2(s2),
  _s3(s3),
  _r(radius)
{
  double rAr = rA() + r();
  double rAr2 = rAr*rAr;
  double rBr = rB() + r();
  double rBr2 = rBr*rBr;
  double rCr = rC() + r();
  double rCr2 = rCr*rCr;
  double A2 = A().dot(A());
  double B2 = B().dot(B());
  double C2 = C().dot(C());
  SCVector3 BA = B()-A();
  double DdotBA = 0.5*(rAr2 - rBr2 + B2 - A2);
  double DAdotBA = DdotBA - A().dot(BA);
  double BA2 = BA.dot(BA);
  SCVector3 CA = C() - A();
  double CAdotBA = CA.dot(BA);
  SCVector3 CA_perpBA = CA - (CAdotBA/BA2)*BA;
  double CA_perpBA2 = CA_perpBA.dot(CA_perpBA);
  if (CA_perpBA2 < 1.0e-23) {
      _solution_exists = 0;
      return;
    }
  double DdotCA_perpBA = 0.5*(rAr2 - rCr2 + C2 - A2)
    - CAdotBA*DdotBA/BA2;
  double DAdotCA_perpBA = DdotCA_perpBA - A().dot(CA_perpBA);
  double rAt2 = DAdotBA*DAdotBA/BA2 + DAdotCA_perpBA*DAdotCA_perpBA/CA_perpBA2;
  double h2 = rAr2 - rAt2;
  if (h2 <= 0.0) {
      _solution_exists = 0;
      return;
    }
  _solution_exists = 1;

  double h = sqrt(h2);
  if (h<r()) {
      _reentrant = 1;
      //ExEnv::outn() << "WARNING: throwing out reentrant shape" << endl;
      //_solution_exists = 0;
      //return;
    }
  else {
      _reentrant = 0;
      //ExEnv::outn() << "WARNING: throwing out nonreentrant shape" << endl;
      //_solution_exists = 0;
      //return;
    }

  // The projection of D into the ABC plane
  SCVector3 MA = (DAdotBA/BA2)*BA + (DAdotCA_perpBA/CA_perpBA2)*CA_perpBA;
  M = MA + A();
  SCVector3 BAxCA = BA.cross(CA);
  D[0] = M + h * BAxCA * ( 1.0/BAxCA.norm() );
  D[1] = M - h * BAxCA * ( 1.0/BAxCA.norm() );

  // The projection of D into the ABC plane, M, does not lie in the
  // ABC triangle, then this shape must be treated carefully and it
  // will be designated as folded.
  SCVector3 MC = M - C();
  if (!(is_in_unbounded_triangle(MA, BA, CA)
        &&is_in_unbounded_triangle(MC, B() - C(), A() - C()))) {
      _folded = 1;
      SCVector3 MB = M - B();
      double MA2 = MA.dot(MA);
      double MB2 = MB.dot(MB);
      double MC2 = MC.dot(MC);
      if (MA2 < MB2) {
          F1 = A();
          if (MB2 < MC2) F2 = B();
          else F2 = C();
        }
      else {
          F1 = B();
          if (MA2 < MC2) F2 = A();
          else F2 = C();
        }
    }
  else _folded = 0;
  
  //ExEnv::outn() << scprintf("r = %14.8f, h = %14.8f\n",r(),h);
  //M.print();
  //D[0].print();
  //D[1].print();
  //A().print();
  //B().print();
  //C().print();

  int i;
  for (i=0; i<2; i++) {
      SCVector3 AD = A() - D[i];
      SCVector3 BD = B() - D[i];
      SCVector3 CD = C() - D[i];
      BDxCD[i] = BD.cross(CD);
      BDxCDdotAD[i] = BDxCD[i].dot(AD);
      CDxAD[i] = CD.cross(AD);
      CDxADdotBD[i] = CDxAD[i].dot(BD);
      ADxBD[i] = AD.cross(BD);
      ADxBDdotCD[i] = ADxBD[i].dot(CD);
    }

  for (i=0; i<2; i++) MD[i] = M - D[i];

  // reentrant surfaces need a whole bunch more to be able to compute
  // the distance to the surface
  if (_reentrant) {
      int i;
      double rMD = MD[0].norm(); // this is the same as rMD[1]
      theta_intersect = M_PI_2 - asin(rMD/r());
      r_intersect = r() * sin(theta_intersect);

      {
        // Does the circle of intersection intersect with BA?
        SCVector3 MA = M - A();
        SCVector3 uBA = B() - A(); uBA.normalize();
        SCVector3 XA = uBA * MA.dot(uBA);
        SCVector3 XM = XA - MA;
        double rXM2 = XM.dot(XM);
        double d_intersect_from_x2 = r_intersect*r_intersect - rXM2;
        if (d_intersect_from_x2 > 0.0) {
            _intersects_AB = 1;
            double tmp = sqrt(d_intersect_from_x2);
            double d_intersect_from_x[2];
            d_intersect_from_x[0] = tmp;
            d_intersect_from_x[1] = -tmp;
            for (i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    IABD[i][j] = XM + d_intersect_from_x[j]*uBA + MD[i];
                  }
              }
          }
        else _intersects_AB = 0;
      }

      {
        // Does the circle of intersection intersect with BC?
        SCVector3 MC = M - C();
        SCVector3 uBC = B() - C(); uBC.normalize();
        SCVector3 XC = uBC * MC.dot(uBC);
        SCVector3 XM = XC - MC;
        double rXM2 = XM.dot(XM);
        double d_intersect_from_x2 = r_intersect*r_intersect - rXM2;
        if (d_intersect_from_x2 > 0.0) {
            _intersects_BC = 1;
            double tmp = sqrt(d_intersect_from_x2);
            double d_intersect_from_x[2];
            d_intersect_from_x[0] = tmp;
            d_intersect_from_x[1] = -tmp;
            for (i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    IBCD[i][j] = XM + d_intersect_from_x[j]*uBC + MD[i];
                  }
              }
          }
        else _intersects_BC = 0;
      }

      {
        // Does the circle of intersection intersect with CA?
        SCVector3 MA = M - A();
        SCVector3 uCA = C() - A(); uCA.normalize();
        SCVector3 XA = uCA * MA.dot(uCA);
        SCVector3 XM = XA - MA;
        double rXM2 = XM.dot(XM);
        double d_intersect_from_x2 = r_intersect*r_intersect - rXM2;
        if (d_intersect_from_x2 > 0.0) {
            _intersects_CA = 1;
            double tmp = sqrt(d_intersect_from_x2);
            double d_intersect_from_x[2];
            d_intersect_from_x[0] = tmp;
            d_intersect_from_x[1] = -tmp;
            for (i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    ICAD[i][j] = XM + d_intersect_from_x[j]*uCA + MD[i];
                  }
              }
          }
        else _intersects_CA = 0;
      }

    }

#if 0 // test code
  ExEnv::outn() << "Uncapped5SphereExclusionShape: running some tests" << endl;
  verbose = 1;

  FILE* testout = fopen("testout.vect", "w");

  const double scalefactor_inc = 0.05;
  const double start = -10.0;
  const double end = 10.0;

  SCVector3 middle = (1.0/3.0)*(A()+B()+C());

  int nlines = 1;
  int nvert = (int) ( (end-start) / scalefactor_inc);
  int ncolor = nvert;

  fprintf(testout, "VECT\n%d %d %d\n", nlines, nvert, ncolor);

  fprintf(testout, "%d\n", nvert);
  fprintf(testout, "%d\n", ncolor);

  double scalefactor = start;
  for (int ii = 0; ii<nvert; ii++) {
      SCVector3 position = (D[0] - middle) * scalefactor + middle;
      double d = distance_to_surface(position);
      fprintf(testout, "%f %f %f # value = %f\n",
              position[0], position[1], position[2], d);
      scalefactor += scalefactor_inc;
    }
  scalefactor = start;
  for (ii = 0; ii<nvert; ii++) {
      SCVector3 position = (D[0] - middle) * scalefactor + middle;
      double d = distance_to_surface(position);
      ExEnv::outn() << scprintf("d = %f\n", d);
      if (d<0.0) fprintf(testout,"1.0 0.0 0.0 0.5\n");
      else fprintf(testout,"0.0 0.0 1.0 0.5\n");
      scalefactor += scalefactor_inc;
    }

  fclose(testout);
  ExEnv::outn() << "testout.vect written" << endl;

  verbose = 0;
#endif // test code
  
}
int
Uncapped5SphereExclusionShape::is_outside(const SCVector3&X) const
{
  SCVector3 Xv(X);

  if (verbose) ExEnv::outn() << scprintf("point %14.8f %14.8f %14.8f\n",X(0),X(1),X(2));

  // The folded case isn't handled correctly here, so use
  // the less efficient distance_to_surface routine.
  if (_folded) {
      return distance_to_surface(X) >= 0.0;
    }

  for (int i=0; i<2; i++) {
      SCVector3 XD = Xv - D[i];
      double rXD = XD.norm();
      if (rXD <= r()) return 1;
      double u = BDxCD[i].dot(XD)/BDxCDdotAD[i];
      if (u <= 0.0) return 1;
      double v = CDxAD[i].dot(XD)/CDxADdotBD[i];
      if (v <= 0.0) return 1;
      double w = ADxBD[i].dot(XD)/ADxBDdotCD[i];
      if (w <= 0.0) return 1;
    }

  if (verbose) ExEnv::outn() << "is_inside" << endl;

  return 0;
}
static int
is_contained_in_unbounded_pyramid(SCVector3 XD,
                                  SCVector3 AD,
                                  SCVector3 BD,
                                  SCVector3 CD)
{
  SCVector3 BDxCD = BD.cross(CD);
  SCVector3 CDxAD = CD.cross(AD);
  SCVector3 ADxBD = AD.cross(BD);
  double u = BDxCD.dot(XD)/BDxCD.dot(AD);
  if (u <= 0.0) return 0;
  double v = CDxAD.dot(XD)/CDxAD.dot(BD);
  if (v <= 0.0) return 0;
  double w = ADxBD.dot(XD)/ADxBD.dot(CD);
  if (w <= 0.0) return 0;
  return 1;
}
double
Uncapped5SphereExclusionShape::
  distance_to_surface(const SCVector3&X,SCVector3*grad) const
{
  SCVector3 Xv(X);

  // Find out if I'm on the D[0] side or the D[1] side of the ABC plane.
  int side;
  SCVector3 XM = Xv - M;
  if (MD[0].dot(XM) > 0.0) side = 1;
  else side = 0;

  if (verbose) {
      ExEnv::outn() << scprintf("distance_to_surface: folded = %d, side = %d\n",
                       _folded, side);
      ExEnv::outn() << "XM = "; XM.print();
      ExEnv::outn() << "MD[0] = "; MD[0].print();
      ExEnv::outn() << "MD[0].dot(XM) = " << MD[0].dot(XM) << endl;
    }

  SCVector3 XD = Xv - D[side];
  double u = BDxCD[side].dot(XD)/BDxCDdotAD[side];
  if (verbose) ExEnv::outn() << scprintf("u = %14.8f\n", u);
  if (u <= 0.0) return shape_infinity;
  double v = CDxAD[side].dot(XD)/CDxADdotBD[side];
  if (verbose) ExEnv::outn() << scprintf("v = %14.8f\n", v);
  if (v <= 0.0) return shape_infinity;
  double w = ADxBD[side].dot(XD)/ADxBDdotCD[side];
  if (verbose) ExEnv::outn() << scprintf("w = %14.8f\n", w);
  if (w <= 0.0) return shape_infinity;
  double rXD = XD.norm();
  if (verbose) ExEnv::outn() << scprintf("r() - rXD = %14.8f\n", r() - rXD);
  if (rXD <= r()) {
      if (!_reentrant) return r() - rXD;
      // this shape is reentrant
      // make sure that we're on the right side
      if ((side == 1) || (u + v + w <= 1.0)) {
          // see if we're outside the cone that intersects
          // the opposite sphere
          double angle = acos(fabs(XD.dot(MD[side]))
                              /(XD.norm()*MD[side].norm()));
          if (angle >= theta_intersect) {
              if (grad) {
                  *grad = (-1.0/rXD)*XD;
                }
              return r() - rXD;
            }
          if (_intersects_AB
              &&is_contained_in_unbounded_pyramid(XD,
                                                  MD[side],
                                                  IABD[side][0],
                                                  IABD[side][1])) {
              //ExEnv::outn() << scprintf("XD: "); XD.print();
              //ExEnv::outn() << scprintf("MD[%d]: ",i); MD[i].print();
              //ExEnv::outn() << scprintf("IABD[%d][0]: ",i); IABD[i][0].print();
              //ExEnv::outn() << scprintf("IABD[%d][1]: ",i); IABD[i][1].print();
              return closest_distance(XD,(SCVector3*)IABD[side],2,grad);
            }
          if (_intersects_BC
              &&is_contained_in_unbounded_pyramid(XD,
                                                  MD[side],
                                                  IBCD[side][0],
                                                  IBCD[side][1])) {
              return closest_distance(XD,(SCVector3*)IBCD[side],2,grad);
            }
          if (_intersects_CA
              &&is_contained_in_unbounded_pyramid(XD,
                                                  MD[side],
                                                  ICAD[side][0],
                                                  ICAD[side][1])) {
              return closest_distance(XD,(SCVector3*)ICAD[side],2,grad);
            }
          // at this point we are closest to the ring formed
          // by the intersection of the two probe spheres
          double distance_to_plane;
          double distance_to_ring_in_plane;
          double MDnorm = MD[side].norm();
          SCVector3 XM = XD - MD[side];
          SCVector3 XM_in_plane;
          if (MDnorm<1.0e-16) {
              distance_to_plane = 0.0;
              XM_in_plane = XD;
            }
          else {
              distance_to_plane = XM.dot(MD[side])/MDnorm;
              XM_in_plane = XM - (distance_to_plane/MDnorm)*MD[side];
            }
          if (grad) {
              double XM_in_plane_norm = XM_in_plane.norm();
              if (XM_in_plane_norm < 1.e-8) {
                  // the grad points along MD
                  double scale = - distance_to_plane
                         /(MDnorm*sqrt(r_intersect*r_intersect
                                       + distance_to_plane*distance_to_plane));
                  *grad = MD[side] * scale;
                }
              else {
                  SCVector3 point_on_ring;
                  point_on_ring = XM_in_plane
                                * (r_intersect/XM_in_plane_norm) + M;
                  SCVector3 gradv = Xv - point_on_ring;
                  gradv.normalize();
                  *grad = gradv;
                }
            }
          distance_to_ring_in_plane =
                         r_intersect - sqrt(XM_in_plane.dot(XM_in_plane));
          return sqrt(distance_to_ring_in_plane*distance_to_ring_in_plane
                      +distance_to_plane*distance_to_plane);
        }
    }

  if (verbose) ExEnv::outn() << "returning -1.0" << endl;
  return -1.0;
}

void
Uncapped5SphereExclusionShape::boundingbox(double valuemin, double valuemax,
                                           SCVector3& p1,
                                           SCVector3& p2)
{
  SCVector3 p11;
  SCVector3 p12;
  SCVector3 p21;
  SCVector3 p22;
  SCVector3 p31;
  SCVector3 p32;

  _s1.boundingbox(valuemin,valuemax,p11,p12);
  _s2.boundingbox(valuemin,valuemax,p21,p22);
  _s3.boundingbox(valuemin,valuemax,p31,p32);

  int i;
  for (i=0; i<3; i++) {
      if (p11[i] < p21[i]) p1[i] = p11[i];
      else p1[i] = p21[i];
      if (p31[i] < p1[i]) p1[i] = p31[i];
      if (p12[i] > p22[i]) p2[i] = p12[i];
      else p2[i] = p22[i];
      if (p32[i] > p2[i]) p2[i] = p32[i];
    }
}

int
Uncapped5SphereExclusionShape::gradient_implemented() const
{
  return 1;
}

/////////////////////////////////////////////////////////////////////
// Unionshape

static ClassDesc UnionShape_cd(
  typeid(UnionShape),"UnionShape",1,"public Shape",
  0, 0, 0);

UnionShape::UnionShape()
{
}

UnionShape::~UnionShape()
{
}

void
UnionShape::add_shape(Ref<Shape> s)
{
  _shapes.insert(s);
}

// NOTE: this underestimates the distance to the surface when
//inside the surface
double
UnionShape::distance_to_surface(const SCVector3&p,SCVector3* grad) const
{
  std::set<Ref<Shape> >::const_iterator imin = _shapes.begin();
  if (imin == _shapes.end()) return 0.0;
  double min = (*imin)->distance_to_surface(p);
  for (std::set<Ref<Shape> >::const_iterator i=imin; i!=_shapes.end(); i++) {
      double d = (*i)->distance_to_surface(p);
      if (min <= 0.0) {
          if (d < 0.0 && d > min) { min = d; imin = i; }
        }
      else {
          if (min > d) { min = d; imin = i; }
        }
    }

  if (grad) {
      (*imin)->distance_to_surface(p,grad);
    }
  return min;
}

int
UnionShape::is_outside(const SCVector3&p) const
{
  for (std::set<Ref<Shape> >::const_iterator i=_shapes.begin();
       i!=_shapes.end(); i++) {
      if (!(*i)->is_outside(p)) return 0;
    }

  return 1;
}

void
UnionShape::boundingbox(double valuemin, double valuemax,
                        SCVector3& p1,
                        SCVector3& p2)
{
  if (_shapes.begin() == _shapes.end()) {
      for (int i=0; i<3; i++) p1[i] = p2[i] = 0.0;
      return;
    }
  
  SCVector3 pt1;
  SCVector3 pt2;
  
  std::set<Ref<Shape> >::iterator j = _shapes.begin();
  int i;
  (*j)->boundingbox(valuemin,valuemax,p1,p2);
  for (j++; j!=_shapes.end(); j++) {
      (*j)->boundingbox(valuemin,valuemax,pt1,pt2);
      for (i=0; i<3; i++) {
          if (pt1[i] < p1[i]) p1[i] = pt1[i];
          if (pt2[i] > p2[i]) p2[i] = pt2[i];
        }
    }
}

int
UnionShape::gradient_implemented() const
{
  for (std::set<Ref<Shape> >::const_iterator j=_shapes.begin();
       j!=_shapes.end(); j++) {
      if (!(*j)->gradient_implemented()) return 0;
    }
  return 1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
