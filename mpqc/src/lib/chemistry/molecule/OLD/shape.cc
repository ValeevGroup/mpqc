
extern "C" {
# include <math.h>
  }

#include "shape.h"
#include "molecule.h"

REF_def(Shape);
ARRAY_def(RefShape);
SET_def(RefShape);
ARRAYSET_def(RefShape);

Shape::~Shape()
{
}

int Shape::is_outside(const Point3&r) const
{
  if (distance_to_surface(r)>0.0) return 1;
  return 0;
}

SphereShape::SphereShape(const Point3&o,double r):
  _origin(o),
  _radius(r)
{
}

SphereShape::~SphereShape()
{
}

REF_def(SphereShape);
ARRAY_def(RefSphereShape);
SET_def(RefSphereShape);
ARRAYSET_def(RefSphereShape);

double SphereShape::distance_to_surface(const Point3&p) const
{
  int i;
  double r2 = 0.0;
  for (i=0; i<3; i++) {
      double tmp = p[i] - _origin[i];
      r2 += tmp*tmp;
    }
  double r = sqrt(r2);
  double d = r - _radius;
  if (d < 0.0) return -1.0;
  return d;
}

void SphereShape::print(FILE*fp) const
{
  fprintf(fp,"SphereShape: r = %8.4f o = (%8.4f %8.4f %8.4f)\n",
          radius(),origin()[0],origin()[1],origin()[2]);
}

////////////////////////////////////////////////////////////////////////
// UncappedTorusHoleShape

UncappedTorusHoleShape::UncappedTorusHoleShape(double r,
                               const SphereShape& s1,
                               const SphereShape& s2):
_r(r),
_s1(s1),
_s2(s2)
{
}

UncappedTorusHoleShape*
UncappedTorusHoleShape::newUncappedTorusHoleShape(double r,
                                                  const SphereShape&s1,
                                                  const SphereShape&s2)
{
  // if the probe sphere fits between the two spheres, then there
  // is no need to include this shape
  DVector3 A(s1.origin());
  DVector3 B(s2.origin());
  DVector3 BA = B - A;
  if (2.0*r <  BA.norm() - s1.radius() - s2.radius()) return 0;

  // check to see if the surface is reentrant
  double rrs1 = r+s1.radius();
  double rrs2 = r+s2.radius();
  DVector3 R12 = ((DVector3)s1.origin()) - ((DVector3)s2.origin());
  double r12 = sqrt(R12.dot(R12));
  if (sqrt(rrs1*rrs1-r*r) + sqrt(rrs2*rrs2-r*r) < r12)
    return new ReentrantUncappedTorusHoleShape(r,s1,s2);

  // otherwise create an ordinary torus hole
  return new NonreentrantUncappedTorusHoleShape(r,s1,s2);
}

// Given a node, finds a sphere in the plane of n and the centers
// of _s1 and _s2 that touches the UncappedTorusHole.  There are two
// candidates, the one closest to n is chosen.
SphereShape UncappedTorusHoleShape::in_plane_sphere(const Point3& n) const
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

  double r_a = _s1.radius();
  double r_b = _s2.radius();
  DVector3 A = _s1.origin();
  DVector3 B = _s2.origin();
  DVector3 N = n;
  DVector3 R_AB = B - A;
  DVector3 R_AN = N - A;

  // vector projection of R_AN onto R_AB
  DVector3 P_AN_AB = R_AB * (R_AN.dot(R_AB)/R_AB.dot(R_AB));
  // the perpendicular vector
  DVector3 U_perp = R_AN - P_AN_AB;

  // if |U| is tiny, then any vector perp to AB will do
  double u = U_perp.dot(U_perp);
  if (u<1.0e-23) {
      DVector3 try = R_AB;
      try[0] += 1.0;
      try = try - R_AB * (try.dot(R_AB)/R_AB.dot(R_AB));
      if (try.dot(try) < 1.0e-23) {
          try = R_AB;
          try[1] += 1.0;
          try = try - R_AB * (try.dot(R_AB)/R_AB.dot(R_AB));
        }
      U_perp = try;
    }

  U_perp.normalize();
  //printf("A: "); A.print();
  //printf("U_perp: "); U_perp.print();
  //printf("R_AB: "); R_AB.print();

  double r_AB = sqrt(R_AB.dot(R_AB));
  double r_A = _s1.radius();
  double r_B = _s2.radius();

  double r_Ar = r_A + _r;
  double r_Br = r_B + _r;
  double a = ((r_Ar*r_Ar - r_Br*r_Br)/(2.0*r_AB*r_AB)) + 0.5;
  double b = sqrt(r_Ar*r_Ar - a*a*r_AB*r_AB);

  //printf("r_Ar = %f, r_AB = %f\n",r_Ar,r_AB);
  //printf("a = %f, b = %f\n",a,b);

  Point3 P = A + a * R_AB + b * U_perp;
  //printf("a*R_AB: "); (a*R_AB).print();
  //printf("b*U_perp: "); (b*U_perp).print();

  SphereShape s(P,_r);

  return s;
}

void UncappedTorusHoleShape::print(FILE*fp) const
{
  fprintf(fp,"UncappedTorusHoleShape:\n");
  fprintf(fp,"  r = %8.5f\n",_r);
  fprintf(fp,"  s1 = ");
  _s1.print(fp);
  fprintf(fp,"  s2 = ");
  _s2.print(fp);
}

/////////////////////////////////////////////////////////////////////
// is in triangle

static int
is_in_triangle(const DVector3&XA,const DVector3&PA,const DVector3&BA)
{
  DVector3 axis = BA.cross(PA);

  DVector3 BA_perp = BA; BA_perp.rotate(M_PI_2,axis);
  double u = BA_perp.dot(XA)/BA_perp.dot(PA);
  if (u<0.0 || u>1.0) return 0;

  DVector3 PA_perp = PA; PA_perp.rotate(M_PI_2,axis);
  double w = PA_perp.dot(XA)/PA_perp.dot(BA);
  if (w<0.0 || w>1.0) return 0;

  double uw = u+w;
  if (uw<0.0 || uw>1.0) return 0;
  return 1;
}

/////////////////////////////////////////////////////////////////////
// ReentrantUncappedTorusHoleShape

ReentrantUncappedTorusHoleShape::ReentrantUncappedTorusHoleShape(double r,
                                                 const SphereShape& s1,
                                                 const SphereShape& s2):
  UncappedTorusHoleShape(r,s1,s2)
{
  rAP = r + s1.radius();
  rBP = r + s2.radius();
  BA = B() - A();

  // Find the points at the ends of the two cone-like objects, Z1 and Z2.
  // they are given by:
  //   Z = z BA, where BA = B-A
  //   r^2 = PZ.PZ, where PZ = P-Z and P is the center of a probe sphere
  // this gives
  //  z^2 BA.BA - 2z P.BA + P.P - r^2 = 0

// I don't need this unless I want distance to surface implemented
//   Point3 arbitrary; 
//   arbitrary[0] = arbitrary[1] = arbitrary[2] = 0.0;
//   SphereShape PS = in_plane_sphere(arbitrary);
//   DVector3 P(PS.origin());
// 
//   double a = BA.dot(BA);
//   double minus_b = 2.0 * P.dot(BA);
//   double c = P.dot(P) - r*r;
//   double b2m4ac = minus_b*minus_b - 4*a*c;
//   double sb2m4ac = sqrt(b2m4ac);
//   if (b2m4ac >= 0.0) {
//       sb2m4ac = sqrt(b2m4ac);
//     }
//   else if (b2m4ac > -1.0e-10) {
//       sb2m4ac = 0.0;
//     }
//   else {
//       fprintf(stderr,"ReentrantUncappedTorusHoleShape:: imaginary point\n");
//       abort();
//     }
//   zA = (minus_b - sb2m4ac)/(2.0*a);
//   zB = (minus_b + sb2m4ac)/(2.0*a);
//   ZA = BA * zA;
//   ZB = BA * zB;
}
ReentrantUncappedTorusHoleShape::~ReentrantUncappedTorusHoleShape()
{
}
int ReentrantUncappedTorusHoleShape::
  is_outside(const Point3&X) const
{
  DVector3 Xv(X);

  DVector3 P = in_plane_sphere(X).origin();
  DVector3 PX = P - Xv;
  double rPX = PX.norm();
  if (rPX > rAP || rPX > rBP) return 1;

  DVector3 PA = P - A();

  DVector3 XA = Xv - A();

  if (!is_in_triangle(XA,PA,BA)) return 1;

  if (rPX < radius()) {
      return 1;
    }

  return 0;
}
double ReentrantUncappedTorusHoleShape::
  distance_to_surface(const Point3&r) const
{
  fprintf(stderr,"ReentrantUncappedTorusHoleShape::"
          "distance_to_surface: can't compute\n");
  abort();
  return 0.0;
}

/////////////////////////////////////////////////////////////////////
// NonreentrantUncappedTorusHoleShape

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
  distance_to_surface(const Point3&X) const
{
  DVector3 Xv(X);

  DVector3 P = in_plane_sphere(X).origin();
  DVector3 PX = P - Xv;
  double rPX = PX.norm();
  if (rPX > rAP || rPX > rBP) return infinity;

  DVector3 PA = P - A();
  DVector3 XA = Xv - A();

  DVector3 axis = BA.cross(PA);

  DVector3 BA_perp = BA; BA_perp.rotate(M_PI_2,axis);
  double u = BA_perp.dot(XA)/BA_perp.dot(PA);
  if (u<0.0 || u>1.0) return infinity;

  DVector3 PA_perp = PA; PA_perp.rotate(M_PI_2,axis);
  double w = PA_perp.dot(XA)/PA_perp.dot(BA);
  if (w<0.0 || w>1.0) return infinity;

  double uw = u+w;
  if (uw<0.0 || uw>1.0) return infinity;

  if (rPX < radius()) return radius() - rPX;

  return -1;
}

/////////////////////////////////////////////////////////////////////
// Uncapped5SphereExclusionShape

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
Uncapped5SphereExclusionShape::
  Uncapped5SphereExclusionShape(double radius,
                                const SphereShape&s1,
                                const SphereShape&s2,
                                const SphereShape&s3):
  _r(radius),
  _s1(s1),
  _s2(s2),
  _s3(s3)
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
  DVector3 BA = B()-A();
  double DdotBA = 0.5*(rAr2 - rBr2 + B2 - A2);
  double DAdotBA = DdotBA - A().dot(BA);
  double BA2 = BA.dot(BA);
  DVector3 CA = C() - A();
  double CAdotBA = CA.dot(BA);
  DVector3 CA_perpBA = CA - (CAdotBA/BA2)*BA;
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
  if (h<r()) _reentrant = 1;
  else _reentrant = 0;

  D[0] = A() + (DAdotBA/BA2)*BA + (DAdotCA_perpBA/CA_perpBA2)*CA_perpBA;
  D[1] = D[0];
  DVector3 BAxCA = BA.cross(CA);
  D[0] = D[0] + h * BAxCA * ( 1.0/BAxCA.norm() );
  D[1] = D[1] - h * BAxCA * ( 1.0/BAxCA.norm() );

  D[0].print();
  D[1].print();
  A().print();
  B().print();
  C().print();

  for (int i=0; i<2; i++) {
      DVector3 AD = A() - D[i];
      DVector3 BD = B() - D[i];
      DVector3 CD = C() - D[i];
      BDxCD[i] = BD.cross(CD);
      BDxCDdotAD[i] = BDxCD[i].dot(AD);
      CDxAD[i] = CD.cross(AD);
      CDxADdotBD[i] = CDxAD[i].dot(BD);
      ADxBD[i] = AD.cross(BD);
      ADxBDdotCD[i] = ADxBD[i].dot(CD);
    }
}
int Uncapped5SphereExclusionShape::
  is_outside(const Point3&X) const
{
  DVector3 Xv(X);

  for (int i=0; i<2; i++) {
      DVector3 XD = Xv - D[i];
      double rXD = XD.norm();
      if (rXD <= r()) return 1;
      double u = BDxCD[i].dot(XD)/BDxCDdotAD[i];
      if (u <= 0.0) return 1;
      double v = CDxAD[i].dot(XD)/CDxADdotBD[i];
      if (v <= 0.0) return 1;
      double w = ADxBD[i].dot(XD)/ADxBDdotCD[i];
      if (w <= 0.0) return 1;
    }

  return 0;
}
double Uncapped5SphereExclusionShape::
  distance_to_surface(const Point3&r) const
{
  fprintf(stderr,"Uncapped5SphereExclusionShape::"
          "distance_to_surface: can't compute\n");
  abort();
  return 0.0;
}


/////////////////////////////////////////////////////////////////////
// Unionshape

UnionShape::UnionShape()
{
}

UnionShape::~UnionShape()
{
}

void UnionShape::add_shape(RefShape s)
{
  _shapes.add(s);
}

double UnionShape::distance_to_surface(const Point3&p) const
{
  if (_shapes.length() == 0) return 0.0;
  double min = _shapes[0]->distance_to_surface(p);
  if (min < 0.0) return -1.0;
  for (int i=1; i<_shapes.length(); i++) {
      double d = _shapes[i]->distance_to_surface(p);
      if (d < 0.0) return -1.0;
      if (min > d) min = d;
    }

  return min;
}

int UnionShape::is_outside(const Point3&p) const
{
  if (_shapes.length() == 0) return 1;
  for (int i=0; i<_shapes.length(); i++) {
      if (!_shapes[i]->is_outside(p)) return 0;
    }

  return 1;
}

VDWShape::VDWShape(Molecule&mol)
{
  for (int i=0; i<mol.natom(); i++) {
      Point3 r;
      for (int j=0; j<3; j++) r[j] = mol[i][j];
      add_shape(new SphereShape(r,mol[i].element().atomic_radius_au()));
    }
}

VDWShape::~VDWShape()
{
}

ConnollyShape::ConnollyShape(Molecule&mol,double probe_radius)
{
  ArraySetRefSphereShape spheres;
  for (int i=0; i<mol.natom(); i++) {
      Point3 r;
      for (int j=0; j<3; j++) r[j] = mol[i][j];
      RefSphereShape
        sphere(new SphereShape(r,mol[i].element().atomic_radius_au()));
      add_shape(sphere.pointer());
      spheres.add(sphere);
    }

  for (i=0; i<spheres.length(); i++) {
      for (int j=0; j<i; j++) {
          RefShape th =
            UncappedTorusHoleShape::newUncappedTorusHoleShape(probe_radius,
                                              *(spheres[i].pointer()),
                                              *(spheres[j].pointer()));
          if (th.null()) continue;
          add_shape(th);
          // now check for excluding volume for groups of three spheres
          for (int k=0; k<j; k++) {
              RefShape e =
                Uncapped5SphereExclusionShape::
              newUncapped5SphereExclusionShape(probe_radius,
                                               *(spheres[i].pointer()),
                                               *(spheres[j].pointer()),
                                               *(spheres[k].pointer()));
              if (e.nonnull()) add_shape(e);
            }
        }
    }
}

ConnollyShape::~ConnollyShape()
{
}
