
#ifndef _chemistry_molecule_shape_h
#define _chemistry_molecule_shape_h

#include <chemistry/molecule/molecule.h>

#include <math/nihmatrix/dvector3.h>
#include <math/topology/point.h>
#include <util/container/ref.h>
#include <util/container/array.h>
#include <util/container/set.h>

class Shape: public RefCount {
 public:
  const double infinity = 1.0e23;

  virtual double distance_to_surface(const Point3&r) const = 0;
  virtual int is_outside(const Point3&r) const;
  virtual ~Shape();
};

REF_dec(Shape);
ARRAY_dec(RefShape);
SET_dec(RefShape);
ARRAYSET_dec(RefShape);

class SphereShape: public Shape {
 private:
  Point3 _origin;
  double _radius;
 public:
  SphereShape(const Point3&,double);
  ~SphereShape();
  inline double radius() const { return _radius; };
  inline const Point3& origin() const { return _origin; };
  double distance_to_surface(const Point3&r) const;
  void print(FILE*fp=stdout) const;
};

REF_dec(SphereShape);
ARRAY_dec(RefSphereShape);
SET_dec(RefSphereShape);
ARRAYSET_dec(RefSphereShape);

class UncappedTorusHoleShape: public Shape
{
 private:
  SphereShape _s1;
  SphereShape _s2;
  double _r;
 protected:
  SphereShape in_plane_sphere(const Point3&) const;
  UncappedTorusHoleShape(double r,const SphereShape&,const SphereShape&);
 public:
  static UncappedTorusHoleShape* newUncappedTorusHoleShape(double r,
                                           const SphereShape&,
                                           const SphereShape&);
  inline ~UncappedTorusHoleShape() {};
  inline const SphereShape& sphere(int i) const { return (i?_s2:_s1); };
  inline const DVector3 A() const { DVector3 v(_s1.origin()); return v; }
  inline const DVector3 B() const { DVector3 v(_s2.origin()); return v; }
  inline double radius() const { return _r; };
  void print(FILE*fp=stdout) const;
};

class NonreentrantUncappedTorusHoleShape: public UncappedTorusHoleShape
{
 private:
  double rAP;
  double rBP;
  DVector3 BA;
 public:
  NonreentrantUncappedTorusHoleShape(double r,
                                     const SphereShape&,
                                     const SphereShape&);
  ~NonreentrantUncappedTorusHoleShape();
  double distance_to_surface(const Point3&r) const;
};

class ReentrantUncappedTorusHoleShape: public UncappedTorusHoleShape
{
 private:
  double rAP;
  double rBP;
  DVector3 BA;
// I don't need these unless I want distance to surface implemented
//   DVector3 ZA;
//   DVector3 ZB;
//   double zA;
//   double zB;
 public:
  ReentrantUncappedTorusHoleShape(double r,
                                  const SphereShape&,
                                  const SphereShape&);
  ~ReentrantUncappedTorusHoleShape();
  int is_outside(const Point3&r) const;
  double distance_to_surface(const Point3&r) const;
};

class Uncapped5SphereExclusionShape: public Shape
{
 private:
  int _solution_exists;
  int _reentrant;
  SphereShape _s1;
  SphereShape _s2;
  SphereShape _s3;
  DVector3 D[2];
  double BDxCDdotAD[2];
  DVector3 BDxCD[2];
  double CDxADdotBD[2];
  DVector3 CDxAD[2];
  double ADxBDdotCD[2];
  DVector3 ADxBD[2];
  double _r;
 protected:
  Uncapped5SphereExclusionShape(double r,
                                const SphereShape&,
                                const SphereShape&,
                                const SphereShape&);
 public:
  static Uncapped5SphereExclusionShape*
  newUncapped5SphereExclusionShape(double r,
                                   const SphereShape&,
                                   const SphereShape&,
                                   const SphereShape&);
  inline ~Uncapped5SphereExclusionShape() {};
  inline const DVector3 A() const { DVector3 v(_s1.origin()); return v; }
  inline const DVector3 B() const { DVector3 v(_s2.origin()); return v; }
  inline const DVector3 C() const { DVector3 v(_s3.origin()); return v; }
  inline double rA() const { return _s1.radius(); };
  inline double rB() const { return _s2.radius(); };
  inline double rC() const { return _s3.radius(); };
  inline double r() const { return _r; };
  inline int solution_exists() const { return _solution_exists; };
  void print(FILE*fp=stdout) const;
  double distance_to_surface(const Point3&r) const;
  int is_outside(const Point3&) const;
};

class UnionShape: public Shape {
 protected:
  ArraySetRefShape _shapes;
 public:
  void add_shape(RefShape);
  UnionShape();
  ~UnionShape();
  double distance_to_surface(const Point3&r) const;
  int is_outside(const Point3&r) const;
};

class VDWShape: public UnionShape {
 public:
  VDWShape(Molecule&);
  ~VDWShape();
};  

class ConnollyShape: public UnionShape {
 public:
  ConnollyShape(Molecule&,double probe_radius = 2.6456173);
  ~ConnollyShape();
};  

#endif
