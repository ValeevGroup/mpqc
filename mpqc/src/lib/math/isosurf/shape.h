
#ifndef _math_isosurf_shape_h
#define _math_isosurf_shape_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/volume.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <util/container/ref.h>
#include <util/container/array.h>
#include <util/container/set.h>

static const double infinity = 1.0e23;

class Shape: public Volume {
#   define CLASSNAME Shape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    
    Shape();
    virtual double distance_to_surface(const SCVector3&r,
                                       double*grad=0) const = 0;
    virtual int is_outside(const SCVector3&r) const;
    virtual ~Shape();
    void compute();
    RefSCVector interpolate(RefSCVector& p1,RefSCVector& p2,double val);
};

SavableState_REF_dec(Shape);
ARRAY_dec(RefShape);
SET_dec(RefShape);
ARRAYSET_dec(RefShape);

class SphereShape: public Shape {
#   define CLASSNAME SphereShape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    SCVector3 _origin;
    double _radius;
  public:
    SphereShape(const SCVector3&,double);
    SphereShape(const SphereShape&);
    ~SphereShape();
    void boundingbox(double minvalue, double maxvalue,
                     RefSCVector& p1, RefSCVector&p2);
    inline double radius() const { return _radius; };
    inline const SCVector3& origin() const { return _origin; };
    double distance_to_surface(const SCVector3&r,double*grad=0) const;
    void print(FILE*fp=stdout) const;
};

REF_dec(SphereShape);
ARRAY_dec(RefSphereShape);
SET_dec(RefSphereShape);
ARRAYSET_dec(RefSphereShape);

class UncappedTorusHoleShape: public Shape
{
#   define CLASSNAME UncappedTorusHoleShape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    SphereShape _s1;
    SphereShape _s2;
    double _r;
  protected:
    SphereShape in_plane_sphere(const SCVector3&) const;
    UncappedTorusHoleShape(double r,const SphereShape&,const SphereShape&);
  public:
    static UncappedTorusHoleShape*
    newUncappedTorusHoleShape(double r,
                              const SphereShape&,
                              const SphereShape&);
    inline ~UncappedTorusHoleShape() {};
    inline const SphereShape& sphere(int i) const { return (i?_s2:_s1); };
    inline const SCVector3 A() const { SCVector3 v(_s1.origin()); return v; }
    inline const SCVector3 B() const { SCVector3 v(_s2.origin()); return v; }
    inline double radius() const { return _r; };
    void print(FILE*fp=stdout) const;
    void boundingbox(double valuemin, double valuemax,
                     RefSCVector& p1, RefSCVector&p2);
};

class NonreentrantUncappedTorusHoleShape: public UncappedTorusHoleShape
{
#   define CLASSNAME NonreentrantUncappedTorusHoleShape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    double rAP;
    double rBP;
    SCVector3 BA;
  public:
    NonreentrantUncappedTorusHoleShape(double r,
                                       const SphereShape&,
                                       const SphereShape&);
    ~NonreentrantUncappedTorusHoleShape();
    double distance_to_surface(const SCVector3&r,double*grad=0) const;
};

class ReentrantUncappedTorusHoleShape: public UncappedTorusHoleShape
{
#   define CLASSNAME ReentrantUncappedTorusHoleShape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    double rAP;
    double rBP;
    SCVector3 BA;
    SCVector3 I[2]; // the intersect points
  public:
    ReentrantUncappedTorusHoleShape(double r,
                                    const SphereShape&,
                                    const SphereShape&);
    ~ReentrantUncappedTorusHoleShape();
    int is_outside(const SCVector3&r) const;
    double distance_to_surface(const SCVector3&r,double*grad=0) const;
};

class Uncapped5SphereExclusionShape: public Shape
{
#   define CLASSNAME Uncapped5SphereExclusionShape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    int _solution_exists;
    int _reentrant;
    SphereShape _s1;
    SphereShape _s2;
    SphereShape _s3;
    SCVector3 D[2];
    double BDxCDdotAD[2];
    SCVector3 BDxCD[2];
    double CDxADdotBD[2];
    SCVector3 CDxAD[2];
    double ADxBDdotCD[2];
    SCVector3 ADxBD[2];
    double _r;
    
    // these are needed for reentrant surfaces to compute distances
    SCVector3 M;   // projection of D onto ABC plane
    SCVector3 MD[2];  // M - D 
    double theta_intersect; // angle M - D - intersect_point
    double r_intersect; // the radius of the intersect circle
    int _intersects_AB;
    SCVector3 IABD[2][2];
    int _intersects_BC;
    SCVector3 IBCD[2][2];
    int _intersects_CA;
    SCVector3 ICAD[2][2];
    
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
    inline const SCVector3 A() const { SCVector3 v(_s1.origin()); return v; }
    inline const SCVector3 B() const { SCVector3 v(_s2.origin()); return v; }
    inline const SCVector3 C() const { SCVector3 v(_s3.origin()); return v; }
    inline double rA() const { return _s1.radius(); };
    inline double rB() const { return _s2.radius(); };
    inline double rC() const { return _s3.radius(); };
    inline double r() const { return _r; };
    inline int solution_exists() const { return _solution_exists; };
    void print(FILE*fp=stdout) const;
    double distance_to_surface(const SCVector3&r,double*grad=0) const;
    int is_outside(const SCVector3&) const;
    void boundingbox(double valuemin, double valuemax,
                     RefSCVector& p1, RefSCVector&p2);
};

class UnionShape: public Shape {
#   define CLASSNAME UnionShape
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    ArraysetRefShape _shapes;
  public:
    void add_shape(RefShape);
    UnionShape();
    ~UnionShape();
    double distance_to_surface(const SCVector3&r,double*grad=0) const;
    int is_outside(const SCVector3&r) const;
    void boundingbox(double valuemin, double valuemax,
                     RefSCVector& p1, RefSCVector& p2);
};

#endif
