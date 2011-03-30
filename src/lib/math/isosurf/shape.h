//
// shape.h
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

#ifndef _math_isosurf_shape_h
#define _math_isosurf_shape_h

#include <set>

#include <math/isosurf/volume.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>

namespace sc {

/** A Shape is a Volume represents an 3D solid.  The value of the Shape at
each point in space is the distance to the surface.  The distance is
negative if the point is inside the solid.  For Shape specializations that
cannot compute the distance to the surface, the value will be 1.0 outside
and -1.0 inside the solid. */
class Shape: public Volume {
  public:
    Shape();
    Shape(const Ref<KeyVal>&keyval);
    virtual double distance_to_surface(const SCVector3&r,
                                       SCVector3*grad=0) const = 0;
    virtual int is_outside(const SCVector3&r) const;
    virtual ~Shape();
    void compute();
    void interpolate(const SCVector3& p1,
                     const SCVector3& p2,
                     double val,
                     SCVector3& result);

    int value_implemented() const;
};



class SphereShape: public Shape {
  private:
    SCVector3 _origin;
    double _radius;
  public:
    SphereShape(const SCVector3&,double);
    SphereShape(const Ref<KeyVal>&);
    SphereShape(const SphereShape&);
    ~SphereShape();
    void boundingbox(double minvalue, double maxvalue,
                     SCVector3& p1, SCVector3&p2);
    double radius() const { return _radius; }
    const SCVector3& origin() const { return _origin; }
    double distance_to_surface(const SCVector3&r,SCVector3*grad=0) const;
    void print(std::ostream&o=ExEnv::out0()) const;

    // these are used to update the parameters describing the sphere
    double radius(double r);
    const SCVector3& origin(const SCVector3& o);

    int gradient_implemented() const;
};

inline double
SphereShape::radius(double r)
{
  obsolete();
  return _radius = r;
}

inline const SCVector3&
SphereShape::origin(const SCVector3& o)
{
  obsolete();
  _origin = o;
  return _origin;
}

class UncappedTorusHoleShape: public Shape
{
  private:
    SphereShape _s1;
    SphereShape _s2;
    double _r;
  protected:
    void in_plane_sphere(const SCVector3& point,
                         SCVector3& origin) const;
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
    void print(std::ostream&o=ExEnv::out0()) const;
    void boundingbox(double valuemin, double valuemax,
                     SCVector3& p1, SCVector3&p2);

    int gradient_implemented() const;
};

class NonreentrantUncappedTorusHoleShape: public UncappedTorusHoleShape
{
  private:
    double rAP;
    double rBP;
    SCVector3 BA;
  public:
    NonreentrantUncappedTorusHoleShape(double r,
                                       const SphereShape&,
                                       const SphereShape&);
    ~NonreentrantUncappedTorusHoleShape();
    double distance_to_surface(const SCVector3&r,SCVector3*grad=0) const;

    int gradient_implemented() const;
};

class ReentrantUncappedTorusHoleShape: public UncappedTorusHoleShape
{
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
    double distance_to_surface(const SCVector3&r,SCVector3*grad=0) const;

    int gradient_implemented() const;
};

class Uncapped5SphereExclusionShape: public Shape
{
  private:
    int _solution_exists;
    int _reentrant;
    int _folded;
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

    // these are needed for folded shapes
    // F1 and F2 are the two points of A, B, and C that are closed to M
    SCVector3 F1;
    SCVector3 F2;
    
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
    double distance_to_surface(const SCVector3&r,SCVector3*grad=0) const;
    int is_outside(const SCVector3&) const;
    void boundingbox(double valuemin, double valuemax,
                     SCVector3& p1, SCVector3&p2);

    int gradient_implemented() const;
};

/** A UnionShape is volume enclosed by a set of Shape's. */
class UnionShape: public Shape {
  protected:
    std::set<Ref<Shape> > _shapes;
  public:
    void add_shape(Ref<Shape>);
    UnionShape();
    ~UnionShape();
    double distance_to_surface(const SCVector3&r,SCVector3*grad=0) const;
    int is_outside(const SCVector3&r) const;
    void boundingbox(double valuemin, double valuemax,
                     SCVector3& p1, SCVector3& p2);

    int gradient_implemented() const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
