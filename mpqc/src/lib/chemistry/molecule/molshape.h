
#ifndef _chemistry_molecule_molshape_h
#define _chemistry_molecule_molshape_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/shape.h>
#include <chemistry/molecule/molinfo.h>
#include <chemistry/molecule/molecule.h>

class VDWShape: public UnionShape {
#   define CLASSNAME VDWShape
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 public:
    VDWShape(const RefMolecule&);
    VDWShape(const RefKeyVal&);
    ~VDWShape();
    void initialize(const RefMolecule&);
};  

class ConnollyShape: public UnionShape {
#   define CLASSNAME ConnollyShape
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefAtomInfo atominfo_;
 public:
    ConnollyShape(const RefKeyVal&);
    ~ConnollyShape();
    void initialize(const RefMolecule&,double probe_radius);
};

#ifndef COUNT_CONNOLLY2
# define COUNT_CONNOLLY2 1
#endif

// This is a utility class needed by ConnollyShape2
class CS2Sphere
{
    SCVector3 _v;
    double _radius;

  public:
#if COUNT_CONNOLLY2
    static int n_no_spheres_;
    static int n_probe_enclosed_by_a_sphere_;
    static int n_probe_center_not_enclosed_;
    static int n_surface_of_s0_not_covered_;
    static int n_plane_totally_covered_;
    static int n_internal_edge_not_covered_;
    static int n_totally_covered_;
#endif

    CS2Sphere(const SCVector3& v, double rad):
    _v(v),_radius(rad){}
    CS2Sphere(double x, double y, double z, double rad):
    _v(x,y,z),_radius(rad){}
    CS2Sphere(void) {};
    void initialize(SCVector3& v, double rad) {
        _v = v; _radius = rad; }

    CS2Sphere& operator=(const CS2Sphere&s) {
        _v = s._v; _radius = s._radius; return *this; }
    
    // Return the distance between the centers of the two
    // spheres
    double distance(CS2Sphere &asphere)
    { return sqrt((_v[0]-asphere._v[0])*(_v[0]-asphere._v[0])+
                  (_v[1]-asphere._v[1])*(_v[1]-asphere._v[1])+
                  (_v[2]-asphere._v[2])*(_v[2]-asphere._v[2]));}  

    // Return the radius of the circle intersecting the two spheres
    // Note that this assumes the spheres do overlap!
    double common_radius(CS2Sphere &asphere);

    // Return the center
    const SCVector3& center(void) const { return _v; }
    double x() const { return _v[0]; }
    double y() const { return _v[1]; }
    double z() const { return _v[2]; }

    // Return the vector3d connecting the two centers
    SCVector3 center_vec(const CS2Sphere &asphere) { return _v - asphere._v; }

    double radius(void) const {return _radius;}

    void recenter(const SCVector3 &v) { _v -= v; }
    void print(FILE* fp=stdout) 
    {
        fprintf(fp,"Rad=%lf, Center=(%lf,%lf,%lf), From origin=%lf\n",
                _radius, _v[0], _v[1], _v[2], _v.norm());
    }

    // Function to determine if there is any portion of this that 
    // is not inside one or more of the spheres in s[].  Returns
    // 1 if the intersection is empty, otherwise 0 is returned.
    // Warning: the spheres in s are modified.
    int intersect(CS2Sphere *s,
                  int n_spheres) const;

    static void print_counts(FILE*fp = stdout);
};

class ConnollyShape2: public Shape {
#   define CLASSNAME ConnollyShape2
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    CS2Sphere* sphere;
    double probe_r;
    int n_spheres;
    RefAtomInfo atominfo_;

#if COUNT_CONNOLLY2
    static int n_total_;
    static int n_inside_vdw_;
    const int n_with_nsphere_dim_ = 10;
    static int n_with_nsphere_[n_with_nsphere_dim_];
#endif

  public:
    ConnollyShape2(const RefKeyVal&);
    ~ConnollyShape2();
    void initialize(const RefMolecule&,double probe_radius);
    double distance_to_surface(const SCVector3&r,
                               double*grad=0) const;
    void boundingbox(double valuemin,
                     double valuemax,
                     RefSCVector& p1, RefSCVector& p2);

    static void print_counts(FILE*fp = stdout);
};

#endif
