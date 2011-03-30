//
// molshape.h
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

#ifndef _chemistry_molecule_molshape_h
#define _chemistry_molecule_molshape_h

#include <util/misc/formio.h>

#include <math/isosurf/shape.h>
#include <chemistry/molecule/atominfo.h>
#include <chemistry/molecule/molecule.h>

namespace sc {

/** The VDWShape class describes the surface of a
    molecule as the union of atom centered spheres, each the
    van der Waals radius of the atom.
*/
class VDWShape: public UnionShape {
 private:
    Ref<AtomInfo> atominfo_;
 public:
    VDWShape(const Ref<Molecule>&);
    VDWShape(const Ref<KeyVal>&);
    ~VDWShape();
    void initialize(const Ref<Molecule>&);
};  

/** DiscreteConnollyShape and ConnollyShape should produce the same result.
    The discrete version is a shape union of discrete subshapes and is
    slower.  These classes describe the solvent accessible surface of a
    molecule.  */
class DiscreteConnollyShape: public UnionShape {
  private:
    double radius_scale_factor_;
    Ref<AtomInfo> atominfo_;
 public:
    DiscreteConnollyShape(const Ref<KeyVal>&);
    ~DiscreteConnollyShape();
    void initialize(const Ref<Molecule>&,double probe_radius);
};

#ifndef COUNT_CONNOLLY
# define COUNT_CONNOLLY 1
#endif

// This is a utility class needed by ConnollyShape2
class CS2Sphere
{
    SCVector3 _v;
    double _radius;

  public:
#if COUNT_CONNOLLY
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
    void print(std::ostream& os=ExEnv::out0()) const
    {
      os << indent
         << scprintf("Rad=%lf, Center=(%lf,%lf,%lf), From origin=%lf\n",
                     _radius, _v[0], _v[1], _v[2], _v.norm());
    }

    // Function to determine if there is any portion of this that 
    // is not inside one or more of the spheres in s[].  Returns
    // 1 if the intersection is empty, otherwise 0 is returned.
    // Warning: the spheres in s are modified.
    int intersect(CS2Sphere *s,
                  int n_spheres) const;

    static void print_counts(std::ostream& = ExEnv::out0());
};

#define CONNOLLYSHAPE_N_WITH_NSPHERE_DIM 10
/** DiscreteConnollyShape and ConnollyShape should produce the same result.
    The discrete version is a shape union of discrete subshapes and is
    slower.  These classes describe the solvent accessible surface of a
    molecule. */
class ConnollyShape: public Shape {
  private:
    CS2Sphere* sphere;
    double probe_r;
    double radius_scale_factor_;
    int n_spheres;
    Ref<AtomInfo> atominfo_;

    std::vector<int> ***box_;
    double l_;
    int xmax_;
    int ymax_;
    int zmax_;
    SCVector3 lower_;

    int get_box(const SCVector3 &v, int &x, int &y, int &z) const;

#if COUNT_CONNOLLY
    static int n_total_;
    static int n_inside_vdw_;
    static int n_with_nsphere_[CONNOLLYSHAPE_N_WITH_NSPHERE_DIM];
#endif

  public:
    ConnollyShape(const Ref<KeyVal>&);
    ~ConnollyShape();
    void initialize(const Ref<Molecule>&,double probe_radius);
    void clear();
    double distance_to_surface(const SCVector3&r,
                               SCVector3*grad=0) const;
    void boundingbox(double valuemin,
                     double valuemax,
                     SCVector3& p1, SCVector3& p2);

    static void print_counts(std::ostream& = ExEnv::out0());
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
