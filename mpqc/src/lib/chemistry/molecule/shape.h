
#ifndef _chemistry_molecule_shape_h
#define _chemistry_molecule_shape_h

#include <math/isosurf/shape.h>
#include <chemistry/molecule/molecule.h>

class VDWShape: public UnionShape {
#   define CLASSNAME VDWShape
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 public:
    VDWShape(Molecule&);
    VDWShape(KeyVal&);
    ~VDWShape();
    void initialize(Molecule&);
};  

class ConnollyShape: public UnionShape {
#   define CLASSNAME ConnollyShape
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 public:
    ConnollyShape(Molecule&,double probe_radius = 2.6456173);
    ConnollyShape(KeyVal&);
    ~ConnollyShape();
    void initialize(Molecule&,double probe_radius);
};  

#endif
