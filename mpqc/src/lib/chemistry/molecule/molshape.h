
#ifndef _chemistry_molecule_molshape_h
#define _chemistry_molecule_molshape_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/shape.h>
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
 public:
    ConnollyShape(const RefMolecule&,double probe_radius = 2.6456173);
    ConnollyShape(const RefKeyVal&);
    ~ConnollyShape();
    void initialize(const RefMolecule&,double probe_radius);
};  

#endif
