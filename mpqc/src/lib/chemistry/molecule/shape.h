
#ifndef _chemistry_molecule_shape_h
#define _chemistry_molecule_shape_h

#include <math/isosurf/shape.h>
#include <chemistry/molecule/molecule.h>

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
