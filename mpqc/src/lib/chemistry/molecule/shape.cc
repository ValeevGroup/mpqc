
extern "C" {
# include <math.h>
  }

#include "shape.h"
#include "molecule.h"

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
  ArraysetRefSphereShape spheres;
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
