
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
# include <math.h>
  }

#include "molshape.h"
#include "molecule.h"

#define CLASSNAME VDWShape
#define PARENTS public UnionShape
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VDWShape::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = UnionShape::_castdown(cd);
  return do_castdowns(casts,cd);
}

VDWShape::VDWShape(const RefMolecule&mol)
{
  initialize(mol);
}

VDWShape::VDWShape(const RefKeyVal&keyval)
{
  RefMolecule mol = keyval->describedclassvalue("molecule");
  initialize(mol);
}

void
VDWShape::initialize(const RefMolecule&mol)
{
  _shapes.clear();
  for (int i=0; i<mol->natom(); i++) {
      SCVector3 r;
      for (int j=0; j<3; j++) r[j] = mol->operator[](i)[j];
      add_shape(
          new SphereShape(r,mol->operator[](i).element().atomic_radius())
          );
    }
}

VDWShape::~VDWShape()
{
}

#define CLASSNAME ConnollyShape
#define PARENTS public UnionShape
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ConnollyShape::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = UnionShape::_castdown(cd);
  return do_castdowns(casts,cd);
}

ConnollyShape::ConnollyShape(const RefMolecule&mol,double probe_radius)
{
  initialize(mol,probe_radius);
}

ConnollyShape::ConnollyShape(const RefKeyVal&keyval)
{
  RefMolecule mol = keyval->describedclassvalue("molecule");
  double probe_radius = keyval->doublevalue("probe_radius");
  if (keyval->error() != KeyVal::OK) {
      probe_radius = 2.6456173;
    }
  initialize(mol,probe_radius);
}

void
ConnollyShape::initialize(const RefMolecule&mol,double probe_radius)
{
  _shapes.clear();
  ArraysetRefSphereShape spheres;
  for (int i=0; i<mol->natom(); i++) {
      SCVector3 r;
      for (int j=0; j<3; j++) r[j] = mol->operator[](i)[j];
      RefSphereShape
        sphere(
            new SphereShape(r,mol->operator[](i).element().atomic_radius())
            );
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
