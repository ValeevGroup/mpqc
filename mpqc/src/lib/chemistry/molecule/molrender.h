
#ifndef _chemistry_molecule_molrender_h
#define _chemistry_molecule_molrender_h

#include <util/render/object.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/molinfo.h>
#include <math/isosurf/surf.h>

class RenderedMolecule: public RenderedObject {
#   define CLASSNAME RenderedMolecule
#   include <util/class/classda.h>
  protected:
    RefRenderedObject object_;
    RefMolecule mol_;
    RefAtomInfo atominfo_;

    virtual void init() = 0;
  public:
    RenderedMolecule(const RefKeyVal& keyval);
    ~RenderedMolecule();

    void render(const RefRender&);
};

class RenderedStickMolecule: public RenderedMolecule {
#   define CLASSNAME RenderedStickMolecule
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    void init();
  public:
    RenderedStickMolecule(const RefKeyVal& keyval);
    ~RenderedStickMolecule();
};

class RenderedBallMolecule: public RenderedMolecule {
#   define CLASSNAME RenderedBallMolecule
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    void init();
  public:
    RenderedBallMolecule(const RefKeyVal& keyval);
    ~RenderedBallMolecule();
};

class RenderedMolecularSurface: public RenderedMolecule {
#   define CLASSNAME RenderedMolecularSurface
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    RefTriangulatedImplicitSurface surf_;
    void init();
  public:
    RenderedMolecularSurface(const RefKeyVal& keyval);
    ~RenderedMolecularSurface();
};

#endif
