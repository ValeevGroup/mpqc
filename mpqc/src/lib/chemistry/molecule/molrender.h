//
// molrender.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_molecule_molrender_h
#define _chemistry_molecule_molrender_h

#include <util/render/object.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/atominfo.h>
#include <math/isosurf/surf.h>

class RenderedMolecule: public RenderedObject {
#   define CLASSNAME RenderedMolecule
#   include <util/class/classda.h>
  protected:
    RefRenderedObject object_;
    RefMolecule mol_;
    RefAtomInfo atominfo_;

  public:
    RenderedMolecule(const RefKeyVal& keyval);
    ~RenderedMolecule();

    RefMolecule molecule() { return mol_; }

    // init must be called if the molecule changes
    virtual void init() = 0;

    void render(const RefRender&);
};
DescribedClass_REF_dec(RenderedMolecule);

class RenderedStickMolecule: public RenderedMolecule {
#   define CLASSNAME RenderedStickMolecule
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    int use_color_;
  public:
    RenderedStickMolecule(const RefKeyVal& keyval);
    ~RenderedStickMolecule();

    void init();
};

class RenderedBallMolecule: public RenderedMolecule {
#   define CLASSNAME RenderedBallMolecule
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
  public:
    RenderedBallMolecule(const RefKeyVal& keyval);
    ~RenderedBallMolecule();

    void init();
};

class MoleculeColorizer: public DescribedClass {
#   define CLASSNAME MoleculeColorizer
#   include <util/class/classda.h>
  protected:
    RefMolecule mol_;
  public:
    MoleculeColorizer(const RefMolecule &);
    MoleculeColorizer(const RefKeyVal&);
    ~MoleculeColorizer();

    virtual void colorize(const RefRenderedPolygons &) = 0;
};
DescribedClass_REF_dec(MoleculeColorizer);

class AtomProximityColorizer: public MoleculeColorizer {
#   define CLASSNAME AtomProximityColorizer
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    RefAtomInfo atominfo_;
  public:
    AtomProximityColorizer(const RefMolecule&, const RefAtomInfo &);
    AtomProximityColorizer(const RefKeyVal &);
    ~AtomProximityColorizer();

    void colorize(const RefRenderedPolygons &);
};

class RenderedMolecularSurface: public RenderedMolecule {
#   define CLASSNAME RenderedMolecularSurface
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    RefTriangulatedImplicitSurface surf_;
    RefMoleculeColorizer colorizer_;
  public:
    RenderedMolecularSurface(const RefKeyVal& keyval);
    ~RenderedMolecularSurface();

    void init();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
