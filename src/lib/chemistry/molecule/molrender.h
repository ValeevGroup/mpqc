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

#ifndef _chemistry_molecule_molrender_h
#define _chemistry_molecule_molrender_h

#include <util/render/object.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/atominfo.h>
#include <math/isosurf/surf.h>

namespace sc {

class RenderedMolecule: public RenderedObject {
  protected:
    Ref<RenderedObject> object_;
    Ref<Molecule> mol_;
    Ref<AtomInfo> atominfo_;

  public:
    RenderedMolecule(const Ref<KeyVal>& keyval);
    ~RenderedMolecule();

    Ref<Molecule> molecule() { return mol_; }

    // init must be called if the molecule changes
    virtual void init() = 0;

    void render(const Ref<Render>&);
};


class RenderedStickMolecule: public RenderedMolecule {
  protected:
    int use_color_;
  public:
    RenderedStickMolecule(const Ref<KeyVal>& keyval);
    ~RenderedStickMolecule();

    void init();
};

class RenderedBallMolecule: public RenderedMolecule {
  protected:
  public:
    RenderedBallMolecule(const Ref<KeyVal>& keyval);
    ~RenderedBallMolecule();

    void init();
};

class MoleculeColorizer: public DescribedClass {
  protected:
    Ref<Molecule> mol_;
  public:
    MoleculeColorizer(const Ref<Molecule> &);
    MoleculeColorizer(const Ref<KeyVal>&);
    ~MoleculeColorizer();

    virtual void colorize(const Ref<RenderedPolygons> &) = 0;
};


class AtomProximityColorizer: public MoleculeColorizer {
  protected:
    Ref<AtomInfo> atominfo_;
  public:
    AtomProximityColorizer(const Ref<Molecule>&, const Ref<AtomInfo> &);
    AtomProximityColorizer(const Ref<KeyVal> &);
    ~AtomProximityColorizer();

    void colorize(const Ref<RenderedPolygons> &);
};

class RenderedMolecularSurface: public RenderedMolecule {
  protected:
    Ref<TriangulatedImplicitSurface> surf_;
    Ref<MoleculeColorizer> colorizer_;
  public:
    RenderedMolecularSurface(const Ref<KeyVal>& keyval);
    ~RenderedMolecularSurface();

    void init(int reinit_surf);
    void init();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
