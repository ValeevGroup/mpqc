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

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
