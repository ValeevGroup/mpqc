//
// atomcent.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _chemistry_molecule_atomcent_h
#define _chemistry_molecule_atomcent_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <chemistry/molecule/chemelem.h>
#include <math/topology/point.h>

//.  The \clsnm{AtomicCenter} class is used to describe atoms as points in
//space.  Thus an \clsnm{AtomicCenter} contains a
//\clsnmref{ChemicalElement} and a \clsnm{Point}.
//
// \clsnm{AtomicCenter} is a SavableState and has a \clsnmref{StateIn}
// constructor.  \clsnm{AtomicCenter} does not have a \clsnmref{KeyVal}
// constructor.
class AtomicCenter: public SavableState
{
#   define CLASSNAME AtomicCenter
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    Point p;
    ChemicalElement element_;
    char *label_;
  public:
    AtomicCenter();
    AtomicCenter(const AtomicCenter&);
    AtomicCenter(StateIn&);
    //. The first argument to this constructor is a string suitable for use
    //in the constructor of a \clsnmref{ChemicalElement} (@ref{The
    //ChemicalElement Class}).  The next three arguments are the x, y, and
    //z Cartesian coordinates of the atom.  The last optional argument is a
    //label for the atom.
    AtomicCenter(const char*symbol,double x,double y,double z, const char* =0);

    ~AtomicCenter();

    AtomicCenter& operator=(const AtomicCenter&ac);

    //. Returns a reference to the i'th cartesian coordinate
    double& operator[](int i) { return p[i]; };
    //. Casts \clsnm{AtomicCenter} to an \srccd{int} which is the atomic
    // number.
    operator int() { return element_.number(); };
    //. Returns the \clsnmref{ChemicalElement} for this atom.
    ChemicalElement& element() { return element_; };
    //. \srccd{const} version of the above.
    const ChemicalElement& element() const { return element_; };
    //. Returns the \clsnm{Point} containing the cartesian coordinates.
    Point& point() { return p; };
    //. \srccd{const} version of the above.
    const Point& point() const { return p; };
    //. Returns either the x, y, or z coordinate.
    const double& r(int xyz) const { return p[xyz]; };
    //. Returns the label for this atom.
    const char * label() const { return label_; }

    //. Returns the mass for this atom.
    double mass() const;

    void save_data_state(StateOut&);

    //. Print information about the atom
    void print(ostream& =cout);
};
DescribedClass_REF_dec(AtomicCenter);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
