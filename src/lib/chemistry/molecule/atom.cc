//
// atom.h
//
// Copyright (C) 2013 Drew Lewis
//
// Author: Drew Lewis <drew90@vt.edu>
// Maintainer: Drew Lewis and Edward Valeev
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

#include <chemistry/molecule/atom.h>

void
sc::ToStateOut(const Atom &a, StateOut &so, int &count) {
    count += so.put_array_double(a.r(), 3);
    count += so.put(a.Z());
    count += so.put(a.have_charge());
    count += so.put(a.have_fragment());
    count += so.put(a.charge());
    count += so.put(a.fragment());
    count += so.put(a.mass());
    count += so.put(a.label());
}

void
sc::FromStateIn(Atom &a, StateIn &si, int &count){
    count += si.get_array_double(a.r_,3);
    count += si.get(a.Z_);
    count += si.get(a.have_charge_);
    count += si.get(a.have_fragment_);
    count += si.get(a.charge_);
    count += si.get(a.fragment_);
    count += si.get(a.mass_);
    count += si.get(a.label_);
}

bool sc::operator ==(const Atom& a, const Atom& b) {
  if (a.Z() != b.Z())
    return false;
  for(int xyz=0; xyz<3; ++xyz)
    if (a.r(xyz) != b.r(xyz) )
      return false;
  if (a.have_charge() != b.have_charge())
    return false;
  if (a.have_charge()) {
    if (a.charge() != b.charge())
      return false;
  }
  if (a.have_fragment() != b.have_fragment())
    return false;
  if (a.have_fragment()) {
    if (a.fragment() != b.fragment())
      return false;
  }
  if (a.mass() != b.mass())
    return false;
  // labels are inconsequential
  return true;
}
