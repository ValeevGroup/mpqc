//
// formula.cc --- implementation of the MolecularFormula class
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <sstream>
#include <map>
#include <chemistry/molecule/formula.h>

using namespace sc;

MolecularFormula::MolecularFormula(const Ref<Molecule>&m):
  form_(0)
{
  compute_form(m.pointer());
  compute_atomtypes(m.pointer());
}

MolecularFormula::MolecularFormula(const Molecule *m):
  form_(0)
{
  compute_form(m);
  compute_atomtypes(m);
}

MolecularFormula::~MolecularFormula()
{
  delete[] form_;
  delete[] Z_;
  delete[] nZ_;
}

void
MolecularFormula::compute_form(const Molecule *m)
{
  const Molecule& mol = *m;

  std::map<std::string, int> count;

  for (int a=0; a < mol.natom(); a++) {
    std::string symbol(mol.atom_symbol(a));
    if (count.find(symbol) == count.end()) count[symbol] = 0;
    count[symbol]++;
  }

  std::ostringstream sstr;
  if (count.find("Q") != count.end()) {
      sstr << "Q";
      if (count["Q"] > 1) sstr << count["Q"];
      count.erase("Q");
    }
  if (count.find("C") != count.end()) {
      sstr << "C";
      if (count["C"] > 1) sstr << count["C"];
      count.erase("C");
    }
  if (count.find("H") != count.end()) {
      sstr << "H";
      if (count["H"] > 1) sstr << count["H"];
      count.erase("H");
    }
  for (std::map<std::string,int>::iterator i = count.begin();
       i != count.end(); i++) {
      sstr << i->first;
      if (i->second > 1) sstr << i->second;
    }
  form_ = strcpy(new char[sstr.str().size()+1],sstr.str().c_str());
}

void
MolecularFormula::compute_atomtypes(const Molecule *m)
{
  std::map<int, int> atomtypeinfo;
  int natoms = m->natom();
  int i, Z;

  for (i=0; i< natoms; i++) {
    // don't include ghost functions or point charges
    if (m->charge(i) == 0.0) continue;
    if (m->atom_symbol(i) == "Q") continue;
    Z = m->Z(i);
    if (atomtypeinfo.find(Z) != atomtypeinfo.end()) atomtypeinfo[Z]++;
    else atomtypeinfo[Z] = 1;
    }

  natomtypes_ = atomtypeinfo.size();

  Z_ = new int[natomtypes_];
  nZ_ = new int[natomtypes_];

  std::map<int, int>::iterator iter;

  for (iter = atomtypeinfo.begin(), i=0; iter != atomtypeinfo.end(); iter++, i++) {
    Z_[i] = iter->first;
    nZ_[i] = iter->second;
    }

}

const char *
MolecularFormula::formula() const
{
  return form_;
}

int
MolecularFormula::natomtypes() 
{
  return(natomtypes_);
}

int
MolecularFormula::Z(int itype) 
{
  return Z_[itype];
}

int
MolecularFormula::nZ(int itype) 
{
  return nZ_[itype];
}
