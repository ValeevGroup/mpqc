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

#include <util/container/avlmap.h>
#include <chemistry/molecule/formula.h>

static const char * symbols[] = {
  "C", "H", "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "B", "Ba", "Be",
  "Bi", "Bk", "Br", "Ca", "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu",
  "Dy", "Er", "Es", "Eu", "F", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "Ha", "He",
  "Hf", "Hg", "Ho", "I", "In", "Ir", "K", "Kr", "La", "Li", "Lr", "Lu", "Md",
  "Mg", "Mn", "Mo", "N", "Na", "Nb", "Nd", "Ne", "Ni", "No", "Np", "O", "Os",
  "P", "Pa", "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rf",
  "Rh", "Rn", "Ru", "S", "Sb", "Sc", "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb",
  "Tc", "Te", "Th", "Ti", "Tl", "Tm", "U", "V", "W", "Xe", "Y", "Yb", "Zn",
  "Zr", 0
};

MolecularFormula::MolecularFormula(const RefMolecule&m):
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

  memset(count_, 0, sizeof(count_));

  int ntype=0;
  int maxcount=0;
  unsigned int maxsym=0;
  for (int a=0; a < mol.natom(); a++) {
    int i=0;
    while(symbols[i]) {
      if (!strcmp(AtomInfo::symbol(mol.Z(a)), symbols[i])) {
        count_[i]++;

        maxcount = (count_[i] > maxcount) ? count_[i] : maxcount;
        maxsym = (strlen(symbols[i]) > maxsym) ? strlen(symbols[i]) : maxsym;
        
        if (count_[i]==1)
          ntype++;
        
        break;
      }
      i++;
    }
  }

  // allocate storage for formula
  int ndigits = ((int) (log((double)maxcount)/log(10.0))) + 1;
  form_ = new char[(ndigits+maxsym)*ntype+1];
  form_[0] = 0;
  
  int c;
  for (int i=0; i < nelem_; i++) {
    if ((c=count_[i])) {
      char *temp = new char[ndigits+maxsym+1];
      if (c > 1)
        sprintf(temp, "%s%d", symbols[i], count_[i]);
      else 
        sprintf(temp, "%s", symbols[i]);

      strcat(form_, temp);
      delete[] temp;
    }
  }
}

void
MolecularFormula::compute_atomtypes(const Molecule *m)
{
  AVLMap<int, int> atomtypeinfo;
  int natoms = m->natom();
  int i, Z;

  for (i=0; i< natoms; i++) {
    if (m->charge(i) == 0.0) continue;
    Z = m->Z(i);
    if (atomtypeinfo.find(Z) != atomtypeinfo.end()) atomtypeinfo[Z]++;
    else atomtypeinfo[Z] = 1;
    }

  natomtypes_ = atomtypeinfo.length();

  Z_ = new int[natomtypes_];
  nZ_ = new int[natomtypes_];

  AVLMap<int, int>::iterator iter;

  for (iter = atomtypeinfo.begin(), i=0; iter != atomtypeinfo.end(); iter++, i++) {
    Z_[i] = iter.key();
    nZ_[i] = iter.data();
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
