//
// formula.h --- class for calculation molecular formulae
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

#ifndef _chemistry_molecule_formula_h
#define _chemistry_molecule_formula_h

#include <chemistry/molecule/molecule.h>

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

/** The MolecularFormula class is used to calculate the molecular
 formula of a Molecule.  There is only one constructor which
 takes Ref<Molecule> as input. */
class MolecularFormula {
  private:
    int natomtypes_;
    int *Z_, *nZ_;
    char *form_;

    void compute_atomtypes(const Molecule *m);
    void compute_form(const Molecule *m);
  public:
    /// Constructors.  The argument must be nonnull.
    MolecularFormula(const Ref<Molecule>&m);
    MolecularFormula(const Molecule *m);

    ~MolecularFormula();

    /// Returns a null terminated string containing the molecular formula.
    const char * formula() const;
    /// Returns the number of atomtypes
    int natomtypes(); 
    /// Returns atomic number of given atomtypeindex
    int Z(int itype); 
    /// Returns number of atoms of given atomtypeindex
    int nZ(int itype);
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif
