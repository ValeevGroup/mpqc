//
// symmetrize.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <iostream.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>

int
main(int argc, char *argv[])
{
  int i;

  char *infile = (argv[1]) ? argv[1] : "mpqc.in";
  RefKeyVal kv(new ParsedKeyVal(infile));

  RefMolecule mol = kv->describedclassvalue("molecule");

  mol->print();

  mol->move_to_com();
  cout << "Molecule at com:\n";
  mol->print();
  
  mol->transform_to_principal_axes(0);
  cout << "Molecule wrt principal axes:\n";
  mol->print();
  mol->point_group()->symm_frame().print();

  mol->symmetrize();
  cout << "symmetrized molecule\n";
  mol->print();

  mol->cleanup_molecule();
  cout << "cleaned molecule\n";
  mol->print();
  
  int nunique = mol->num_unique_atoms();
  int * unique_atoms = mol->find_unique_atoms();

  cout << scprintf("\nnunique=%d: ",nunique);
  for (i=0; i < nunique; i++) cout << scprintf(" %d",unique_atoms[i]+1);
  cout << endl;

  RefMolecule unique = new Molecule;
  for (i=0; i < nunique; i++)
    unique->add_atom(mol->Z(unique_atoms[i]),
                     mol->r(unique_atoms[i])[0],
                     mol->r(unique_atoms[i])[1],
                     mol->r(unique_atoms[i])[2]
                     );

  cout << "unique atoms\n";
  unique->print();
  
  exit(0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
