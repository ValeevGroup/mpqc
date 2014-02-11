//
// tascf_test.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <chemistry/qc/scf/tascf.hpp>
#define BOOST_TEST_MODULE test_tascf
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;
using namespace sc;
using namespace mpqc;

BOOST_AUTO_TEST_CASE( construct_scf_programmatically ){

    // Make a molecule H2
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(1,0,1,-1);
    mol->add_atom(1,0,1,1);

    // Make keyval
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", mol.pointer());
    Ref<GaussianBasisSet> basis =
                    new GaussianBasisSet(Ref<KeyVal>(akv));
    akv->assign("basis", basis.pointer());

    //Construct object
    Ref<v3::SCF> tscf = new v3::SCF(akv);
    tscf->print();
}

BOOST_AUTO_TEST_CASE( construct_scf_txtkeyval ){

  const char *input =      SRCDIR "/tascf_test.in";
  Ref<KeyVal> kv = new ParsedKeyVal(input);
  Ref<v3::SCF> tscf; tscf << kv->describedclassvalue("rhf");

  tscf->print();

}

