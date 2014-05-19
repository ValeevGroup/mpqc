//
// tawfn_test.cpp
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

#include <iostream>
#include "tawfn.hpp"
#define BOOST_TEST_MODULE test_mock_wfn
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;
using namespace mpqc;

class Mock_wfn : public Wavefunction {
public:
    Mock_wfn(sc::Ref<sc::KeyVal> &kval) : Wavefunction(kval) {}
    virtual ~Mock_wfn() {};
    int nelectron() override { return 2.0; }
    void compute() override {}
};
BOOST_AUTO_TEST_CASE( test_mock_wfn_constructor ){

    BOOST_MESSAGE("Testing TiledArrayWavefunction");
    // Make a molecule H2
    sc::Ref<sc::Molecule> mol = new sc::Molecule;
    mol->add_atom(1,0,1,-1);
    mol->add_atom(1,0,1,1);

    // Make keyval
    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", mol.pointer());
    sc::Ref<sc::GaussianBasisSet> basis =
                    new sc::GaussianBasisSet(sc::Ref<sc::KeyVal>(akv));
    akv->assign("basis", basis.pointer());
    sc::Ref<sc::KeyVal> kval = sc::Ref<sc::KeyVal>(akv);
    Mock_wfn twfn(kval);

}
