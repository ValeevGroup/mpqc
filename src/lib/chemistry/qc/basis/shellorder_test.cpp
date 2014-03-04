//
// shellorder_test.cpp
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


#include <chemistry/qc/basis/shellorder.hpp>

#define BOOST_TEST_MODULE test_shellorder
#include <boost/test/included/unit_test.hpp>

using namespace mpqc;
using namespace sc;
using namespace boost::unit_test;

// Ensure the constructor works
BOOST_AUTO_TEST_CASE(shellorder_constructor_test){
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", mol.pointer());
    Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

    TA::ShellOrder sh(bs);

    BOOST_REQUIRE(&sh != nullptr);
}

// Test that the shell range function returns something reasonable.
BOOST_AUTO_TEST_CASE( shellorder_shellrange_test ){
    using ShellR = mpqc::TA::ShellOrder::ShellRange;
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", mol.pointer());
    Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

    std::vector<ShellR> correct_ranges;
    correct_ranges.push_back(ShellR{0,15}); // one tile
    correct_ranges.push_back(ShellR{0,14,15}); // two tiles
    correct_ranges.push_back(ShellR{0,9,10,15}); // three tiles
    correct_ranges.push_back(ShellR{0,6,7,12,15}); // four tiles
    correct_ranges.push_back(ShellR{0,4,5,10,13,15}); // five tiles
    correct_ranges.push_back(ShellR{0,1,2,7,10,12,15}); // six tiles


    for(size_t i = 1; i <= 6; ++i){
        TA::ShellOrder sh(bs);

        // Not sure how to test this yet, but make clusters
        std::vector<mpqc::TA::ShellOrder::Shell> shells = sh.ordered_shells(i);

        ShellR ranges = sh.shell_ranges();
        BOOST_CHECK_EQUAL_COLLECTIONS(ranges.begin(),ranges.end(),
                      correct_ranges[i-1].begin(), correct_ranges[i-1].end());
    }

}
