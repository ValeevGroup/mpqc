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
#include <utility>

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

// Test that the shell range function doesn't depend on the ordering by trying
// all 720 permuations.
BOOST_AUTO_TEST_CASE( shellorder_shellrange_atom_order_test ){
    using ShellR = mpqc::TA::ShellOrder::ShellRange;
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    std::vector<sc::Atom> atoms = mol->atoms();

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

    for(auto i = 0; i < 6; ++i){
        std::swap(atoms[0], atoms[i]);
    for(auto j = 1; j < 6; ++j){
        std::swap(atoms[1], atoms[j]);
    for(auto k = 2; k < 6; ++k){
        std::swap(atoms[2], atoms[k]);
    for(auto l = 3; l < 6; ++l){
        std::swap(atoms[3], atoms[l]);
    for(auto m = 4; m < 6; ++m){
        std::swap(atoms[4], atoms[m]);
       Ref<Molecule> mol_t = new Molecule;
       mol_t->add_atom(atoms[0].Z(),atoms[0].r(0),atoms[0].r(1),atoms[0].r(2));
       mol_t->add_atom(atoms[1].Z(),atoms[1].r(0),atoms[1].r(1),atoms[1].r(2));
       mol_t->add_atom(atoms[2].Z(),atoms[2].r(0),atoms[2].r(1),atoms[2].r(2));
       mol_t->add_atom(atoms[3].Z(),atoms[3].r(0),atoms[3].r(1),atoms[3].r(2));
       mol_t->add_atom(atoms[4].Z(),atoms[4].r(0),atoms[4].r(1),atoms[4].r(2));
       mol_t->add_atom(atoms[5].Z(),atoms[5].r(0),atoms[5].r(1),atoms[5].r(2));

       Ref<AssignedKeyVal> akv = new AssignedKeyVal;
       akv->assign("name", "STO-3G");
       akv->assign("molecule", mol.pointer());
       Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

        for(size_t q = 1; q <= 6; ++q){
            TA::ShellOrder sh(bs);

            // Not sure how to test this yet, but make clusters
            std::vector<mpqc::TA::ShellOrder::Shell> shells =
                                                        sh.ordered_shells(i);

            ShellR ranges = sh.shell_ranges();
            BOOST_CHECK_EQUAL_COLLECTIONS(ranges.begin(),ranges.end(),
                      correct_ranges[i-1].begin(), correct_ranges[i-1].end());
        }

        std::swap(atoms[m], atoms[4]);
    }
        std::swap(atoms[l], atoms[3]);
    }
        std::swap(atoms[k], atoms[2]);
    }
        std::swap(atoms[j], atoms[1]);
    }
        std::swap(atoms[i], atoms[0]);
    }
}
