//
// trange1_test.cpp
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

#define BOOST_TEST_MODULE TRange1 generators

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/split.h>
#include <vector>
#include <utility>
#include "trange1.hpp"

BOOST_AUTO_TEST_CASE(tile_by_atom) {
    using namespace mpqc;
    namespace TA = TiledArray;

    // C1 Molecule with a different atom in every direction
    sc::Ref<sc::Molecule>  mol = new sc::Molecule;
    mol->add_atom( 6,     0,     0,    0);
    mol->add_atom( 9,    -1,    -1,    0);
    mol->add_atom( 1,   0.6,  -0.1,  0.9);
    mol->add_atom(17, -0.75,   1.5,    0);
    mol->add_atom(35,   1.1, -0.18, -1.5);

    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("name", "3-21G");
    akv->assign("molecule", mol.pointer());

    sc::Ref<sc::SplitBasisSet> basis =
                               new sc::SplitBasisSet(sc::Ref<sc::KeyVal>(akv));

    BOOST_TEST_MESSAGE("Testing MPQC TiledArray tiling by atom.");
    {
        TA::TiledRange1 tr1 = tiling::tile_by_atom(basis);
        std::vector<std::pair<std::size_t, std::size_t> > ranges;
        ranges.push_back(std::make_pair(0,9));
        ranges.push_back(std::make_pair(9,18));
        ranges.push_back(std::make_pair(18,20));
        ranges.push_back(std::make_pair(20,33));
        ranges.push_back(std::make_pair(33,56));

        for(auto i = 0; i < 5; ++i){
            BOOST_CHECK_EQUAL(tr1.tile(i).first, ranges[i].first);
            BOOST_CHECK_EQUAL(tr1.tile(i).second, ranges[i].second);
        }
    }
}

/*
BOOST_AUTO_TEST_CASE(tile_by_shell) {
    using namespace mpqc;
    namespace TA = TiledArray;

    // C1 Molecule with a different atom in every direction
    sc::Ref<sc::Molecule>  mol = new sc::Molecule;
    mol->add_atom( 6,     0,     0,    0);
    mol->add_atom( 9,    -1,    -1,    0);
    mol->add_atom( 1,   0.6,  -0.1,  0.9);
    mol->add_atom(17, -0.75,   1.5,    0);
    mol->add_atom(35,   1.1, -0.18, -1.5);

    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("name", "3-21G");
    akv->assign("molecule", mol.pointer());

    sc::Ref<sc::SplitBasisSet> basis =
                               new sc::SplitBasisSet(sc::Ref<sc::KeyVal>(akv));

    BOOST_TEST_MESSAGE("Testing MPQC TiledArray tiling by atom.");
    {
        TA::TiledRange1 tr1 = tiling::tile_by_shell(basis);
        std::vector<std::pair<std::size_t, std::size_t> > ranges;

}
*/





