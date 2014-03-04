//
// tiledbasisset_test.cpp
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

#include <util/madness/init.h>
#include <util/misc/regtime.h>
#include <util/madness/world.h>
#include <iostream>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/qc/basis/integral.h>
#include <mpqc/integrals/integralenginepool.hpp>
#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#define BOOST_TEST_MODULE test_tabasis
#include <boost/test/included/unit_test.hpp>



using namespace boost::unit_test;
using namespace sc;
using namespace mpqc;

struct MADConfig {
    MADConfig() {
      mpqc::MADNESSRuntime::initialize();
    }
    ~MADConfig() {
      mpqc::MADNESSRuntime::finalize();
    }
};

struct MolFixture {
    MolFixture(){
        mol = new Molecule;
        mol->add_atom(2, 0.0, 0.0, 0.0);
        mol->add_atom(2, 0.0, 0.0, 3.7);
        mol->add_atom(10, 2.8, 0.0, 1.8);
        mol->add_atom(18, 0.0, 2.8, 1.8);
        mol->add_atom(36, 0.0, -2.8, 1.8);
        mol->add_atom(18, -2.8, 0.0, 1.8);
    }
    Ref<Molecule> mol;
};

BOOST_GLOBAL_FIXTURE( MADConfig );

BOOST_AUTO_TEST_CASE( tiledbasisset_constructor_test ){
    MolFixture MF;

    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", MF.mol.pointer());
    Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

    // Construct via gbs
    for(size_t i = 1; i <= 6; ++i){
        Ref<TA::TiledBasisSet> tbs_from_gbs = new TA::TiledBasisSet(bs,i);
        BOOST_REQUIRE(!tbs_from_gbs.null());
        BOOST_CHECK(tbs_from_gbs->trange1().tiles().second == i);
    }
    // Construct via kval
    for(size_t i = 1; i <= 6; ++i){
        akv->assign("ntiles", int(i));
        Ref<TA::TiledBasisSet> tbs_from_kval = new TA::TiledBasisSet(
                                                Ref<KeyVal>(akv));
        BOOST_REQUIRE(!tbs_from_kval.null());
        BOOST_CHECK(tbs_from_kval->trange1().tiles().second == i);
    }
}

BOOST_AUTO_TEST_CASE( tiledbasisset_trange1_values_test ) {

    MolFixture MF;
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", MF.mol.pointer());
    Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

     // Test Single tile case
    Ref<TA::TiledBasisSet> tbs1 = new TA::TiledBasisSet(bs,1);
    double first = tbs1->trange1().tile(0).first;
    double second = tbs1->trange1().tile(0).second;
    BOOST_CHECK(first == 0);
    BOOST_CHECK(second == 43);

     // Test 3 tile case
    Ref<TA::TiledBasisSet> tbs3 = new TA::TiledBasisSet(bs,3);
    double t1_first = tbs3->trange1().tile(0).first;
    double t1_second = tbs3->trange1().tile(0).second;
    double t2_second = tbs3->trange1().tile(1).second;
    double t3_second = tbs3->trange1().tile(2).second;
    BOOST_CHECK(t1_first == 0);
    BOOST_CHECK(t1_second == 37);
    BOOST_CHECK(t2_second == 38);
    BOOST_CHECK(t3_second == 43);
}

BOOST_AUTO_TEST_CASE( tiledbasisset_deterministic_tiliing_test ){
    std::array<std::array<double, 3>, 6> c{{
                                             {{0.0,0.0,0.0}},
                                             {{0.0,0.0,3.7}},
                                             {{2.8,0.0,1.8}},
                                             {{0.0,2.8,1.8}},
                                             {{0.0,-2.8,1.8}},
                                             {{-2.8,0.0,1.8}}
                                          }};
    for(size_t i = 0; i < c.size(); ++i){
        // Scramble molecule input
        for(size_t j = 0; j < c.size(); ++j){
            c[i].swap(c[j]);

            Ref<Molecule> mol = new Molecule;
            mol->add_atom(2,  c[0][0], c[0][1], c[0][2]);
            mol->add_atom(2,  c[1][0], c[1][1], c[1][2]);
            mol->add_atom(10, c[2][0], c[2][1], c[2][2]);
            mol->add_atom(18, c[3][0], c[3][1], c[3][2]);
            mol->add_atom(36, c[4][0], c[4][1], c[4][2]);
            mol->add_atom(18, c[5][0], c[5][1], c[5][2]);

            Ref<AssignedKeyVal> akv = new AssignedKeyVal;
            akv->assign("name", "STO-3G");
            akv->assign("molecule", mol.pointer());
            Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

             // Test 3 tile case and make sure tiling is the same.
            Ref<TA::TiledBasisSet> tbs3 = new TA::TiledBasisSet(bs,3);
            double t1_second = tbs3->trange1().tile(0).second;
            double t2_second = tbs3->trange1().tile(1).second;
            BOOST_CHECK_EQUAL(t1_second, 37);
            BOOST_CHECK_EQUAL(t2_second, 38);

            c[j].swap(c[i]);
        }
    }
}











