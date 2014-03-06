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
#include <utility>



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

BOOST_GLOBAL_FIXTURE( MADConfig );

BOOST_AUTO_TEST_CASE( tiledbasisset_gbs_constructor_test ){

    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    std::vector<sc::Atom> atoms = mol->atoms();

    for(size_t i = 0; i < 6; ++i){
        std::swap(atoms[0], atoms[i]);
    for(size_t j = 1; j < 6; ++j){
        std::swap(atoms[1], atoms[j]);
    for(size_t k = 2; k < 6; ++k){
        std::swap(atoms[2], atoms[k]);
    for(size_t l = 3; l < 6; ++l){
        std::swap(atoms[3], atoms[l]);
    for(size_t m = 4; m < 6; ++m){
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
        akv->assign("molecule", mol_t.pointer());
        Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

        // Construct via gbs
        for(size_t z = 1; z <= 6; ++z){
            std::cout << "i = " << i << " z = " << z << std::endl;
            Ref<TA::TiledBasisSet> tbs_from_gbs = new TA::TiledBasisSet(bs,z);
            BOOST_REQUIRE(!tbs_from_gbs.null());
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

BOOST_AUTO_TEST_CASE( tiledbasisset_kval_constructor_test ){

    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    std::vector<sc::Atom> atoms = mol->atoms();

    for(size_t i = 0; i < 6; ++i){
        std::swap(atoms[0], atoms[i]);
    for(size_t j = 1; j < 6; ++j){
        std::swap(atoms[1], atoms[j]);
    for(size_t k = 2; k < 6; ++k){
        std::swap(atoms[2], atoms[k]);
    for(size_t l = 3; l < 6; ++l){
        std::swap(atoms[3], atoms[l]);
    for(size_t m = 4; m < 6; ++m){
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
        akv->assign("molecule", mol_t.pointer());
        Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

        // Construct via kval
        for(size_t i = 1; i <= 6; ++i){
            akv->assign("ntiles", int(i));
            Ref<TA::TiledBasisSet> tbs_from_kval = new TA::TiledBasisSet(
                                                    Ref<KeyVal>(akv));
            BOOST_REQUIRE(!tbs_from_kval.null());
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

BOOST_AUTO_TEST_CASE( tiledbasisset_trange_test ){

    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    std::vector<sc::Atom> atoms = mol->atoms();

    for(size_t i = 0; i < 6; ++i){
        std::swap(atoms[0], atoms[i]);
    for(size_t j = 1; j < 6; ++j){
        std::swap(atoms[1], atoms[j]);
    for(size_t k = 2; k < 6; ++k){
        std::swap(atoms[2], atoms[k]);
    for(size_t l = 3; l < 6; ++l){
        std::swap(atoms[3], atoms[l]);
    for(size_t m = 4; m < 6; ++m){
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
        akv->assign("molecule", mol_t.pointer());
        Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

        // Construct via gbs
        for(size_t i = 1; i <= 6; ++i){
            Ref<TA::TiledBasisSet> tbs_from_gbs = new TA::TiledBasisSet(bs,i);
            BOOST_REQUIRE(!tbs_from_gbs.null());
            ::TiledArray::TiledRange1 trange1 = tbs_from_gbs->trange1();
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

BOOST_AUTO_TEST_CASE( tiledbasisset_trange1_values_test ) {

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

     // Test Single tile case
    Ref<TA::TiledBasisSet> tbs1 = new TA::TiledBasisSet(bs,1);
    double first = tbs1->trange1().tile(0).first;
    double second = tbs1->trange1().tile(0).second;
    BOOST_CHECK(first == 0);
    BOOST_CHECK(second == 43);

     // Test 2 tile case
    Ref<TA::TiledBasisSet> tbs2 = new TA::TiledBasisSet(bs,2);
    double t_21_f = tbs2->trange1().tile(0).first;
    double t_21_s = tbs2->trange1().tile(0).second;
    double t_22_f = tbs2->trange1().tile(1).first;
    double t_22_s = tbs2->trange1().tile(1).second;
    BOOST_CHECK(t_21_f == 0);
    BOOST_CHECK(t_21_s == 42);
    BOOST_CHECK(t_22_f == 42);
    BOOST_CHECK(t_22_s == 43);

     // Test 3 tile case
    Ref<TA::TiledBasisSet> tbs3 = new TA::TiledBasisSet(bs,3);
    double t1_first = tbs3->trange1().tile(0).first;
    double t1_second = tbs3->trange1().tile(0).second;
    double t2_first = tbs3->trange1().tile(1).first;
    double t2_second = tbs3->trange1().tile(1).second;
    double t3_first = tbs3->trange1().tile(2).first;
    double t3_second = tbs3->trange1().tile(2).second;
    BOOST_CHECK(t1_first == 0);
    BOOST_CHECK(t1_second == 24);
    BOOST_CHECK(t2_first == 24);
    BOOST_CHECK(t2_second == 25);
    BOOST_CHECK(t3_first == 25);
    BOOST_CHECK(t3_second == 43);

     // Test 4 tile case
    Ref<TA::TiledBasisSet> tbs4 = new TA::TiledBasisSet(bs,4);
    double t_41_f = tbs4->trange1().tile(0).first;
    double t_41_s = tbs4->trange1().tile(0).second;
    double t_42_f = tbs4->trange1().tile(1).first;
    double t_42_s = tbs4->trange1().tile(1).second;
    double t_43_f = tbs4->trange1().tile(2).first;
    double t_43_s = tbs4->trange1().tile(2).second;
    double t_44_f = tbs4->trange1().tile(3).first;
    double t_44_s = tbs4->trange1().tile(3).second;
    BOOST_CHECK(t_41_f == 0);
    BOOST_CHECK(t_41_s == 15);
    BOOST_CHECK(t_42_f == 15);
    BOOST_CHECK(t_42_s == 16);
    BOOST_CHECK(t_43_f == 16);
    BOOST_CHECK(t_43_s == 34);
    BOOST_CHECK(t_44_f == 34);
    BOOST_CHECK(t_44_s == 43);

     // Test 5 tile case
    Ref<TA::TiledBasisSet> tbs5 = new TA::TiledBasisSet(bs,5);
    double t_51_f = tbs5->trange1().tile(0).first;
    double t_51_s = tbs5->trange1().tile(0).second;
    double t_52_f = tbs5->trange1().tile(1).first;
    double t_52_s = tbs5->trange1().tile(1).second;
    double t_53_f = tbs5->trange1().tile(2).first;
    double t_53_s = tbs5->trange1().tile(2).second;
    double t_54_f = tbs5->trange1().tile(3).first;
    double t_54_s = tbs5->trange1().tile(3).second;
    double t_55_f = tbs5->trange1().tile(4).first;
    double t_55_s = tbs5->trange1().tile(4).second;
    BOOST_CHECK(t_51_f == 0);
    BOOST_CHECK(t_51_s == 10);
    BOOST_CHECK(t_52_f == 10);
    BOOST_CHECK(t_52_s == 11);
    BOOST_CHECK(t_53_f == 11);
    BOOST_CHECK(t_53_s == 29);
    BOOST_CHECK(t_54_f == 29);
    BOOST_CHECK(t_54_s == 38);
    BOOST_CHECK(t_55_f == 38);
    BOOST_CHECK(t_55_s == 43);

     // Test 6 tile case
    Ref<TA::TiledBasisSet> tbs6 = new TA::TiledBasisSet(bs,6);
    double t_61_f = tbs6->trange1().tile(0).first;
    double t_61_s = tbs6->trange1().tile(0).second;
    double t_62_f = tbs6->trange1().tile(1).first;
    double t_62_s = tbs6->trange1().tile(1).second;
    double t_63_f = tbs6->trange1().tile(2).first;
    double t_63_s = tbs6->trange1().tile(2).second;
    double t_64_f = tbs6->trange1().tile(3).first;
    double t_64_s = tbs6->trange1().tile(3).second;
    double t_65_f = tbs6->trange1().tile(4).first;
    double t_65_s = tbs6->trange1().tile(4).second;
    double t_66_f = tbs6->trange1().tile(5).first;
    double t_66_s = tbs6->trange1().tile(5).second;
    BOOST_CHECK(t_61_f == 0);
    BOOST_CHECK(t_61_s == 1);
    BOOST_CHECK(t_62_f == 1);
    BOOST_CHECK(t_62_s == 2);
    BOOST_CHECK(t_63_f == 2);
    BOOST_CHECK(t_63_s == 20);
    BOOST_CHECK(t_64_f == 20);
    BOOST_CHECK(t_64_s == 29);
    BOOST_CHECK(t_65_f == 29);
    BOOST_CHECK(t_65_s == 34);
    BOOST_CHECK(t_66_f == 34);
    BOOST_CHECK(t_66_s == 43);

}

BOOST_AUTO_TEST_CASE( tiledbasisset_deterministic_tiliing_test ){

    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 0.0);
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 2.8, 0.0, 1.8);
    mol->add_atom(18, 0.0, 2.8, 1.8);
    mol->add_atom(36, 0.0, -2.8, 1.8);
    mol->add_atom(18, -2.8, 0.0, 1.8);

    std::vector<sc::Atom> atoms = mol->atoms();



    for(size_t i = 0; i < 6; ++i){
        std::swap(atoms[0], atoms[i]);
    for(size_t j = 1; j < 6; ++j){
        std::swap(atoms[1], atoms[j]);
    for(size_t k = 2; k < 6; ++k){
        std::swap(atoms[2], atoms[k]);
    for(size_t l = 3; l < 6; ++l){
        std::swap(atoms[3], atoms[l]);
    for(size_t m = 4; m < 6; ++m){
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
        akv->assign("molecule", mol_t.pointer());
        Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

        // Test 3 tile case and make sure tiling is the same.
        Ref<TA::TiledBasisSet> tbs3 = new TA::TiledBasisSet(bs,3);
        double t1_first = tbs3->trange1().tile(0).first;
        double t1_second = tbs3->trange1().tile(0).second;
        double t2_first = tbs3->trange1().tile(1).first;
        double t2_second = tbs3->trange1().tile(1).second;
        double t3_first = tbs3->trange1().tile(2).first;
        double t3_second = tbs3->trange1().tile(2).second;
        BOOST_CHECK_EQUAL(t1_first, 0);
        BOOST_CHECK_EQUAL(t1_second, 24);
        BOOST_CHECK_EQUAL(t2_first, 24);
        BOOST_CHECK_EQUAL(t2_second, 25);
        BOOST_CHECK_EQUAL(t3_first, 25);
        BOOST_CHECK_EQUAL(t3_second, 43);

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

BOOST_AUTO_TEST_CASE( tiledbasisset_deterministic_overlap_test ){

    Ref<mpqc::World> world = new mpqc::World;

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

    Ref<TA::TiledBasisSet> tbs = new TA::TiledBasisSet(bs,3);

    Ref<Integral> integral = Integral::get_default_integral()->clone();
    integral->set_basis(tbs);

    IntegralEnginePool<Ref<OneBodyInt> > overlap_pool(integral->overlap());

    ::TiledArray::Array<double,2> S = Integrals(*world->madworld(),
                                                overlap_pool, tbs);

    for(size_t i = 0; i < 6; ++i){
        std::swap(atoms[0], atoms[i]);
    for(size_t j = 1; j < 6; ++j){
        std::swap(atoms[1], atoms[j]);
    for(size_t k = 2; k < 6; ++k){
        std::swap(atoms[2], atoms[k]);
    for(size_t l = 3; l < 6; ++l){
        std::swap(atoms[3], atoms[l]);
    for(size_t m = 4; m < 6; ++m){
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
        akv->assign("molecule", mol_t.pointer());
        Ref<GaussianBasisSet> bs = new GaussianBasisSet(akv);

        // Test 3 tile case and make sure tiling is the same.
        Ref<TA::TiledBasisSet> tbs3 = new TA::TiledBasisSet(bs,3);

        integral->set_basis(tbs3);
        IntegralEnginePool<Ref<OneBodyInt> > overlap_pool2(integral->overlap());


        ::TiledArray::Array<double,2> S_perm = Integrals(*world->madworld(),
                                                    overlap_pool2, tbs3);
        ::TiledArray::Array<double,2> diff = S("i,j") - S_perm("i,j");
        double infnorm =  ::TiledArray::expressions::norminf(diff("i,j"));
        world->madworld()->gop.fence();

        double zero = 0.0;

        BOOST_CHECK_CLOSE(infnorm, zero, 1e-15);

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











