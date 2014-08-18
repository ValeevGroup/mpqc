//
// taclhf_test.cpp
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
#include <chemistry/qc/libint2/libint2.h>
#include <iostream>
#include <chemistry/qc/scf/taclhf.hpp>
#include <chemistry/qc/scf/cldfgengine.hpp>
#define BOOST_TEST_MODULE test_taclhf
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;
using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

struct Mock_CLHF : public CLHF{
  Mock_CLHF(const sc::Ref<sc::KeyVal> &kval) : CLHF(kval) {}
  ~Mock_CLHF() {}

  CLHF::TAMatrix data_fock(){return scf_ao_fock_();}
  CLHF::TAMatrix data_Gmat(){return G("i,j");}
  double energy(){return iter_energy() +
                         molecule()->nuclear_repulsion_energy();}
  };

struct MADConfig {
    MADConfig() {
      int argc = boost::unit_test::framework::master_test_suite().argc;
      char** argv = boost::unit_test::framework::master_test_suite().argv;
      ExEnv::init(argc, argv);
      mpqc::MADNESSRuntime::initialize();
    }
    ~MADConfig() {
      mpqc::MADNESSRuntime::finalize();
    }
};

BOOST_GLOBAL_FIXTURE( MADConfig );

BOOST_AUTO_TEST_CASE( construct_clhf_programmatically ){

    // Make a molecule H2
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(1,0,0,0);
    mol->add_atom(1,0,0,1.4);

    // Make world
    Ref<World> world = new World;

    // Make keyval
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "cc-pVDZ");
    akv->assign("ntiles", 2);
    akv->assign("molecule", mol.pointer());
    Ref<TiledBasisSet> basis =
                    new TiledBasisSet(Ref<KeyVal>(akv));
    akv->assign("basis", basis.pointer());

    akv->assign("name", "cc-pVDZ-JK/FIT");
    Ref<TiledBasisSet> dfbasis =
                    new TiledBasisSet(Ref<KeyVal>(akv));
    akv->assign("dfbasis", dfbasis.pointer());
    akv->assign("world", world.pointer());

    Ref<Integral> int_fac = new IntegralLibint2;
    akv->assign("integrals", int_fac.pointer());

    Ref<GEngineBase> Geng = new ClDFGEngine(static_cast<Ref<KeyVal> >(akv));
    akv->assign("GEngine", Geng.pointer());

    //Construct object
    Ref<Mock_CLHF> thf = new Mock_CLHF(akv);
    thf->print();
    std::cout << "\tEnergy for initial Fock = " << thf->energy() << std::endl;
}




