//
// cldfgfactory_test.cpp
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

#include <chemistry/qc/scf/cldfgengine.hpp>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <util/madness/world.h>
#include <util/madness/init.h>
#include <chemistry/qc/libint2/libint2.h>
#include <iostream>
#include <util/madness/init.h>
#define BOOST_TEST_MODULE test_gfactory
#include <boost/test/included/unit_test.hpp>

using namespace mpqc;
using namespace sc;
using namespace mpqc::TA;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(cldfgengine_test){
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  sc::ExEnv::init(argc, argv);
  mpqc::MADNESSRuntime::initialize();

  sc::Ref<World> world = new World();

  sc::Ref<sc::Molecule>  mol = new sc::Molecule();
  mol->add_atom(1,0,0,0);
  mol->add_atom(3,0,0,1.4);

  sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
  akv->assign("name", "cc-pVDZ");
  akv->assign("molecule", mol.pointer());
  akv->assign("world", world.pointer());
  akv->assign("ntiles", 2);

  sc::Ref<TiledBasisSet> tbs = new TiledBasisSet(sc::Ref<sc::KeyVal>(akv));
  sc::Ref<sc::AssignedKeyVal> akv2 = new sc::AssignedKeyVal;
  akv2->assign("name", "cc-pVDZ-RI");
  akv2->assign("molecule", mol.pointer());
  akv2->assign("world", world.pointer());
  akv2->assign("ntiles", 2);
  sc::Ref<TiledBasisSet> dftbs = new TiledBasisSet(sc::Ref<sc::KeyVal>(akv2));

  // H2 only has 1 occupied orbital
  std::array<std::size_t,2> short_tile{{0,2}};
  TiledArray::TiledRange1 short_t1(begin(short_tile), end(short_tile));

  std::array<TiledArray::TiledRange1, 2> blocking{{ tbs->trange1(), short_t1}};

  TiledArray::TiledRange trange(blocking.begin(), blocking.end());

  ClDFGEngine::TAMatrix C(*world->madworld(), trange);
  C.set_all_local(0.77);
  std::cout << "C = \n" << C << std::endl;

  ClDFGEngine::TAMatrix D; D("i,j") = C("i,k") * C("j,k");
  std::cout << "D = \n" << D << std::endl;

  sc::Ref<sc::IntegralLibint2> ints = new sc::IntegralLibint2(sc::Ref<sc::KeyVal>(akv));

  akv->assign("basis", tbs.pointer());
  akv->assign("dfbasis", dftbs.pointer());
  akv->assign("integrals", ints.pointer());

  ClDFGEngine G(static_cast<sc::Ref<sc::KeyVal> >(akv));
  if(!G.coefficients_set()){ G.set_coefficients({&C}); }
  world->madworld()->gop.fence();
  std::cout << "Coefficients have been set." << std::endl;

  ClDFGEngine GD(static_cast<sc::Ref<sc::KeyVal> >(akv));
  if(!GD.densities_set()){ GD.set_densities({&D}); }
  world->madworld()->gop.fence();

  ClDFGEngine::TAMatrix Gmat = G("i,j");
  world->madworld()->gop.fence();
  std::cout << "G = \n" << Gmat << std::endl;

  ClDFGEngine::TAMatrix GDmat = GD("i,j");
  world->madworld()->gop.fence();
  std::cout << "GD = \n" << GDmat << std::endl;
}
