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

#include <chemistry/qc/scf/cldffgfactory.hpp>
#define BOOST_TEST_MODULE test_gfactory
#include <boost/test/included/unit_test.hpp>

using namespace mpqc;
using namespace mpqc::TA;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(test_dfgfactory){
  sc::ExEnv::init(argc,argv);
  MADNESSRuntime::initialize();

  sc::Ref<World> world = new World();

  sc::Ref<sc::Integral> ints = sc::Integral::get_default_integral();

  sc::Ref<sc::Molecule>  mol = new sc::Molecule();
  mol->add_atom(1,0,0,0);
  mol->add_atom(1,0,0,1.4);

  sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
  akv->assign("name", "STO-3G");
  akv->assign("molecule", mol.pointer());
  akv->assign("world", world.pointer());
  akv->assign("ntiles", 1);

  sc::Ref<TiledBasisSet> tbs = new TiledBasisSet(akv);
  akv-assign("name", "3-21G");
  sc::Ref<TiledBasisSet> dftbs = new TiledBasisSet(akv);


  TiledArray::TiledRange trange(tbs->trange1().begin(), tbs->trange1().end());

  ClDfGFactory::TAMatrix dens(world->madworld(), trange);

  ClDfGFactory G(ints, tbs, dftbs, dens, world);

  std::cout << G("i,j") << std::endl;

  MADNESSRuntime::finalize();
}

